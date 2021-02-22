#include "lfd.h"

#include "FreeRTOS.h"
#include "queue.h"
//#include "task.h"
#include "semphr.h"
#include "debug.h"

#include "stabilizer.h"

#include "sensors.h"

#include "log.h"
#include "param.h"

#include "arm_math.h"

#include "math.h"

#include "lfd_expertk.h"


#define STATE_SIZE (3)
#define BUFF_SIZE_U (3)
#define BUFF_SIZE (5)
#define TS (0.004f)
#define TS2 ((TS) * (TS))
#define PWM_THRESHOLD (60000.0f)
#define SAFETY_HEIGHT (2.0f)

// Declare functions to use:

void LfDcontroller_Step(float y, float stamp);
void LfDcontroller_Init(void);
static void insert_newmeas_circ(float y, float stamp);
static void update_O(float t[BUFF_SIZE]);
static void finalize_data_circ();
static void LfDestimator_State();
static void init_O();
static float compute_ctrl(int K_column);
static void LfDcontroller_ScaleSet_Control(float u);
void LfDestimator_Set_Ready();
bool LfDestimator_Check_NewMeasurement();

// ===================================
// MEMORY BUFFERS

// On line i: [1, -sum(ts(k)), 1/2 * (sum(ts(k))^2]
float O[BUFF_SIZE * STATE_SIZE]; /*{
								   1, 0, 0,
								   1, -TS, 0.5f * TS2,
								   1, -2.0f * TS, 2.0f * TS2,
								   1, -3.0f * TS, 9.0f/2.0f * TS2,
								   1, -4.0f * TS, 8.0f * TS2};
								   */


// Pseudo inverse
float O_inv[STATE_SIZE * BUFF_SIZE];

//
// Measurement Buffer
static float Ybuff[BUFF_SIZE];
static float X[STATE_SIZE];
static float Z[STATE_SIZE];
static float X_old[3];
static float X_old_old[3];

static float Tbuff[BUFF_SIZE];
static float DTbuff[BUFF_SIZE];
static float TotalTime;

// Temp Buffers for the evaluation of the pseudoinverse
float TempNyNx[BUFF_SIZE * STATE_SIZE];
float TempNxNy[STATE_SIZE * BUFF_SIZE];
float TempNxNx[STATE_SIZE * STATE_SIZE];
float TempNxNx2[STATE_SIZE * STATE_SIZE];

// Filter Data
static int Nmeas = 0;

arm_matrix_instance_f32 Om = {BUFF_SIZE, STATE_SIZE, O};
arm_matrix_instance_f32 O_invm = {STATE_SIZE, BUFF_SIZE, O_inv};

arm_matrix_instance_f32 TempNyNxm = {BUFF_SIZE, STATE_SIZE, TempNyNx};
arm_matrix_instance_f32 TempNxNym = {STATE_SIZE, BUFF_SIZE, TempNxNy};
arm_matrix_instance_f32 TempNxNxm = {STATE_SIZE, STATE_SIZE, TempNxNx};
arm_matrix_instance_f32 TempNxNx2m = {STATE_SIZE, STATE_SIZE, TempNxNx2};

arm_matrix_instance_f32 Ybuffm = {BUFF_SIZE, 1, Ybuff};
arm_matrix_instance_f32 Xm = {STATE_SIZE, 1, X};

void eval_pseudoinv(arm_matrix_instance_f32 *Pseudo)
{
    arm_mat_trans_f32(&Om, &TempNxNym);                // O'
    arm_mat_mult_f32(&TempNxNym, &Om, &TempNxNxm);     // (O' x O)
    arm_mat_inverse_f32(&TempNxNxm, &TempNxNx2m);      // (O' x O)^-1 = Pseudo inverse
    arm_mat_mult_f32(&TempNxNx2m, &TempNxNym, Pseudo); // (O' x O)^-1 x O'
}

// ====================================
//

static bool isInit = false;
static xSemaphoreHandle mutex;

// Timestamps
static float timestamp;
static float timestamp_old;
static int exp_traj_l = 567;
static float exp_traj_timestep = 0.003f;
//static uint64_t us_timestamp_old;
//static uint64_t timestamp_ctrl;

float dt_ms;
float t_s;
float Time0;
static uint32_t msg_counter = 0;

// ====================================
// Estimator State
static float state_z;
static float alpha_ = -9.81f;
static float beta_ = 18.7291f;
static float ctrl_lfd_;
static float unscaled_ctrl_lfd_;

static float Tracking[] = {1.0f, 0.0f};

// Control Placeholder
static float u = 0;
static float Klfd[2];

// Land Mode
static float Land = 0;

static bool ctrl_lfd_active = false;
static bool ctrl_lfd_start = false;
static bool updated = false;

// =============================================
// Getters and Setters

void LfDcontroller_Set_Control(const float u)
{
    ctrl_lfd_ = u;
}

void LfDcontroller_Set_UnscaledControl(const float u)
{
    unscaled_ctrl_lfd_ = u;
}

// Other modules use these
float LfDcontroller_Get_Control()
{
    return ctrl_lfd_;
}

void LfDcontroller_Set_ControllerReady(void)
{
    ctrl_lfd_active = true;
}
void LfDcontroller_Set_ControllerStart(float stamp)
{
    ctrl_lfd_start = true;
    Time0 = stamp;
    u = compute_ctrl(0);
    LfDcontroller_ScaleSet_Control(u);
}
float LfDestimator_Get_StateEstimate()
{
    return state_z;
}
// ===================================

/**
 * This function is triggered by the arrival of new measurements
 */
bool LfDcontroller_NewMeasurement(const positionMeasurement_t *pos)
{
    msg_counter = msg_counter + 1;

    // Measure the timestamp
    timestamp = pos->t; // Time in second
    dt_ms = (timestamp - timestamp_old) * 1e3f;
    timestamp_old = timestamp;
    state_z = pos->z;

    // Data Driven Controller
    LfDcontroller_Step(state_z, timestamp);

    return true;
}

/**
 * Controller Step
*/
void LfDcontroller_Step(float y, float stamp)
{
    if (!isInit)
    {
        LfDcontroller_Init();
    }

    insert_newmeas_circ(y, stamp);
    Nmeas++;
    if (Nmeas >= BUFF_SIZE)
    {
        // State Estimation
        finalize_data_circ();
        LfDestimator_State();
        if (ctrl_lfd_start)
        {
            // Control
            float Time = stamp;
            int K_column = (int)round((Time-Time0) / exp_traj_timestep ) % exp_traj_l;
            u = compute_ctrl(K_column);
            LfDcontroller_ScaleSet_Control(u);
        }
        else{
            LfDcontroller_Set_ControllerStart(stamp);
        }
        LfDestimator_Set_Ready();
    }
}

/**
 * Initialize controller
 */
void LfDcontroller_Init(void)
{
    if (isInit)
    {
        return;
    }

    
    init_O();

    eval_pseudoinv(&O_invm);

    mutex = xSemaphoreCreateMutex();
    alpha_ = -9.81f;
    beta_ = 18.7291f;

    isInit = true;
    return;
}

/**
 * Insert the k-th measurement in the buffer
 */


static void insert_newmeas_circ(float y, float stamp)
{

    for (int index = 1; index < BUFF_SIZE; index++)
    {
        Ybuff[BUFF_SIZE - index] = Ybuff[BUFF_SIZE - index - 1];
        Tbuff[BUFF_SIZE - index] = Tbuff[BUFF_SIZE - index - 1];
    }
    Ybuff[0] = y;
    Tbuff[0] = stamp;
}


void init_O()
{
    int i;
    for (i = 0; i < BUFF_SIZE; i++)
    {
        O[(i * STATE_SIZE) + 0] = 1;
        O[(i * STATE_SIZE) + 1] = -(i * TS);
        O[(i * STATE_SIZE) + 2] = 0.5f * (i * i * TS2);
    }
}

static void update_O(float t[BUFF_SIZE])
{
    int i;
    for (i = 0; i < BUFF_SIZE; i++)
    {
        O[(i * STATE_SIZE) + 1] = -(t[i]);
        O[(i * STATE_SIZE) + 2] = 0.5f * (t[i] * t[i]);
    }
}

static void finalize_data_circ()
{

    // Finalize the DT vector, computing the differences
    // [0, dt1, dt1 + dt2, ...]
    DTbuff[0] = Tbuff[0];
    for (int i = 0; i < BUFF_SIZE; i++)
    {
        DTbuff[i] = Tbuff[0] - Tbuff[i];
    }
    TotalTime = DTbuff[1];
    // Update the Observability matrix
    update_O(DTbuff);
    // Update the pseduoinverse matrix
    eval_pseudoinv(&O_invm);
}


/**
 * Estimate state
 */
static void LfDestimator_State()
{
    // Save the old state before the update
    for (int i = 0; i < 3; i++)
    {
        X_old_old[i] = X_old[i];
        X_old[i] = X[i];
    }

    arm_mat_mult_f32(&O_invm, &Ybuffm, &Xm);
    for (int index = 0; index < (BUFF_SIZE_U - 1); index++)
    {
        Z[(BUFF_SIZE_U - index) - 1] = Z[(BUFF_SIZE_U - index) - 2];
    }
    Z[0] = X[2];
    return;
}
/**
 * Compute control value
 */
static float compute_ctrl(int K_column)
{
    float u_fb = 0;
    static float u_p;
    static float u_d;

    // Failsafe for Fly Away
    if (X[0] > SAFETY_HEIGHT)
    {
        Land = 1;
    }
    // Update the control gain

   
    expert_K(K_column, Klfd);

    u_p = Klfd[0] * (X[0] - Tracking[0]);
    u_d = Klfd[1] * (X[1] - Tracking[1]);

    u_fb = u_p + u_d;
    u = (1.0f / beta_) * (-alpha_ + u_fb);
    if (u < 0.0f)
    {
        u = 0.0f;
    }
    if (u > 1.0f)
    {
        u = 1.0f;
    }
    return u;
}

/**
 * Scale and Push control value
*/

static void LfDcontroller_ScaleSet_Control(float u)
{
    LfDcontroller_Set_UnscaledControl(u);
    /** The Following Line Changes U to fake a time varying Beta 
	

    fakeBeta = (18.7291f - Factor * (X[0] * X[0] * X[0] * X[0]));
    U = (fakeBeta / 18.7291f) * (u * PWM_THRESHOLD);
    */
    static float U = u * PWM_THRESHOLD;
    if (!Land)
    {
        LfDcontroller_Set_Control(U);
    }
    else
    {
        LfDcontroller_Set_Control(0.0);
    }
}

/**
 * Flag control value ready
 */
void LfDestimator_Set_Ready()
{
    updated = true;
    return;
}

/**
 * Check whether the estimator is ready. In that case
 * reset the flag and return 'true'. Otherwise, return 'false'.
 */
bool LfDestimator_Check_NewMeasurement()
{
    bool out;

    //xSemaphoreTake(mutex, portMAX_DELAY);
    out = updated;
    if (out)
        updated = false; // Reset the flag
    //xSemaphoreGive(mutex);

    return out;
}

bool LfDestimator_Check_Init(void)
{
    return isInit;
}

/**
 * Feeds direct messages for testing. Read requirements below.
 */
/**
 * Feed the LfD controller with the state estimate.
 * (It will be necessary to disable the callback from the sensor to avoid
 * overlapping the estimation (and race conditions on the variables)
 *
 * If this is a feature that would be embedded in the final version we
 * can work on a logic for all these workarounds...
 */

/*void LfDestimator_Feed_State(float z, float zd, uint64_t us_timestamp)
{

	// Measure the time to check whether the trigger is really periodic

	dt_ms = (float)(us_timestamp - us_timestamp_old) / 1e6f;
	us_timestamp_old = us_timestamp;
	for (int i = 0; i < 2; i++) {
		X_old[i] = X[i];
	}

	X[0] = z;
	X[1] = zd;
	// The main loop is supposed to spin at 1Khz. Clearly, the information from the Z
	// is not only provided by the camera, but takes into account the filter prediction
	// capabilities.
	TotalTime = dt_ms;
	// Estimate Parameters
	float cdd = LfDcontroller_Get_UnscaledControlOld();
	float cddd = LfDcontroller_Get_UnscaledControlOldOld();
	float u = LfDestimator_AlphaBeta(TotalTime,cdd,cddd);

	if (ctrl_lfd_active && Step > 1) {
		set_Control(u);
	}

	LfDestimator_Set_Ready();

}*/

// Logging variables
//
LOG_GROUP_START(estimator_lfd)
LOG_ADD(LOG_FLOAT, est_x, &X[0])
LOG_ADD(LOG_FLOAT, est_xd, &X[1])
LOG_ADD(LOG_FLOAT, LfDK0, &Klfd[0])
LOG_ADD(LOG_FLOAT, LfDK1, &Klfd[1])
LOG_ADD(LOG_FLOAT, fblu, &U)
LOG_GROUP_STOP(estimator_lfd)

PARAM_GROUP_START(controller_lfd)
PARAM_ADD(PARAM_FLOAT, ctrl_exptl, &exp_traj_l)
PARAM_ADD(PARAM_FLOAT, ctrl_lfdA, &alpha_)
PARAM_ADD(PARAM_FLOAT, ctrl_lfdB, &beta_)
PARAM_ADD(PARAM_FLOAT, ctrl_lfdTr, &Tracking[0])
PARAM_ADD(PARAM_FLOAT, ctrl_lfdLd, &Land)
PARAM_ADD(PARAM_UINT8, ctrl_0, &ctrl_lfd_start)
//PARAM_ADD(PARAM_FLOAT, ctrl_bfac, &Factor)
PARAM_GROUP_STOP(controller_lfd)
