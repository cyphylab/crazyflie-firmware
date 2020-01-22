#include "DataDriven.h"

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

#define STATE_SIZE (3)
#define BUFF_SIZE (5)
#define TS (0.004f)
#define TS2 ((TS) * (TS))
#define PWM_THRESHOLD (60000.0f)
#define SAFETY_HEIGHT (2.0f)



// Declare functions to use:

void DDcontroller_Step(float y, float stamp);
void DDcontroller_Init(void);
static void insert_newmeas_batch(float y, float stamp, int k);
static void insert_newmeas_circ(float y, float stamp);
static void update_O(float t[BUFF_SIZE]);
static void finalize_data_circ();
static void finalize_data_batch();
static void DDestimator_State();
static void init_O();
static float DDestimator_AlphaBeta(float TotalTime, float ctrl_dd, float ctrl_ddd);
static float compute_ctrl(float alpha, float beta);
static void DDcontroller_ScaleSet_Control(float u);
void DDestimator_Set_Ready();
bool DDestimator_Check_NewMeasurement();
// ===================================
// MEMORY BUFFERS

// A matrix
float A[STATE_SIZE * STATE_SIZE] = {
	1.0,	-(TS), (0.5f * TS * TS),
	0.0,	1.0, -(TS),
	0.0, 	0.0,	1.0,
};


// On line i: [1, -sum(ts(k)), 1/2 * (sum(ts(k))^2]
float O[BUFF_SIZE * STATE_SIZE]; /*{
								   1, 0, 0,
								   1, -TS, 0.5f * TS2,
								   1, -2.0f * TS, 2.0f * TS2,
								   1, -3.0f * TS, 9.0f/2.0f * TS2,
								   1, -4.0f * TS, 8.0f * TS2};
								   */

// C matrix
float C[STATE_SIZE] = {1, 0, 0};

// Pseudo inverse
float O_inv[STATE_SIZE * BUFF_SIZE];


//
// Measurement Buffer
static float Ybuff[BUFF_SIZE];
static float X[STATE_SIZE];
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

arm_matrix_instance_f32 Am = {STATE_SIZE, STATE_SIZE, A};
arm_matrix_instance_f32 Cm = {1, STATE_SIZE, C};
arm_matrix_instance_f32 Om = {BUFF_SIZE, STATE_SIZE, O};
arm_matrix_instance_f32 O_invm = {STATE_SIZE, BUFF_SIZE, O_inv};

arm_matrix_instance_f32 TempNyNxm = {BUFF_SIZE, STATE_SIZE, TempNyNx};
arm_matrix_instance_f32 TempNxNym = {STATE_SIZE, BUFF_SIZE, TempNxNy};
arm_matrix_instance_f32 TempNxNxm = {STATE_SIZE, STATE_SIZE, TempNxNx};
arm_matrix_instance_f32 TempNxNx2m = {STATE_SIZE, STATE_SIZE, TempNxNx2};

arm_matrix_instance_f32 Ybuffm = {BUFF_SIZE, 1, Ybuff};
arm_matrix_instance_f32 Xm = {STATE_SIZE, 1, X};

void eval_pseudoinv(arm_matrix_instance_f32* Pseudo)
{
	arm_mat_trans_f32(&Om, &TempNxNym); // O'
	arm_mat_mult_f32(&TempNxNym, &Om, &TempNxNxm); // (O' x O)
	arm_mat_inverse_f32(&TempNxNxm, &TempNxNx2m); // (O' x O)^-1 x O' = Pseudo inverse
	arm_mat_mult_f32(&TempNxNx2m, &TempNxNym, Pseudo);
}


// ====================================
//

static bool isInit = false;
static xSemaphoreHandle mutex;

// Timestamps
static float timestamp;
static float timestamp_old;
//static uint64_t us_timestamp_old;
//static uint64_t timestamp_ctrl;

float dt_ms;
float t_s;
static uint32_t msg_counter = 0;

// ====================================
// Estimator State
static float state_z;
static float alpha_ = 0.0f;
static float beta_ = 1.0f;
static int Method = 0;

// Estimator Parameters
static float gamma1 = 0.0f;
static float gamma2 = 0.0f;

// Control gain
static float P1 = 1.0f;
static float P2 = 1.0f;
static float Kdd[2];
static float ctrl_dd_;
static float ctrl_ddd_;
static float unscaled_ctrl_dd_;
static float unscaled_ctrl_ddd_;
static float Tracking[] = {1.0f, 0.0f};
static float excitation_Threshold = 0.8f;

// Control Placeholder
static float u = 0;
static float U = 0;

// Land Mode
static float Land = 0;

// Step Counter
static int Step = 0;

static bool ctrl_dd_active = false;
static bool ctrl_dd_start = false;
static bool updated = false;
static int WaitBetweenLearning = 0;
static int SecondStep = 0;

// =============================================
// Getters and Setters

static void DDestimator_Set_Alpha(float a)
{
	alpha_ = a;
}

static float DDestimator_Get_Alpha()
{
	return alpha_;
}

static void DDestimator_Set_Beta(float b)
{
	beta_ = b;
}
static float DDestimator_Get_Beta()
{
	return beta_;
}

void DDcontroller_Set_Control(const float u)
{
	ctrl_ddd_ = ctrl_dd_;
	ctrl_dd_ = u;
}

void DDcontroller_Set_UnscaledControl(const float u)
{
	unscaled_ctrl_ddd_ = ctrl_dd_;
	unscaled_ctrl_dd_ = u;
}
static float DDcontroller_Get_UnscaledControlOld()
{
	return unscaled_ctrl_dd_;
}

static float DDcontroller_Get_UnscaledControlOldOld()
{
	return unscaled_ctrl_ddd_;
}

// Other modules use these
float DDcontroller_Get_Control()
{
	return ctrl_dd_;
}

void DDcontroller_Set_ControllerReady(void)
{
	ctrl_dd_active = true;
}
void DDcontroller_Set_ControllerStart(void)
{
	ctrl_dd_start = true;
	DDcontroller_ScaleSet_Control(0.6f);
}
float DDestimator_Get_StateEstimate()
{
	return state_z;
}
// ===================================




/**
 * This function is triggered by the arrival of new measurements
 */
bool DDcontroller_NewMeasurement(const positionMeasurement_t *pos)
{
	msg_counter = msg_counter + 1;

	// Measure the timestamp
	timestamp = pos->t; // Time in second
	dt_ms  = (timestamp - timestamp_old)* 1e3f;
	timestamp_old = timestamp;
	state_z = pos->z;

	// Data Driven Controller

	DDcontroller_Step(state_z, timestamp);

	return true;
}

/**
 * Controller Step
*/
void DDcontroller_Step(float y, float stamp)
{
	if (!isInit) {
		DDcontroller_Init();
	}


	if (Method==0) {

		insert_newmeas_batch(y, stamp, Nmeas);
		Nmeas++;
		if (Nmeas == BUFF_SIZE) {
			Nmeas = 0;
			finalize_data_batch();

			DDestimator_State();
			if (ctrl_dd_start) {
				// Estimate Parameters
				float cdd = DDcontroller_Get_UnscaledControlOld();
				float cddd = DDcontroller_Get_UnscaledControlOldOld();
				float u = DDestimator_AlphaBeta(TotalTime,cdd,cddd);

				// Control

				DDcontroller_ScaleSet_Control(u);
			} else if (!ctrl_dd_start && X[1]>1.0f) {
				DDcontroller_Set_ControllerStart();
			} else if (!ctrl_dd_start && X[1]<-1.0f) {
				DDcontroller_Set_ControllerStart();
			}
			DDestimator_Set_Ready();
		}

	} else {

		insert_newmeas_circ(y, stamp);
		if (Nmeas >= BUFF_SIZE) {
			// State Estimation
			finalize_data_circ();
			DDestimator_State();
			if (ctrl_dd_start) {
				// Estimate Parameters
				float cdd = DDcontroller_Get_UnscaledControlOld();
				float cddd = DDcontroller_Get_UnscaledControlOldOld();
				float u = DDestimator_AlphaBeta(TotalTime,cdd,cddd);

				// Control
				DDcontroller_ScaleSet_Control(u);


			} else if (!ctrl_dd_start && X[1]>1.0f) {
				DDcontroller_Set_ControllerStart();
			} else if (!ctrl_dd_start && X[1]<-1.0f) {
				DDcontroller_Set_ControllerStart();
			}
			DDestimator_Set_Ready();
		} else {
			Nmeas++;
		}
	}
}


/**
 * Initialize controller
 */
void DDcontroller_Init(void)
{
	if (isInit)  {
		return;
	}

	Kdd[0] = -P1 * P2;
	Kdd[1] = P1 + P2;
	Kdd[2] = 0.0f;

	DEBUG_PRINT("DD Controller Gain: [%.3f, %.3f] \n", (double)Kdd[0], (double)Kdd[1]);
	init_O();

	eval_pseudoinv(&O_invm);

	mutex = xSemaphoreCreateMutex();
	alpha_ = 0.0f;
	beta_ = 18.7291f;
	// Initialize the DD Library
	DataDriven_initialize();

	isInit = true;
	return;
}


/**
 * Insert the k-th measurement in the buffer
 */
static void insert_newmeas_batch(float y, float stamp, int k)
{
	if (k < 0) {
		// Error
	}
	int index = (BUFF_SIZE - 1) - (k % BUFF_SIZE);
	Ybuff[index] = y;
	Tbuff[index] = stamp;
}

static void insert_newmeas_circ(float y, float stamp)
{

	for (int index = 1; index < BUFF_SIZE; index++) {
		Ybuff[BUFF_SIZE-index] = Ybuff[BUFF_SIZE-index-1];
		Tbuff[BUFF_SIZE-index] = Tbuff[BUFF_SIZE-index-1];
	}
	Ybuff[0] = y;
	Tbuff[0] = stamp;
}

static void update_O(float t[BUFF_SIZE])
{
	int i;
	for (i = 0; i < BUFF_SIZE; i++) {
		O[(i*STATE_SIZE) + 1] = -(t[i]);
		O[(i*STATE_SIZE) + 2] = 0.5f * (t[i] * t[i]);
	}
}

static void finalize_data_circ()
{

	// Finalize the DT vector, computing the differences
	// [0, dt1, dt1 + dt2, ...]
	DTbuff[0] = Tbuff[0];
	for (int i = 0; i < BUFF_SIZE; i++) {
		DTbuff[i] = Tbuff[0] - Tbuff[i];
	}
	TotalTime = DTbuff[1];
	// Update the Observability matrix
	update_O(DTbuff);
	// Update the pseduoinverse matrix
	eval_pseudoinv(&O_invm);
}

static void finalize_data_batch()
{

	// Finalize the DT vector, computing the differences
	// [0, dt1, dt1 + dt2, ...]
	DTbuff[0] = Tbuff[0];
	for (int i = 0; i < BUFF_SIZE; i++) {
		DTbuff[i] = Tbuff[0] - Tbuff[i];
	}
	//TotalTime = DTbuff[BUFF_SIZE-1];
	TotalTime = DTbuff[4];


	// Update the Observability matrix
	update_O(DTbuff);

	// Update the pseduoinverse matrix
	eval_pseudoinv(&O_invm);
}

void init_O()
{
	int i;
	for (i = 0; i < BUFF_SIZE; i++) {
		O[(i*STATE_SIZE) + 0] = 1;
		O[(i*STATE_SIZE) + 1] = -(i * TS);
		O[(i*STATE_SIZE) + 2] = 0.5f * (i * i * TS2);
	}
}

/**
 * Estimate state
 */
static void DDestimator_State()
{
	// Save the old state before the update
	for (int i = 0; i < 3; i++) {
		X_old_old[i] = X_old[i];
		X_old[i] = X[i];
	}

	arm_mat_mult_f32(&O_invm, &Ybuffm, &Xm);
	return;
}


/**
 * Estimate alpha and beta, compute u
 */
static float DDestimator_AlphaBeta(float TotalTime, float ctrl_dd, float ctrl_ddd)
{

	// Get the current value from the global variables
	float alpha = DDestimator_Get_Alpha();
	float beta = DDestimator_Get_Beta();

	// Local variables
	float alpha_new = 0.0f;
	float beta_new = 0.0f;

	float diff_1 = (X[1] - X_old[1]);
	float udiff;
	float u = 0.0f;

	if(SecondStep) {
		float GainObserverAlpha2 = - ctrl_ddd / (TotalTime * (ctrl_dd - ctrl_ddd));
		float GainObserverBeta2 = 1.0f / (TotalTime * (ctrl_dd * ctrl_dd - ctrl_ddd * ctrl_dd));

		alpha_new = alpha + GainObserverAlpha2 * (diff_1 - TotalTime * (alpha + beta * ctrl_dd));
		beta_new = beta + GainObserverBeta2 * ctrl_dd * (diff_1 - TotalTime * (alpha + beta * ctrl_dd));
		DEBUG_PRINT("Recompute0 [ %.6f, %.6f]\n", (double)alpha_new, (double)beta_new);
		if (beta_new < 3.0f) {
			beta_new = 3.0f;
			alpha_new = -2.0f;
			DEBUG_PRINT("Recompute1 [ %.6f, %.6f]\n", (double)alpha_new, (double)beta_new);
		} else if (beta_new > 25.0f) {
			beta_new = 25.0f;
			alpha_new = -15.0f;
			DEBUG_PRINT("Recompute2 [ %.6f, %.6f]\n", (double)alpha_new, (double)beta_new);
		} else {
			DEBUG_PRINT("Recompute3 [ %.6f, %.6f]\n", (double)alpha_new, (double)beta_new);
		}
		WaitBetweenLearning = 0;
		SecondStep = 0;

		u = compute_ctrl(alpha_new, beta_new);

	} else {
		float GainObserverBeta1 = 0.0f;
		float GainObserverAlpha1 = - (TotalTime * GainObserverBeta1 * ctrl_dd * ctrl_dd - 1.0f) / TotalTime;


		alpha_new = alpha + GainObserverAlpha1 * (diff_1 - TotalTime * (alpha + beta * ctrl_dd));
		beta_new = beta + GainObserverBeta1 * ctrl_dd * (diff_1 - TotalTime * (alpha + beta * ctrl_dd));

		u = compute_ctrl(alpha_new, beta_new);

		udiff = u - ctrl_dd;

		if (udiff < 0.0f) {
			udiff = -udiff;
		}

		if (!(udiff >= excitation_Threshold && (WaitBetweenLearning >= 10 || Step < 4) && u != 0.0f)) {

			alpha_new = alpha + gamma1 * (diff_1 - TotalTime * (alpha + beta * ctrl_dd));
			beta_new = beta + gamma2 * ctrl_dd * (diff_1 - TotalTime * (alpha + beta * ctrl_dd));


			if (beta_new < 3.0f) {
				beta_new = 3.0f;
				alpha_new = -2.0f;
				/*DEBUG_PRINT("Input 14 [ %.3f, %.3f]\n", (double)alpha_new, (double)beta_new);*/
			}

			if (beta_new > 25.0f) {
				beta_new = 25.0f;
				alpha_new = -15.0f;
				/*DEBUG_PRINT("Input 15 [ %.3f, %.3f]\n", (double)alpha_new, (double)beta_new);*/
			}

			WaitBetweenLearning ++;
			u = compute_ctrl(alpha_new, beta_new);

		} else {
			SecondStep = 1;
		}
	}


	Step++;
	if (Step == 200) {
		DEBUG_PRINT("Alpha Beta u [ %.3f, %.3f, %.3f]\n",
		            (double)alpha_new, (double)beta_new,(double)u);
		Step=4;
	}

// Update the parameters
	DDestimator_Set_Alpha(alpha_new);
	DDestimator_Set_Beta(beta_new);

	return u;
}

/**
 * Compute control value
 */
static float compute_ctrl(float alpha, float beta)
{
	float u_fb = 0;
	static float u_p;
	static float u_d;

	// Failsafe for Fly Away
	if (X[0]>SAFETY_HEIGHT) {
		Land=1;
	}
	// Update the control gain
	Kdd[0] = - P1 * P2;
	Kdd[1] = P1 + P2;

	u_p = Kdd[0] * (X[0] - Tracking[0]);
	u_d = Kdd[1] * (X[1] - Tracking[1]);

	u_fb = u_p + u_d;
	u = (1.0f / beta) * (-alpha + u_fb);
	if (u < 0.0f) {
		u = 0.0f;
	}
	if (u > 1.0f) {
		u = 1.0f;
	}
	return u;
}

/**
 * Scale and Push control value
*/

static void DDcontroller_ScaleSet_Control(float u)
{
	DDcontroller_Set_UnscaledControl(u);
	U = u * PWM_THRESHOLD;
	if (!Land) {
		DDcontroller_Set_Control(U);
	} else {
		DDcontroller_Set_Control(0.0);
	}
}

/**
 * Flag control value ready
 */
void DDestimator_Set_Ready()
{
	updated = true;
	return;
}

/**
 * Check whether the estimator is ready. In that case
 * reset the flag and return 'true'. Otherwise, return 'false'.
 */
bool DDestimator_Check_NewMeasurement()
{
	bool out;

	//xSemaphoreTake(mutex, portMAX_DELAY);
	out = updated;
	if (out)
		updated = false; // Reset the flag
	//xSemaphoreGive(mutex);

	return out;
}

bool DDestimator_Check_Init(void)
{
	return isInit;
}




/**
 * Feeds direct messages for testing. Read requirements below.
 */
/**
 * Feed the DD controller with the state estimate.
 * (It will be necessary to disable the callback from the sensor to avoid
 * overlapping the estimation (and race conditions on the variables)
 *
 * If this is a feature that would be embedded in the final version we
 * can work on a logic for all these workarounds...
 */

/*void DDestimator_Feed_State(float z, float zd, uint64_t us_timestamp)
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
	float cdd = DDcontroller_Get_UnscaledControlOld();
	float cddd = DDcontroller_Get_UnscaledControlOldOld();
	float u = DDestimator_AlphaBeta(TotalTime,cdd,cddd);

	if (ctrl_dd_active && Step > 1) {
		set_Control(u);
	}

	DDestimator_Set_Ready();

}*/

// Logging variables
//
LOG_GROUP_START(estimator_dd)
LOG_ADD(LOG_FLOAT, est_x, &X[0])
LOG_ADD(LOG_FLOAT, est_xd, &X[1])
LOG_ADD(LOG_FLOAT, est_alpha, &alpha_)
LOG_ADD(LOG_FLOAT, est_beta, &beta_)
LOG_ADD(LOG_FLOAT, sens_dt_ms, &dt_ms)
LOG_GROUP_STOP(estimator_dd)

PARAM_GROUP_START(controller_dd)
PARAM_ADD(PARAM_FLOAT, ctrl_ddP1, &P1)
PARAM_ADD(PARAM_FLOAT, ctrl_ddP2, &P2)
PARAM_ADD(PARAM_FLOAT, ctrl_ddg1, &gamma1)
PARAM_ADD(PARAM_FLOAT, ctrl_ddg2, &gamma2)
PARAM_ADD(PARAM_FLOAT, ctrl_ddA, &alpha_)
PARAM_ADD(PARAM_FLOAT, ctrl_ddB, &beta_)
PARAM_ADD(PARAM_FLOAT, ctrl_ddTr, &Tracking[0])
PARAM_ADD(PARAM_FLOAT, ctrl_ddLd, &Land)
PARAM_ADD(PARAM_UINT8, ctrl_0, &ctrl_dd_start)
PARAM_ADD(PARAM_FLOAT, ctrl_thr, &excitation_Threshold)
PARAM_ADD(PARAM_FLOAT, ctrl_Mtd, &Method)
PARAM_GROUP_STOP(controller_dd)
