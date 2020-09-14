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

#define STATE_SIZE (2)
#define BUFF_SIZE (10)
#define TS (0.004f)
#define TS2 ((TS) * (TS))
#define PWM_THRESHOLD (60000.0f)
#define SAFETY_HEIGHT (1.8f)



// Declare functions to use:

void DDcontroller_Step(float y, float stamp, float TotalTime);
void DDcontroller_Init(void);
static void DDestimator_State(float y, float TotalTime);
static float DDestimator_AlphaBeta(float TotalTime, float ctrl_dd, float ctrl_ddd);
static float compute_ctrl(float alpha, float beta);
static void DDcontroller_ScaleSet_Control(float u);
void DDestimator_Set_Ready();
bool DDestimator_Check_NewMeasurement();
// ===================================
// MEMORY BUFFERS

// A matrix
float A[4] = {
	1.0f,	(TS),
	0.0f,	1.0f,
};

// C matrix
float C[2] = {1.0f, 0.0f};

// Measurement Buffer
static float X[STATE_SIZE];
static float X_old[STATE_SIZE];

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

// Estimator Parametrs
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

// State Observer
static float OP1 = 0.70f;
static float OP2 = 0.70f;
float L[2]= {40.0f,20.0f*20.0f};

// Control Placeholder
static float u = 0;
static float U = 0;

// Land Mode
static float Land = 0;

// Step Counter
static int Step = 0;
static bool ctrl_dd_active = false;
static bool start_control = false;
static bool updated = false;
static int WaitBetweenLearning = 0;
static int SecondStep = 0;
float TotalTime;

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

void DDcontroller_Get_UnscaledControl(const float u)
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
	TotalTime = (timestamp - timestamp_old);
	dt_ms = TotalTime * 1e3f;
	timestamp_old = timestamp;
	state_z = pos->z;

	// Data Driven Controller

	DDcontroller_Step(state_z, timestamp, TotalTime);

	return true;
}

/**
 * Controller Step
*/
void DDcontroller_Step(float y, float stamp, float TotalTime)
{
	if (!isInit) {
		DDcontroller_Init();
	}

	// State Estimation

	DDestimator_State(y, TotalTime);

	// Estimate Parameters
	float cdd = DDcontroller_Get_UnscaledControlOld();
	float cddd = DDcontroller_Get_UnscaledControlOldOld();
	if (ctrl_dd_active) {
		float u = DDestimator_AlphaBeta(TotalTime,cdd,cddd);
		DDcontroller_ScaleSet_Control(u);
		/*if ((X[1]<-1) && start_control == false) {
			start_control = true;
			DDcontroller_ScaleSet_Control(u);
		} else if ((X[1]>1) && start_control == false) {
			start_control = true;
			DDcontroller_ScaleSet_Control(u);
		}
		if (start_control == true) {
			DDcontroller_ScaleSet_Control(u);
		}*/
	}
	DDestimator_Set_Ready();
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
	L[0] =  OP1 + OP2;
	L[1] = OP1 * OP2;
	DEBUG_PRINT("DD Controller Gain: [%.3f, %.3f] \n", (double)Kdd[0], (double)Kdd[1]);

	mutex = xSemaphoreCreateMutex();

	alpha_ = 0.0f;
	beta_ = 1.0f;
	// Initialize the DD Library
	isInit = true;
	return;
}

/**
 * Estimate state
 */
static void DDestimator_State(float y, float TotalTime)
{
	// Save the old state before the update
	for (int i = 0; i < 2; i++) {
		X_old[i] = X[i];
	}
	if (Step == 0) {
		X[0] = y;
		X[1] = 0.0f;
		Step++;
	}

		L[0] =  (OP1 + OP2)/TotalTime;
		L[1] = (OP1 * OP2)/(TotalTime * TotalTime);
/*	L[0] = -2.0f * 0.8f / TotalTime;
	L[1] = 0.8f * 0.8f / (TotalTime * TotalTime);*/
	/*	if (Step == 4) {
			DEBUG_PRINT("Print1 [ %.3f, %.3f, %.3f, %.3f]\n", (double)(y), (double)(TotalTime * L[0]), (double)(TotalTime * L[1]), (double)(X[1]));
		}*/
	A[1] = TotalTime;
	X[0] = A[0] * X_old[0] + A[1] * X_old[1] + TotalTime * L[0] * (y - (C[0] * X_old[0] + C[1] * X_old[1]));
	X[1] = A[2] * X_old[0] + A[3] * X_old[1] + TotalTime * L[1] * (y - (C[0] * X_old[0] + C[1] * X_old[1]));
	/*	if (Step == 4) {
			DEBUG_PRINT("Print2 [ %.3f, %.3f, %.3f, %.3f]\n", (double)(y), ((double)y - (double)(C[0] * X[0] + C[1] * X[1])), (double)(X[0]), (double)(X[1]));
		}*/
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
		} else{
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

		if (!(udiff >= excitation_Threshold && (WaitBetweenLearning >= 4 || Step < 4) && u != 0.0f)) {

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
			WaitBetweenLearning ++;
			SecondStep = 1;
		}
	}


	Step++;
	if (Step == 200) {
		DEBUG_PRINT("Alpha Beta  [ %.3f, %.3f]\n",
		            (double)alpha_new, (double)beta_new);
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
	DDcontroller_Get_UnscaledControl(u);
	U = u * PWM_THRESHOLD;
	if (!Land) {
		DDcontroller_Set_Control(U);
		if (Step == 4) {
			DEBUG_PRINT("[ %.6f, %.6f, %.6f]\n",
			            (double)alpha_, (double)beta_, (double)u);
		}
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
PARAM_ADD(PARAM_FLOAT, ctrl_0, &start_control)
PARAM_ADD(PARAM_FLOAT, ctrl_thr, &excitation_Threshold)
PARAM_ADD(PARAM_FLOAT, ctrl_OP1, &OP1)
PARAM_ADD(PARAM_FLOAT, ctrl_OP2, &OP2)
PARAM_GROUP_STOP(controller_dd)
