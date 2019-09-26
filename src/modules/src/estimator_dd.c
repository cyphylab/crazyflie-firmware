#include "estimator_dd.h"
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
#define BUFF_SIZE (10)
#define TS (0.004f)
#define TS2 ((TS) * (TS))


#define  droneMass (0.032f)

// ===================================
// MEMORY BUFFERS 

// A matrix
float A[STATE_SIZE * STATE_SIZE] = 
{
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

static float Tbuff[BUFF_SIZE];
static float DTbuff[BUFF_SIZE];
static float TotalTime;

// Temp Buffers for the evaluation of the pseudoinverse
float TempNyNx[BUFF_SIZE * STATE_SIZE];
float TempNxNy[STATE_SIZE * BUFF_SIZE];
float TempNxNx[STATE_SIZE * STATE_SIZE];
float TempNxNx2[STATE_SIZE * STATE_SIZE];




// ====================================	
//


static bool isInit = false;
static xSemaphoreHandle mutex;

// Timestamps
static float timestamp;
static float timestamp_old;
static uint64_t us_timestamp_old;
//static uint64_t timestamp_ctrl;

float dt_ms;
float t_s;
float dt_ms_cum = 0;
static uint32_t msg_counter = 0;

// ====================================
// Estimator State 
static float state_z;
static float alpha_ = 0.0f;
static float alpha_new = 0.0f;
static float beta_ = 2.8577f*1e-4f;
// /droneMass;

// Estimator Parametrs
static float gamma1 = 1.0f;
static float gamma2 = 1.0f;

// Control gain
static float P1 = 1.0f;
static float P2 = 1.0f;
static float Kdd[3];
static float ctrl_dd_;
static float ctrl_ddd_;
static float Tracking[] = {1.5f, 0.0, 0.0};

// Control Placeholder
static float u = 0;
static float U = 0;

// Land Mode
static float Land = 0;


// Step Counter
static int Step = 4.0f;
static bool ctrl_dd_active = false;
static bool estimate_least_squares = false;
static bool updated = false;

// ====================================
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

void eval_pseudoinv(arm_matrix_instance_f32* Pseudo) {
	arm_mat_trans_f32(&Om, &TempNxNym); // O'
	arm_mat_mult_f32(&TempNxNym, &Om, &TempNxNxm); // (O' x O)
	arm_mat_inverse_f32(&TempNxNxm, &TempNxNx2m); // (O' x O)^-1 x O' = Pseudo inverse
	arm_mat_mult_f32(&TempNxNx2m, &TempNxNym, Pseudo);
}




// ===================================
// Estimator methods

// Set the flag
void set_estimator_ready() {
	updated = true;

	return;
}

/**
 * Estimate state
 */
static void estimate_state() {
	// Save the old state before the update
	for (int i = 0; i < 3; i++) {
		X_old[i] = X[i];
	}

	arm_mat_mult_f32(&O_invm, &Ybuffm, &Xm);
	return;
}


// =============================================
// Getters and Setters
static float estimatorDD_getAlpha_() {
	return alpha_;
}

static void estimatorDD_setAlpha_(float a) {
	alpha_ = a;
}

static float estimatorDD_getBeta_() {
	return beta_;
}

static void estimatorDD_setBeta_(float b) {
	beta_ = b;
}

static float estimatorDD_getCtrl_dd_() {
	return ctrl_dd_;
}

static float estimatorDD_getCtrl_ddd_() {
	return ctrl_ddd_;
}


/**
 * Estimate params
 */
static void estimate_params(float TotalTime, float c_dd, float c_ddd) {
	
	// Get the current value from the global variables
	float alpha = estimatorDD_getAlpha_();
	float beta = estimatorDD_getBeta_();

	// Local variables
	float alpha_new = alpha;
	float beta_new = beta;

	float ctrl_dd_scaled = c_dd / 65535.0f;
	float ctrl_ddd_scaled = c_ddd / 65535.0f;
	
	if (c_dd > 1.0f){
		if (ctrl_dd_active) {
			switch (Step) {
				case 1:
					estimate_least_squares = false;
					alpha_new = alpha + 1.0f/TotalTime * (X[1] - X_old[1])- (alpha + ctrl_dd_scaled * beta); 
					beta_new = beta;
					//u = 41000.0f/65535.0f;
					//estimatorDDSetControl(u);
					Step++;
					break;
				case 2:
					alpha_new = alpha + ctrl_ddd_scaled / (TotalTime * (ctrl_ddd_scaled - ctrl_dd_scaled))* ((X[1] - X_old[1]) - TotalTime * (alpha + ctrl_dd_scaled * beta)); 
					beta_new = beta + 1 / (TotalTime * (ctrl_dd_scaled - ctrl_ddd_scaled)) * ((X[1] - X_old[1]) - TotalTime * (alpha + ctrl_dd_scaled * beta));
					estimate_least_squares = false;
					Step++;
					break;  
				default:
					alpha_new = alpha + gamma1 * ((X[1] - X_old[1])- TotalTime * (alpha + ctrl_dd_scaled * beta));
					//beta_new = beta + gamma2 * ctrl_dd_scaled * ((X[1] - X_old[1])- TotalTime * (alpha + ctrl_dd_scaled * beta));
					break;  
			}
		}
		else {
			float a_new_part0 = gamma1 * TotalTime * (alpha + ctrl_dd_scaled * beta);
			float a_new_part1 =  gamma1 * (X[1] - X_old[1]);
			alpha_new = alpha - a_new_part0  +  a_new_part1;

			//float b_new_part0 = gamma2 * ctrl_dd_scaled * TotalTime * (alpha + ctrl_dd_scaled * beta);
			//float b_new_part1 = gamma2 * ctrl_dd_scaled * (X[1] - X_old[1]);
			//beta_new = beta - b_new_part0 +  b_new_part1;
		} 
	}
	if (beta_new < 1e-6f){
		beta_new = 1e-6f;
	}
	Step++;
	if (Step == 3000){
		DEBUG_PRINT("Input PID [ %.3f, %.3f, %.3f]\n", 
				(double)alpha_new, (double)beta_new, (double)ctrl_dd_scaled);
		Step=4;  
	}

	// Update the state
	estimatorDD_setAlpha_(alpha_new);
	estimatorDD_setBeta_(beta_new);

	return;
}


/**
 * Compute control value
 */
static float u_p;
static float u_d;
static float u_a;
static void compute_ctrl() {
	float u_fb = 0;

	float alpha = estimatorDD_getAlpha_();
	float beta = estimatorDD_getBeta_();

	// Update the control gain
	Kdd[0] = - P1 * P2;
	Kdd[1] = P1 + P2;
	Kdd[2] = 0.0f; 

	u_p = Kdd[0] * (X[0] - Tracking[0]);
	u_d = Kdd[1] * (X[1] - Tracking[1]);
	u_a = Kdd[2] * (X[2] - Tracking[2]);

	u_fb = u_p + u_d + u_a;	
	//alpha=-12.0f;
	u = (1.0f / beta) * (-alpha + u_fb);

	if (u < 0.0f) {
		u = 0.0f;
	}
	if (u > 1.0f) {
		u = 1.0f;
	}
    //if (Step == 4){
	//	DEBUG_PRINT("Input DD [ %.6f, %.6f, %.6f]\n", 
	//			(double)alpha, (double)beta, (double)u);
	//}
	U = u * 65535.0f;
	if (!estimate_least_squares){
		if (!Land){
			estimatorDDSetControl(U);
			if (Step == 4){
				DEBUG_PRINT("[ %.6f, %.6f, %.6f]\n", 
						(double)alpha, (double)beta, (double)u);
			}
		} else{
			estimatorDDSetControl(0.0); 
		} 
	}
}

/** 
 * Insert the k-th measurement in the buffer
 */
static void insert_newmeas_batch(float y, float stamp, int k) {
	if (k < 0) {
		// Error
	}
	int index = (BUFF_SIZE - 1) - (k % BUFF_SIZE);
	Ybuff[index] = y;
	Tbuff[index] = stamp;
}

static void insert_newmeas_circ(float y, float stamp) {

	for (int index = BUFF_SIZE-1; index > 0; index--) { 
		Ybuff[index] = Ybuff[index-1];
		Tbuff[index] = Tbuff[index-1];       
	}
	Ybuff[0] = y;
	Tbuff[0] = stamp;
}

static void update_O(float t[BUFF_SIZE]) {
	int i;
	for (i = 0; i < BUFF_SIZE; i++) {
		O[(i*STATE_SIZE) + 1] = -(t[i]);
		O[(i*STATE_SIZE) + 2] = 0.5f * (t[i] * t[i]);
	}
}


static void finalize_data() {

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
	// TODO: Either make everything static with void calls,
	// 	either pass the values inside all the chain of calls
	eval_pseudoinv(&O_invm);

	/*
	   static int counter = 0;
	   if (counter == 150) {
	   DEBUG_PRINT("[ %.3f, %.3f, %.3f, %.3f , %.3f, %.3f, %.3f]\n", 
	   (double)*O_invm.pData, (double)*(O_invm.pData + 1), (double)*(O_invm.pData + 2), 
	   (double)*(O_invm.pData + 3), (double)*(O_invm.pData + 4), (double)*(O_invm.pData + 5),
	   (double)*(O_invm.pData + 6));
	   counter = 0;
	   }
	   counter++;
	   */
}

/**
 * Estimator step function
 */
void DDEstimator_step_circ(float y, float stamp) {
	if (!isInit) {
		estimatorDDInit();
	}

	/*
	   static int counter = 0;
	   if (counter == 150 || counter == 151) {
	   DEBUG_PRINT("[ %.3f, %.3f, %.3f, %.3f , %.3f]\n", (double)Tbuff[0], (double)Tbuff[1], (double)Tbuff[2], (double)Tbuff[3], (double)Tbuff[4]);
	   if (counter == 151){    
	   counter = 0;
	   }
	   }
	   counter++;
	   */



	// Update the buffer
	insert_newmeas_circ(y, stamp);


	if (Nmeas >= BUFF_SIZE) {
		// State Estimation
		finalize_data();
		estimate_state();

		// Estimate Parameters
		float cdd = estimatorDD_getCtrl_dd_();
		float cddd = estimatorDD_getCtrl_ddd_();
		estimate_params(TotalTime, cdd, cddd);

		// Control
		if (ctrl_dd_active && Step > 1) {	
			compute_ctrl();
		}
		set_estimator_ready();
	}
	else{
		Nmeas++;
	}    
}

void DDEstimator_step_batch(float y, float stamp) {
	if (!isInit) {
		estimatorDDInit();
	}

	// Update the buffer
	insert_newmeas_batch(y, stamp, Nmeas);
	Nmeas++;

	if (Nmeas == BUFF_SIZE) {
		Nmeas = 0;
		finalize_data();

		estimate_state();
		
		// Estimate Parameters
		float cdd = estimatorDD_getCtrl_dd_();
		float cddd = estimatorDD_getCtrl_ddd_();
		estimate_params(TotalTime, cdd, cddd);

		if (ctrl_dd_active && Step > 1) {	
			compute_ctrl();
		}
		set_estimator_ready();
	}
}


/**
 * Feed the DD controller with the state estimate.
 * (It will be necessary to disable the callback from the sensor to avoid
 * overlapping the estimation (and race conditions on the variables)
 *
 * If this is a feature that would be embedded in the final version we 
 * can work on a logic for all these workarounds...
 */
void estimatorDDFeedState(float z, float zd, uint64_t us_timestamp) {

	// Measure the time to check whether the trigger is really periodic

	dt_ms = (float)(us_timestamp - us_timestamp_old) / 1e6f;
	us_timestamp_old = us_timestamp;
	for (int i = 0; i < 3; i++) {
		X_old[i] = X[i];
	}

	X[0] = z;
	X[1] = zd;
	X[2] = 0.0;

	// The main loop is supposed to spin at 1Khz. Clearly, the information from the Z
	// is not only provided by the camera, but takes into account the filter prediction
	// capabilities.
	float T = dt_ms;

	// Estimate Parameters
	float cdd = estimatorDD_getCtrl_dd_();
	float cddd = estimatorDD_getCtrl_ddd_();
	estimate_params(T, cdd, cddd);

	if (ctrl_dd_active && Step > 1) {
		compute_ctrl();
	}

	set_estimator_ready();

}


// ====================================
void init_O() {
	int i;
	for (i = 0; i < BUFF_SIZE; i++) {
		O[(i*STATE_SIZE) + 0] = 1;
		O[(i*STATE_SIZE) + 1] = -(i * TS);
		O[(i*STATE_SIZE) + 2] = 0.5f * (i * i * TS2);
	}  
}

/**
 * Initialization function 
 */
void estimatorDDInit(void) {
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
	alpha_new = 0.0f;
	beta_ = 18.7291;
	// Initialize the DD Library
	DataDriven_initialize();

	isInit = true;
	return;
}



bool estimatorDDTest(void) {
	return isInit;
}

// This function is triggered by the arrival of new measurements
bool estimatorDDNewMeasurement(const positionMeasurement_t *pos) {

	msg_counter = msg_counter + 1;

	// Measure the timestamp
	timestamp = pos->t; // Time in second
	dt_ms = (timestamp - timestamp_old) * 1e3f;
	t_s = timestamp;
	/*
	   static int counter = 0;

	   if (counter == 400) {
	   DEBUG_PRINT("Timestamp = %.6f s \n", (double)timestamp);
	   DEBUG_PRINT("Old = %.6f \n", (double)timestamp_old);
	   DEBUG_PRINT("DT = %.6f ms \n", (double)dt_ms);
	   counter = 0;
	   }
	   counter++;
	   */
	timestamp_old = timestamp;

	state_z = pos->z;

	// Do something with the new measurement 
	DDEstimator_step_circ(state_z, t_s);

	//	if (msg_counter == 1000) {
	//		DEBUG_PRINT("\n");
	//		DEBUG_PRINT("Execution Time = %llu us\n", ctrl_exetime);
	//
	//		DEBUG_PRINT("Y_BUFF	= [");
	//		for (int i = 0; i < 5; i++) 
	//			DEBUG_PRINT("%.3f, ", (double)Ybuff[i]);
	//		DEBUG_PRINT("]\n");
	//
	//		DEBUG_PRINT("X_est = [%.3f, %.3f, %.3f] \n", (double)X[0], (double)X[1], (double)X[2]);
	//		msg_counter = 0;
	//	}


	return true;
}



float estimatorDDGetEstimatedZ() {
	return state_z;
}

/**
 * Check whether the estimator is ready. In that case
 * reset the flag and return 'true'. Otherwise, return 'false'.
 */
bool estimatorDDHasNewEstimate() {
	bool out;

	//xSemaphoreTake(mutex, portMAX_DELAY);
	out = updated;
	if (out)
		updated = false; // Reset the flag
	//xSemaphoreGive(mutex);

	return out;
}

void estimatorDDSetControl(const float u) {
	ctrl_ddd_ = ctrl_dd_;
	ctrl_dd_ = u;
}

float estimatorDDGetControl() {
	return ctrl_dd_;
}

void estimatorDDParamLeastSquares(void) {   
	ctrl_dd_active = true;
} 



// Logging variables
//
	LOG_GROUP_START(estimator_dd)
	LOG_ADD(LOG_FLOAT, est_x, &X[0])
	LOG_ADD(LOG_FLOAT, est_xd, &X[1])
	LOG_ADD(LOG_FLOAT, est_xdd, &X[2])
	LOG_ADD(LOG_FLOAT, est_alpha, &alpha_)
	LOG_ADD(LOG_FLOAT, est_beta, &beta_)
	LOG_ADD(LOG_FLOAT, sens_dt_ms, &dt_ms)
	LOG_ADD(LOG_FLOAT, est_dt_s, &TotalTime)
LOG_GROUP_STOP(estimator_dd)
	/*
	   LOG_GROUP_START(controller_dd)
	   LOG_ADD(LOG_FLOAT, up, &u_p)
	   LOG_ADD(LOG_FLOAT, ud, &u_d)
	   LOG_GROUP_STOP(controller_dd)
	   */
	PARAM_GROUP_START(controller_dd)
	PARAM_ADD(PARAM_FLOAT, ctrl_ddP1, &P1)
	PARAM_ADD(PARAM_FLOAT, ctrl_ddP2, &P2)
	PARAM_ADD(PARAM_FLOAT, ctrl_ddg1, &gamma1)
	PARAM_ADD(PARAM_FLOAT, ctrl_ddg2, &gamma2)  
	PARAM_ADD(PARAM_FLOAT, ctrl_ddA, &alpha_)
	PARAM_ADD(PARAM_FLOAT, ctrl_ddB, &beta_)   
	PARAM_ADD(PARAM_FLOAT, ctrl_ddTr, &Tracking[0])
	PARAM_ADD(PARAM_FLOAT, ctrl_ddLd, &Land)
PARAM_GROUP_STOP(controller_dd)

