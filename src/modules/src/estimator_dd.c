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


#define  droneMass (0.032f);

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
//static uint64_t timestamp_ctrl;

float dt_ms;
float t_s;
float dt_ms_cum = 0;
static uint32_t msg_counter = 0;

// ====================================
// Estimator State 
static float state_z;
static float alpha = 0.0f;
static float alpha_new = 0.0f;
static float beta = 2.8577f*1e-4f;
// /droneMass;

// Estimator Parametrs
static float gamma1 = 1.0f;

// Control gain
static float P1 = 1.0f;
static float P2 = 1.0f;
static float Kdd[3];
static float ctrl_dd;
static float ctrl_ddd;
static float Tracking[] = {1.2f, 0.0, 0.0};

// Control Placeholder
static float u = 0;

// Step Counter
static int Step=3;
static bool ctrl_dd_active = false;
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
	xSemaphoreTake(mutex, portMAX_DELAY);
	updated = true;
	xSemaphoreGive(mutex);

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


/**
 * Estimate params
 */
static void estimate_params() {
	if (ctrl_dd_active) {
		Step++;
		switch (Step) {
			case 1:
				alpha_new = alpha - (alpha + ctrl_dd * beta) +  1.0f/TotalTime * (X[1] - X_old[1]); 
				beta = beta;
				alpha = alpha_new; 
				u = 41000.0f;
				estimatorDDSetControl(u);
				break;
			case 2:
				alpha_new = alpha + ctrl_ddd / (TotalTime * (ctrl_ddd - ctrl_dd))* ((X[1] - X_old[1]) - TotalTime * (alpha + ctrl_dd * beta)); 
				beta = beta + ctrl_ddd / (TotalTime * (ctrl_dd - ctrl_ddd)) * ((X[1] - X_old[1]) - TotalTime * (alpha + ctrl_dd * beta));
				alpha = alpha_new;
				break;  
			default:
				alpha_new = alpha - gamma1 * (TotalTime) * (alpha + ctrl_dd * beta) +  gamma1 * (X[1] - X_old[1]);
				// beta = beta - gamma1 * ctrl_dd * TotalTime * (alpha + ctrl_dd * beta) +  gamma1 * ctrl_dd* (X[1] - X_old[1]);
				alpha = alpha_new;
				break;  
		}
	}
	else {
		alpha_new = alpha - gamma1 * (TotalTime) * (alpha + ctrl_dd * beta) +  gamma1 * (X[1] - X_old[1]);
		// beta = beta - gamma1 * ctrl_dd * TotalTime * (alpha + ctrl_dd * beta) +  gamma1 * ctrl_dd* (X[1] - X_old[1]);
		alpha = alpha_new; 
	}
	if (beta<1e-6f){
		beta = 1e-6f;
	}
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

	// Update the control gain
	Kdd[0] = -P1 * P2;
	Kdd[1] = P1 + P2;
	Kdd[2] = 0.0f; 

	u_p = Kdd[0] * (X[0] - Tracking[0]);
	u_d = Kdd[1] * (X[1] - Tracking[1]);
	u_a = Kdd[2] * (X[2] - Tracking[2]);

	u_fb = u_p + u_d + u_a;	
	alpha=-12.100f;
	u = (1.0f / beta) * (-alpha + u_fb);

	if (u < 0.0f) {
		u = 0.0f;
	}
	if (u > 65535.0f) {
		u = 65535.0f;
	}

	estimatorDDSetControl(u);
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
		O[(i*STATE_SIZE) + 1] = -(i * t[i]);
		O[(i*STATE_SIZE) + 2] = 0.5f * (i * i * t[i] * t[i]);
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
	TotalTime = DTbuff[1];


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
		finalize_data();
		estimate_state();
		estimate_params();
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
		estimate_params();
		if (ctrl_dd_active && Step > 1) {	
			compute_ctrl();
		}
		set_estimator_ready();
	}
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
	Kdd[1] = -P1 + P2;
	Kdd[2] = 0.0f; 

	DEBUG_PRINT("DD Controller Gain: [%.3f, %.3f] \n", (double)Kdd[0], (double)Kdd[1]);
	init_O();	

	eval_pseudoinv(&O_invm);

	mutex = xSemaphoreCreateMutex();

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
	DDEstimator_step_batch(state_z, t_s);

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
	ctrl_ddd=ctrl_dd;
	ctrl_dd = u;
}

float estimatorDDGetControl() {
	return ctrl_dd;
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
	LOG_ADD(LOG_FLOAT, est_alpha, &alpha)
	LOG_ADD(LOG_FLOAT, est_beta, &beta)
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
	PARAM_ADD(PARAM_FLOAT, ctrl_dd_g, &gamma1)
PARAM_GROUP_STOP(controller_dd)

