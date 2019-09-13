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
#define BUFF_SIZE (5)
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
float O[BUFF_SIZE * STATE_SIZE] = 
{	
	1.0, 	0.0, 	0.0,
	1.0, 	-TS, 	TS2/2,
	1.0, 	-2*TS, 	2*TS2,
	1.0, 	-3*TS,	9/2*TS2,
	1.0, 	-4*TS, 	8*TS2,
};

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
static uint64_t timestamp;
static uint64_t timestamp_old;
//static uint64_t timestamp_ctrl;

float dt_ms;
float t_s;
float dt_ms_cum = 0;
static uint32_t msg_counter = 0;

// ====================================
// Estimator State 
static float state_z;
static float alpha = 0.0;
static float beta;

// Estimator Parametrs
static float gamma1 = 3.0f;

// Control gain
static float L = 900;
static float Kdd[3];
static float ctrl_dd;
static float Tracking[] = {0.0, 0.0, 0.0};

static bool updated = false;

// ====================================
// Filter Data
static int Nmeas;

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
	beta = 1.0f / droneMass;
	alpha = alpha - gamma1 * (Tbuff[0] - Tbuff[BUFF_SIZE - 1]) * (alpha + ctrl_dd * beta) +  gamma1 * (X[1] - X_old[1]); 
	return;
}


/**
 * Compute control value
 */
static void compute_ctrl() {
	int i;
	float u_fb = 0;
	for (i = 0; i < 3; i++) {
		u_fb += Kdd[i] * (X[i] - Tracking[i]);
	}
	ctrl_dd = (1.0f / beta) * (-alpha + u_fb);
}

/** 
 * Insert the k-th measurement in the buffer
 */
static void insert_newmeas(float y, float stamp, int k) {
	if (k < 0) {
		// Error
	}
	int index = (BUFF_SIZE - 1) - (k % BUFF_SIZE);
	Ybuff[index] = y;
	Tbuff[index] = stamp;
}


static void update_O(float t[BUFF_SIZE]) {
	int i;
	for (i = 0; i < BUFF_SIZE; i++) {
		O[(i*STATE_SIZE) + 1] = -(i * t[i]);
		O[(i*STATE_SIZE) + 2] = 0.5f * (t[i] * t[i]);
	}

}

static void finalize_data() {
	int i;

	// Finalize the DT vector, computing the differences
	// [0, dt1, dt1 + dt2, ...]
	for (i = 0; i < BUFF_SIZE; i++) {
		Tbuff[i] = Tbuff[0] - Tbuff[i];
	}

	// Update the Observability matrix
	update_O(Tbuff);

	// Update the pseduoinverse matrix
	// TODO: Either make everything static with void calls,
	// 	either pass the values inside all the chain of calls
	eval_pseudoinv(&O_invm);
}

/**
 * Estimator step function
 */
void DDEstimator_step(float y, float stamp) {
	if (!isInit) {
		estimatorDDInit();
	}

	// Update the buffer
	insert_newmeas(y, stamp, Nmeas);
	Nmeas++;

	if (Nmeas == BUFF_SIZE) {
		Nmeas = 0;
		finalize_data();

		estimate_state();
		estimate_params();	
		compute_ctrl();
		set_estimator_ready();
	}
}

// ====================================


/**
 * Initialization function 
 */
void estimatorDDInit(void) {
	if (isInit)  {
		return;
	}

	Kdd[0] = -L*L;
	Kdd[1] = 0.5f * TS * L * L  - 2.0f * L;
	Kdd[2] = 0.0f; 

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
	timestamp = usecTimestamp(); // Time in microseconds
	dt_ms = (timestamp - timestamp_old) / 1e3;
	t_s = timestamp / 1e6;
	timestamp_old = timestamp;

	state_z = pos->z;

	// Do something with the new measurement 
	DDEstimator_step(state_z, t_s);

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
	ctrl_dd = u;
}

float estimatorDDGetControl() {
	return ctrl_dd;
}



// Logging variables
//
LOG_GROUP_START(estimator_dd)
LOG_ADD(LOG_FLOAT, est_x, &X[0])
LOG_ADD(LOG_FLOAT, est_xd, &X[1])
LOG_ADD(LOG_FLOAT, est_xdd, &X[2])
LOG_ADD(LOG_FLOAT, est_alpha, &alpha)
LOG_ADD(LOG_FLOAT, sens_dt_ms, &dt_ms)
LOG_GROUP_STOP(estimator_dd)
