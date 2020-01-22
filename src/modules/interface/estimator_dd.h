#ifndef __ESTIMATOR_DD_H__
#define __ESTIMATOR_DD_H__

#include <stdint.h>
#include "stabilizer_types.h"

// I should add a structure to clean the interface...

void DDcontroller_Init(void);

bool DDestimator_Check_Init(void);

// Push the new sensor reading in the estimator
bool DDcontroller_NewMeasurement(const positionMeasurement_t *pos);

// Get the estimated value
float DDestimator_Get_StateEstimate();

// Get the control computed by the DD
float DDcontroller_Get_Control();

// Feed the state into the estimator
void DDestimator_Feed_State(float z, float zd, uint64_t t);

void DDcontroller_Set_Control(const float v);

// Check whether the estimator has a finished the estimation cycle
bool DDestimator_Check_NewMeasurement();

// Start alpha and beta estimator
void DDcontroller_Set_ControllerReady(void);

#endif
