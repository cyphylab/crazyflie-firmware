#ifndef __LFD_H__
#define __LFD_H__

#include <stdint.h>
#include "stabilizer_types.h"


// I should add a structure to clean the interface...

void LfDcontroller_Init(void);

bool LfDestimator_Check_Init(void);

// Push the new sensor reading in the estimator
bool LfDcontroller_NewMeasurement(const positionMeasurement_t *pos);

// Get the estimated value
float LfDestimator_Get_StateEstimate();

// Get the control computed by the LfD
float LfDcontroller_Get_Control();

// Feed the state into the estimator
void LfDestimator_Feed_State(float z, float zd, uint64_t t);

void LfDcontroller_Set_Control(const float v);

// Check whether the estimator has a finished the estimation cycle
bool LfDestimator_Check_NewMeasurement();

// Start alpha and beta estimator
void LfDcontroller_Set_ControllerReady(void);

#endif
