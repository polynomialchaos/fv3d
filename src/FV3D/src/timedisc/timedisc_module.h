/*******************************************************************************
 * @file timedisc_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef TIMEDISC_MODULE_H
#define TIMEDISC_MODULE_H

#include "fv3d/fv3d_module.h"

typedef void (*void_timestep_fp_t)(int iter, double t, double dt);
extern void_timestep_fp_t time_step_function_pointer;

typedef double (*double_calc_timestep_fp_t)();
extern double_calc_timestep_fp_t calc_time_step_function_pointer;

extern int is_viscous_dt;
extern int is_transient;

void timedisc_define();
void timedisc();

#endif /* TIMEDISC_MODULE_H */