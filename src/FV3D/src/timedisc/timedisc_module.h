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

typedef void (*void_timestep_ft)(int iter, double t, double dt);
extern void_timestep_ft time_step_function_pointer;

typedef double (*double_calc_timestep_ft)();
extern double_calc_timestep_ft calc_time_step_function_pointer;

extern int is_viscous_dt;
extern int is_transient;

/*******************************************************************************
 * @brief Time discretizazion routine
 ******************************************************************************/
void timedisc();

/*******************************************************************************
 * @brief Define timedisc
 ******************************************************************************/
void timedisc_define();

/*******************************************************************************
 * @brief Finalize timedisc
 ******************************************************************************/
void timedisc_finalize();

/*******************************************************************************
 * @brief Initialize timedisc
 ******************************************************************************/
void timedisc_initialize();

#endif /* TIMEDISC_MODULE_H */