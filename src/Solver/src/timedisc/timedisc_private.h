/*******************************************************************************
 * @file timedisc_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef TIMEDISC_PRIVATE_H
#define TIMEDISC_PRIVATE_H

#include "solver/solver_module.h"

typedef void (*void_timestep_ft)(int iter, double t, double dt);
extern void_timestep_ft time_step_function_pointer;

extern int solver_is_viscous_dt;
extern int solver_is_transient;

extern int n_iter_inner;
extern int n_iter_lsoe;
extern int n_bdf_stages;
extern double **phi_old;

/*******************************************************************************
 * @brief Numerical jacobian routine
 * @param n_var
 * @param n_cells
 ******************************************************************************/
void calc_jacobian_numerical(int n_var, int n_cells);

/*******************************************************************************
 * @brief Implicit time discretizazion routine (Newton)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void time_step_newton(int iter, double t, double dt);

/*******************************************************************************
 * @brief Matrix vector routine (called by solver)
 * @param x
 * @param b
 * @param n_var
 * @param n_cells
 * @return int
 ******************************************************************************/
int matrix_vector_numerical(double *x, double *b, size_t n_var, size_t n_cells);

/*******************************************************************************
 * @brief Explicit time discretizazion routine (LSERKW2)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void time_step_lserkw2(int iter, double t, double dt);

#endif /* TIMEDISC_PRIVATE_H */