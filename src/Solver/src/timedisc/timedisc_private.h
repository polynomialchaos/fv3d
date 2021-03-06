/*******************************************************************************
 * @file timedisc_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef TIMEDISC_PRIVATE_H
#define TIMEDISC_PRIVATE_H

#include "solver/solver_module.h"

typedef void (*void_timestep_ft)(int iter, double t, double dt);
extern void_timestep_ft timestep_function_pointer;

extern double_calc_timestep_ft calc_timestep_function_pointer;

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
 * @brief Matrix vector routine (called by solver)
 * @param x
 * @param b
 * @param n_var
 * @param n_cells
 * @return int
 ******************************************************************************/
int matrix_vector_numerical(double *x, double *b, size_t n_var, size_t n_cells);

/*******************************************************************************
 * @brief Implicit time discretizazion routine (Newton)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void timestep_newton(int iter, double t, double dt);

/*******************************************************************************
 * @brief Explicit time discretizazion routine (LSERKW2)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void timestep_lserkw2(int iter, double t, double dt);

#endif /* TIMEDISC_PRIVATE_H */