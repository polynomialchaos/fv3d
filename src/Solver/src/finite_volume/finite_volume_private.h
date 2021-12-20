/*******************************************************************************
 * @file finite_volume_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef FINITE_VOLUME_PRIVATE_H
#define FINITE_VOLUME_PRIVATE_H

#include "solver/solver_module.h"

typedef double (*double_limiter_ft)(int i_cell, int i_var, double slope);
extern double_limiter_ft limiter_function_pointer;

typedef void (*void_reconstruction_ft)(double t);
extern void_reconstruction_ft reconstruction_function_pointer;

extern void_update_ft update_function_pointer;
extern void_update_gradients_ft update_gradients_function_pointer;

/*******************************************************************************
 * @brief Allocate the solver data
 ******************************************************************************/
data_t *allocate_data();

/*******************************************************************************
 * @brief Calculate (reconstruct) gradients
 ******************************************************************************/
void calc_gradients();

/*******************************************************************************
 * @brief Deallocate the solver data
 * @param data
 ******************************************************************************/
void deallocate_data(data_t *data);

/*******************************************************************************
 * @brief Barth-Jespersenn limiter calculation (return 0-1)
 * @param i_cell
 * @param i_var
 * @param slope
 * @return double
 ******************************************************************************/
double limiter_barth_jespersenn(int i_cell, int i_var, double slope);

/*******************************************************************************
 * @brief None limiter calculation (return 1)
 * @param i_cell
 * @param i_var
 * @param slope
 * @return double
 ******************************************************************************/
double limiter_none(int i_cell, int i_var, double slope);

/*******************************************************************************
 * @brief First-order reconstruction
 ******************************************************************************/
void reconstruction_first_order();

/*******************************************************************************
 * @brief Second-order (linear) reconstruction
 ******************************************************************************/
void reconstruction_linear();

/*******************************************************************************
 * @brief Set/initialize the solution pointer
 * @param t
 ******************************************************************************/
void set_solution();

/*******************************************************************************
 * @brief Update the given variable across all domains
 * @param phi_local
 ******************************************************************************/
void update_parallel(double *phi_local);

#endif /* FINITE_VOLUME_PRIVATE_H */