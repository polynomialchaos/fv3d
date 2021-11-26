/*******************************************************************************
 * @file fv_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef FV_MODULE_H
#define FV_MODULE_H

#include "fv3d/fv3d_module.h"

typedef void (*void_reconstruction_ft)();
extern void_reconstruction_ft reconstruction_function_pointer;

typedef double (*double_limiter_ft)(int i_cell, int i_var, double slope);
extern double_limiter_ft limiter_function_pointer;

typedef void (*void_calc_flux_ft)();
extern void_calc_flux_ft calc_flux_function_pointer;

typedef void (*void_update_ft)(double t);
extern void_update_ft update_function_pointer;

typedef void (*void_update_gradients_ft)();
extern void_update_gradients_ft update_gradients_function_pointer;

typedef void (*void_calc_exact_ft)(int id, double t, double *x, double *phi);
extern void_calc_exact_ft calc_exact_function_pointer;

extern double *phi_total;
extern double *grad_phi_total_x;
extern double *grad_phi_total_y;
extern double *grad_phi_total_z;

extern double *phi_total_left;
extern double *phi_total_right;

extern double *phi_dt;
extern double *flux;

void fv_time_derivative(double t);

/*******************************************************************************
 * @brief Calculate (reconstruct) gradients
 ******************************************************************************/
void calc_gradients();

/*******************************************************************************
 * @brief Define fv
 ******************************************************************************/
void fv_define();

/*******************************************************************************
 * @brief Finalize fv
 ******************************************************************************/
void fv_finalize();

/*******************************************************************************
 * @brief Initialize fv
 ******************************************************************************/
void fv_initialize();

/*******************************************************************************
 * @brief Barth-Jespersenn limiter calculation (return 0-1)
 * @param i_cell
 * @param i_var
 * @param slope
 * @return double
 ******************************************************************************/
double limiter_barth_jespersenn(int i_cell, int i_var, double slope);

/*******************************************************************************
 * @brief Define limiter
 ******************************************************************************/
void limiter_define();

/*******************************************************************************
 * @brief Finalize limiter
 ******************************************************************************/
void limiter_finalize();

/*******************************************************************************
 * @brief Initialize limiter
 ******************************************************************************/
void limiter_initialize();

/*******************************************************************************
 * @brief None limiter calculation (return 1)
 * @param i_cell
 * @param i_var
 * @param slope
 * @return double
 ******************************************************************************/
double limiter_none(int i_cell, int i_var, double slope);

/*******************************************************************************
 * @brief Define reconstruction
 ******************************************************************************/
void reconstruction_define();

/*******************************************************************************
 * @brief Finalize reconstruction
 ******************************************************************************/
void reconstruction_finalize();

/*******************************************************************************
 * @brief First-order reconstruction
 ******************************************************************************/
void reconstruction_first_order();

/*******************************************************************************
 * @brief Initialize reconstruction
 ******************************************************************************/
void reconstruction_initialize();

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

#endif /* FV_MODULE_H */