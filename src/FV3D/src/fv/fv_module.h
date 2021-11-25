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

typedef void (*void_reconstruction_fp_t)();
extern void_reconstruction_fp_t reconstruction_function_pointer;

typedef double (*double_limiter_fp_t)(int i_cell, int i_var, double slope);
extern double_limiter_fp_t limiter_function_pointer;

typedef void (*void_calc_flux_fp_t)();
extern void_calc_flux_fp_t calc_flux_function_pointer;

typedef void (*void_update_fp_t)(double t);
extern void_update_fp_t update_function_pointer;

typedef void (*void_update_gradients_fp_t)();
extern void_update_gradients_fp_t update_gradients_function_pointer;

typedef void (*void_calc_exact_fp_t)(int id, double t, double *x, double *phi);
extern void_calc_exact_fp_t calc_exact_function_pointer;

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
 * @brief Define reconstruction
 ******************************************************************************/
void reconstruction_define();

/*******************************************************************************
 * @brief Finalize reconstruction
 ******************************************************************************/
void reconstruction_finalize();

/*******************************************************************************
 * @brief Initialize reconstruction
 ******************************************************************************/
void reconstruction_initialize();

#endif /* FV_MODULE_H */