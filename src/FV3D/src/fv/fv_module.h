//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef FV_MODULE_H
#define FV_MODULE_H

#include "fv3d/fv3d_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
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

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void fv_define();
void reconstruction_define();
void limiter_define();

void fv_time_derivative(double t);

#endif /* FV_MODULE_H */