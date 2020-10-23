//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef EQUATION_MODULE_H
#define EQUATION_MODULE_H

#include "fv3d_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
typedef struct Variable
{
    string_t name;
} Variable_t;

typedef struct Variables
{
    int n_sol_variables;
    int n_dep_variables;
    int n_tot_variables;

    Variable_t *sol_variables;
    Variable_t *dep_variables;
} Variables_t;

typedef void (*void_update_fp_t)( double );
typedef void (*void_update_gradients_fp_t)();
typedef void (*void_calc_exact_fp_t)( int, double, double, double );
typedef void (*void_calc_timestep_fp_t)( double );
typedef void (*void_calc_flux_fp_t)();

extern Variables_t *all_variables;

extern void_update_fp_t update_function_pointer;
extern void_update_gradients_fp_t update_gradients_function_pointer;
extern void_calc_exact_fp_t calc_exact_function_pointer;
extern void_calc_timestep_fp_t calc_time_step_function_pointer;
extern void_calc_flux_fp_t calc_flux_function_pointer;

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_define();

int add_sol_variable( Variables_t *variables, string_t name );
int add_dep_variable( Variables_t *variables, string_t name );

#endif /* EQUATION_MODULE_H */