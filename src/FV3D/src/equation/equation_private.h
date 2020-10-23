//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef EQUATION_PRIVATE_H
#define EQUATION_PRIVATE_H

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

    Variable_t *sol_variables;
    Variable_t *dep_variables;
} Variables_t;

extern Variables_t *all_variables;

// procedure(),    pointer             :: exact_func_routine => null()             !< calculate exact function
// procedure(),    pointer             :: update_routine => null()                 !< update cell and boundary values
// procedure(),    pointer             :: update_gradients_routine => null()       !< update gradients at boundaries
// procedure(),    pointer             :: calc_time_step_routine => null()         !< time step routine
// procedure(),    pointer             :: calc_flux_routine => null()              !< flux calculation routine

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_define();

#endif /* EQUATION_PRIVATE_H */