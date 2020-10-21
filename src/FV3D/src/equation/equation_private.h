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
extern int n_variables;

    // character(len=_STRLEN_)             :: equation = 'Navier-Stokes'               !< The equation to solve

    // type variable_t
    //     character(len=:),   allocatable :: name                                     !< variable name
    //     logical                         :: is_active = .true.                       !< variable active flag
    // end type variable_t

    // integer                             :: n_variables = 0                          !< number of variables (only active)
    // integer                             :: n_tot_variables = 0                      !< number of total variables
    // type(variable_t),       allocatable :: variables(:)                             !< variables array

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