//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef TIMEDISC_PRIVATE_H
#define TIMEDISC_PRIVATE_H

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
typedef void (*void_timedisc_fp_t)( double t, double dt, int iter );
extern void_timedisc_fp_t time_step;

    // character(len=_STRLEN_) :: time_step= 'Explicit'        !< The timestep mehtod
    // integer                 :: max_iter = 100000            !< The maximum number of iterations
    // logical                 :: transient = .true.           !< The flag wheter to be transient or steady-state
    // real                    :: abort_residual = 1e-10       !< The abort residual
    // real                    :: t_start = 0.0                !< The start time
    // real                    :: t_end = 0.2                  !< The end time
    // real                    :: dt = 1e-10                   !< The timestep

    // procedure(),    pointer :: time_step_routine => null()  !< TimeStep routine
    // logical                 :: is_explicit = .true.         !< Flag if explicit or implicit

    // real                    :: t                            !< Current time
    // integer                 :: iter                         !< Current iteration number

    // logical                 :: is_viscous_dt = .false.      !< viscous timestep flag

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void timedisc_define();
void timedisc();

#endif /* TIMEDISC_PRIVATE_H */