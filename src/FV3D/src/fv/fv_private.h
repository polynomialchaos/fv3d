//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef FV_PRIVATE_H
#define FV_PRIVATE_H

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
    // real,   allocatable,    target  :: phi_total(:,:)           !< solution vector of all variables

    // real,   allocatable,    target  :: grad_phi_total_x(:,:)    !< gradient vector of all variables (x)
    // real,   allocatable,    target  :: grad_phi_total_y(:,:)    !< gradient vector of all variables (y)
    // real,   allocatable,    target  :: grad_phi_total_z(:,:)    !< gradient vector of all variables (z)

    // real,   allocatable,    target  :: phi_total_left(:,:)      !< reconstruction vector of all variables (master)
    // real,   allocatable,    target  :: phi_total_right(:,:)     !< reconstruction vector of all variables (slave)

    // real,   pointer                 :: phi(:,:)                 !< solution vector of active variables
    // real,   allocatable             :: phi_dt(:,:)              !< time derivative solution vector of active variables
    // real,   allocatable             :: flux(:,:)                !< flux vector of variables

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void fv_define();
void fv_time_derivative( double t );

#endif /* FV_PRIVATE_H */