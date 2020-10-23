//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef RESTART_MODULE_H
#define RESTART_MODULE_H

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
// ! logical                 :: use_restart = .false.    !< The flag to start from restart

// ! integer                 :: iter_restart = 0         !< Current iteration number
// ! real                    :: time_restart = 0.0       !< Current time

// ! real,   allocatable     :: phi_restart(:,:)         !< Restart solution vector
// ! real,   allocatable     :: phi_total_restart(:,:)   !< Restart total solution vector
// ! real,   allocatable     :: phi_dt_restart(:,:)      !< Restart temporal derivative solution vector

// ! integer                 :: n_old_stages = 0         !< Number of old stages found in file
// ! real,   allocatable     :: phi_old_restart(:,:,:)   !< Restart solution array (old timesteps for implicit)

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void restart_define();

#endif /* RESTART_MODULE_H */