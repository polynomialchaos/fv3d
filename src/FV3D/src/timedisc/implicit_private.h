//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef IMPLICIT_PRIVATE_H
#define IMPLICIT_PRIVATE_H

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
extern int implicit_active;

    // character(len=_STRLEN_)     :: scheme = 'BDF-2'                 !< The implicit scheme
    // character(len=_STRLEN_)     :: method = 'Newton'                !< The method to solve the non-linear system of equations
    // integer                     :: max_iter = 100                   !< The maximum number of inner iterations
    // character(len=_STRLEN_)     :: jacobian_type = 'Numerical'      !< The type of jacobian generation
    // character(len=_STRLEN_)     :: solver = 'GMRes'                 !< The linear system of equations solver
    // real                        :: tolerance_lsoe = 1e-6            !< The linear solver tolerance
    // integer                     :: max_iter_lsoe = 100              !< The linear solver maximum number of iterations
    // integer                     :: max_krylov_dims = 10             !< The maximum Krylov space dimension
    // integer                     :: max_krylov_restarts = 3          !< The maximum Krylov space restarts

    // procedure(GMRes_d), pointer :: solve_lsoe => null()             !< Linear system of equations solver

    // integer                     :: n_iter_inner                     !< Global value of Newton iterations
    // integer                     :: n_iter_lsoe                      !< Global value of LSoE iterations

    // real,   allocatable         :: phi_old(:,:,:)                   !< The old solution vector

    // integer                     :: n_stages = 0                     !< Number of time update stages
    // real,   allocatable         :: bdf_a(:)                         !< Values of the bdf scheme (aStage)
    // real                        :: bdf_b = 0.0                      !< Values of the bdf scheme (bStage)

    // real,   allocatable         :: Y_n(:,:)                         !< Solver solution array
    // real,   allocatable         :: Y_dt_n(:,:)                      !< Solver resiudal array

    // real,   parameter           :: bdf_a_euler(*) = [0.0]
    // real,   parameter           :: bdf_b_euler = 1.0

    // real,   parameter           :: bdf_a_bdf2(*) = [-1./3.,1/3.]
    // real,   parameter           :: bdf_b_bdf2 = 2./3.

    // real                        :: loc_t                            !< Local time
    // real                        :: loc_dt                           !< Local timestep
    // real,   allocatable         :: loc_bdf_a(:)                     !< Local aStage factor
    // real                        :: loc_bdf_b                        !< Local bStage factor

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void implicit_define();
void time_step_newton( double t, double dt, int iter );

#endif /* IMPLICIT_PRIVATE_H */