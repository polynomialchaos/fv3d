//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef NAVIER_STOKES_PRIVATE_H
#define NAVIER_STOKES_PRIVATE_H

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
extern int navier_stokes_active;
    // real,   parameter               :: RM = 8.31446261815324            !< Universal gas constant, J mol-1 K-1
    // real,   parameter               :: molar_mass_air = 28.96e-3        !< Molar mass of dry air, kg mol-1

    // real                            :: cfl_scale = 1.00                 !< The CFL scale factor
    // real                            :: dfl_scale = 1.00                 !< The DFL scale factor
    // real                            :: mu_mix = 18.13e-6                !< The dynamic viscosity, N s m-2
    // real                            :: R_mix = RM / molar_mass_air      !< The specific gas constant, J kg-1 K-1
    // real                            :: Pr = 0.718                       !< The Prandtl number
    // real                            :: kappa = 1.4                      !< The isentropic exponent

    //

    // real                            :: kappa_m1, kappa_p1               !< kappa-1, kappa+1
    // real                            :: s_kappa, s_kappa_m1, s_kappa_p1  !< 1/kappa, 1/(kappa-1), 1/(kappa+1)
    // real                            :: cp, cv, kappa_pr, lambda         !< cp, cv, kappa/Pr, lambda

    // integer,    parameter           :: n_con            = 5
    // integer                         :: IC(n_con)        = -1
    // integer                         :: IC_RHO           = -1
    // integer                         :: IC_RHO_U         = -1
    // integer                         :: IC_RHO_V         = -1
    // integer                         :: IC_RHO_W         = -1
    // integer                         :: IC_RHO_E         = -1
    // integer                         :: IC_RHO_UVW(3)    = -1

    // integer,    parameter           :: n_prim       = 5
    // integer                         :: IP(n_prim)   = -1
    // integer                         :: IP_U         = -1
    // integer                         :: IP_V         = -1
    // integer                         :: IP_W         = -1
    // integer                         :: IP_P         = -1
    // integer                         :: IP_T         = -1
    // integer                         :: IP_UVW(3)    = -1

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void navier_stokes_define();

#endif /* NAVIER_STOKES_PRIVATE_H */