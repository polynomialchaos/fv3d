//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef NAVIER_STOKES_MODULE_H
#define NAVIER_STOKES_MODULE_H

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
const double RM;
const double molar_mass_air;

extern double cfl_scale;
extern double dfl_scale;
extern double mu_mix;
extern double R_mix;
extern double Pr;
extern double kappa;
extern double kappa_m1;
extern double kappa_p1;
extern double s_kappa;
extern double s_kappa_m1;
extern double s_kappa_p1;
extern double cp;
extern double cv;
extern double kappa_pr;
extern double lambda;

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

    // enum, bind( C )
    //     enumerator :: BND_NULL, BND_FLOW, BND_INFLOW, BND_OUTFLOW, BND_AD_WALL, BND_IT_WALL, &
    //         BND_SLIP_WALL, BND_SYMM, BND_STATE, BND_FUNC
    // end enum

    // character(len=_STRLEN_),    parameter   :: bnd_type_names(BND_FUNC) = [character(len=_STRLEN_)  :: &
    //     'flow', 'inflow', 'outflow', 'wall-adiabatic', 'wall-isothermal', 'wall-slip', 'symmetry', 'state', 'function']
    // character(len=_STRLEN_) :: flux_scheme = 'AUSM'             !< The Riemann solver

    // procedure(),    pointer :: flux_scheme_routine => null()    !< Riemann solver routine

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void navier_stokes_define();
void boundary_define();
void flux_define();

void update_boundaries( double t );
void update_gradients_boundaries( double t );

void calc_flux();

#endif /* NAVIER_STOKES_MODULE_H */