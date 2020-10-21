//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "navier_stokes_private.h"
#include "equation/navier_stokes/boundary_private.h"
#include "equation/navier_stokes/flux_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int navier_stokes_active = 0;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void navier_stokes_initialize();
void navier_stokes_finalize();

const double RM = 8.31446261815324;
const double molar_mass_air = 28.96e-3;

double cfl_scale = 1.00;
double dfl_scale = 1.00;
double mu_mix = 18.13e-6;
double R_mix = RM / molar_mass_air;
double Pr = 0.718;
double kappa = 1.4;

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void navier_stokes_define()
{
    register_initialize_routine( navier_stokes_initialize );
    register_finalize_routine( navier_stokes_finalize );

    set_parameter( "Equation/Navier-Stokes/cfl_scale", ParameterNumber, &cfl_scale, "The cfl_loc scale factor (0 ... 1)", NULL, 0 );
    set_parameter( "Equation/Navier-Stokes/dfl_scale", ParameterNumber, &dfl_scale, "The DFL scale factor (0 ... 1)", NULL, 0 );
    set_parameter( "Equation/Navier-Stokes/mu_mix", ParameterNumber, &mu_mix, "The dynamic viscosity, N s m-2", NULL, 0 );
    set_parameter( "Equation/Navier-Stokes/R_mix", ParameterNumber, &R_mix, "The specific gas constant, J kg-1 K-1", NULL, 0 );
    set_parameter( "Equation/Navier-Stokes/Pr", ParameterNumber, &Pr, "The Prandtl number", NULL, 0 );
    set_parameter( "Equation/Navier-Stokes/kappa", ParameterNumber, &kappa, "The isentropic exponent", NULL, 0 );

    boundary_define();
    flux_define();
}

void navier_stokes_initialize()
{
    if (navier_stokes_active == 0) return;

    get_parameter( "Equation/Navier-Stokes/cfl_scale", ParameterNumber, &cfl_scale );
    get_parameter( "Equation/Navier-Stokes/dfl_scale", ParameterNumber, &dfl_scale );
    get_parameter( "Equation/Navier-Stokes/mu_mix", ParameterNumber, &mu_mix );
    get_parameter( "Equation/Navier-Stokes/R_mix", ParameterNumber, &R_mix );
    get_parameter( "Equation/Navier-Stokes/Pr", ParameterNumber, &Pr );
    get_parameter( "Equation/Navier-Stokes/kappa", ParameterNumber, &kappa );

        // IC_RHO      = add_variable( 'rho', is_active=.true. )
        // IC_RHO_U    = add_variable( 'rho_u', is_active=.true. )
        // IC_RHO_V    = add_variable( 'rho_v', is_active=.true. )
        // IC_RHO_W    = add_variable( 'rho_w', is_active=.true. )
        // IC_RHO_E    = add_variable( 'rho_e', is_active=.true. )
        // IC_RHO_UVW  = [IC_RHO_U, IC_RHO_V, IC_RHO_W]
        // IC          = [IC_RHO, IC_RHO_U, IC_RHO_V, IC_RHO_W, IC_RHO_E]

        // IP_U        = add_variable( 'u (prim)', is_active=.false. )
        // IP_V        = add_variable( 'v (prim)', is_active=.false. )
        // IP_W        = add_variable( 'w (prim)', is_active=.false. )
        // IP_P        = add_variable( 'p (prim)', is_active=.false. )
        // IP_T        = add_variable( 'T (prim)', is_active=.false. )
        // IP_UVW      = [IP_U, IP_V, IP_W]
        // IP          = [IP_U, IP_V, IP_W, IP_P, IP_T]

        // kappa_m1    = kappa - 1.0
        // kappa_p1    = kappa + 1.0
        // s_kappa     = 1.0 / kappa
        // s_kappa_m1  = 1.0 / kappa_m1
        // s_kappa_p1  = 1.0 / kappa_p1

        // cv          = R_mix * s_kappa_m1
        // cp          = kappa * cv
        // kappa_pr    = kappa / Pr
        // lambda      = mu_mix * cp / Pr

        // exact_func_routine => exact_func
        // update_routine => update
        // update_gradients_routine => update_bc_gradients

        // calc_time_step_routine => calc_time_step
        // calc_flux_routine => calc_flux

}

void navier_stokes_finalize()
{
}

void exact_func( int id, double t, double *x, double *phi )
{

        // select case ( id )
        //     case ( 0 )
        //         phi = regions(flow_region)%phi
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown initial function selected (got: ' // set_string( id ) // ')!' )
        // end select

}

void update( double t )
{

        // do i = 1, n_cells
        //     phi_total(:,i) = con_2_prim( phi_total(:,i) )
        // end do

        // call update_bc( t )

}

void calc_time_step( double dt_min )
{
        // lambda_conv = huge( lambda_conv )
        // lambda_visc = huge( lambda_visc )

        // do i = 1, n_cells
        //     s_rho   = 1 / phi(IC_RHO,i)
        //     c       = sqrt( kappa * phi_total(IP_P,i) * s_rho ) ! speed of

        //     ds          = dot_product( abs( phi_total(IP_UVW,i) ) + c, cells(i)%ds )
        //     lambda_conv = min( lambda_conv, cells(i)%volume / ds )

        //     ds          = dot_product( cells(i)%ds, cells(i)%ds )
        //     lambda_visc = min( lambda_visc, s_rho * ds / cells(i)%volume )
        // end do

        // dt_loc(1)   = cfl_scale * lambda_conv
        // dt_loc(2)   = dfl_scale * mu_mix * kappa_pr * lambda_visc

        // is_viscous_dt = (dt_loc(2) .lt. dt_loc(1))
        // dt_min = minval( dt_loc )
}