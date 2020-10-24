//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "navier_stokes_module.h"
#include "equation/equation_module.h"
#include "fv/fv_module.h"
#include "timedisc/timedisc_module.h"

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

void update( double t );
void update_gradients();
void calc_exact_func( int id, double t, double *x, double *phi );
void calc_time_step( double dt_min );

const double RM = 8.31446261815324;
const double molar_mass_air = 28.96e-3;

double cfl_scale    = 1.00;
double dfl_scale    = 1.00;
double mu_mix       = 18.13e-6;
double R_mix        = RM / molar_mass_air;
double Pr           = 0.718;
double kappa        = 1.4;

double kappa_m1;
double kappa_p1;
double s_kappa;
double s_kappa_m1;
double s_kappa_p1;
double cv;
double cp;
double kappa_pr;
double lambda;

const int n_cons    = 5;
const int n_prins   = 5;
int ic_rho          = -1;
int ic_rho_u        = -1;
int ic_rho_v        = -1;
int ic_rho_w        = -1;
int ic_rho_e        = -1;
int ip_u            = -1;
int ip_v            = -1;
int ip_w            = -1;
int ip_p            = -1;
int ip_T            = -1;

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

    kappa_m1    = kappa - 1.0;
    kappa_p1    = kappa + 1.0;
    s_kappa     = 1.0 / kappa;
    s_kappa_m1  = 1.0 / kappa_m1;
    s_kappa_p1  = 1.0 / kappa_p1;
    cv          = R_mix * s_kappa_m1;
    cp          = kappa * cv;
    kappa_pr    = kappa / Pr;
    lambda      = mu_mix * cp / Pr;

    ic_rho      = add_sol_variable( all_variables, "rho" );
    ic_rho_u    = add_sol_variable( all_variables, "rho_u" );
    ic_rho_v    = add_sol_variable( all_variables, "rho_v" );
    ic_rho_w    = add_sol_variable( all_variables, "rho_w" );
    ic_rho_e    = add_sol_variable( all_variables, "rho_e" );

    ip_u        = add_dep_variable( all_variables, "u" );
    ip_v        = add_dep_variable( all_variables, "v" );
    ip_w        = add_dep_variable( all_variables, "w" );
    ip_p        = add_dep_variable( all_variables, "p" );
    ip_T        = add_dep_variable( all_variables, "T" );

    update_function_pointer             = update;
    update_gradients_function_pointer   = update_gradients;
    calc_exact_function_pointer         = calc_exact_func;
    calc_time_step_function_pointer     = calc_time_step;
    calc_flux_function_pointer          = calc_flux;
}

void navier_stokes_finalize()
{
}

void update( double t )
{

        // do i = 1, n_cells
        //     phi_total(:,i) = con_to_prim( phi_total(:,i) )
        // end do

        // call update_bc( t )

}

void update_gradients()
{

        // do i = 1, n_cells
        //     phi_total(:,i) = con_to_prim( phi_total(:,i) )
        // end do

        // call update_bc( t )

}

void calc_exact_func( int id, double t, double *x, double *phi )
{

        // select case ( id )
        //     case ( 0 )
        //         phi = regions(flow_region)%phi
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown initial function selected (got: ' // set_string( id ) // ')!' )
        // end select

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