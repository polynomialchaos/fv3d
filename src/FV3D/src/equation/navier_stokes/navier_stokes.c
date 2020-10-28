//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <math.h>
#include "navier_stokes_module.h"
#include "mesh/mesh_module.h"
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
double calc_time_step();

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
    Cells_t *cells      = global_mesh->cells;
    int n_domain_cells  = cells->n_domain_cells;
    int n_tot_variables = all_variables->n_tot_variables;

    for ( int i = 0; i < n_domain_cells; i++ )
        con_to_prim( &phi_total[i*n_tot_variables] );

    update_boundaries( t );
}

void update_gradients()
{
    update_gradients_boundaries();
}

void calc_exact_func( int id, double t, double *x, double *phi )
{
#ifdef DEBUG
    u_unused( t );
    u_unused( x );
#endif /* DEBUG */
    Regions_t *regions  = global_mesh->regions;
    int n_tot_variables = all_variables->n_tot_variables;

    switch (id)
    {
        case BoundaryFlow:
            copy_n( &regions->phi_total[regions->flow_region*n_tot_variables], phi, n_tot_variables );
            break;
        default:
            check_error( 0 );
            break;
    }
}

double calc_time_step()
{
    Cells_t *cells      = global_mesh->cells;
    int n_cells         = cells->n_domain_cells;
    int n_tot_variables = all_variables->n_tot_variables;

    double lambda_conv = DOUBLE_MAX;
    double lambda_visc = DOUBLE_MAX;

    for ( int i = 0; i < n_cells; i++ )
    {
        double *phi_total_i = &phi_total[i*n_tot_variables];
        double *dx          = &cells->dx[i*DIM];
        double s_rho        = 1.0 / phi_total_i[ic_rho];
        double c            = sqrt( kappa * phi_total_i[ip_p] * s_rho );

        double ds1  = (u_abs( phi_total_i[ip_u] ) + c) * dx[0] +
            (u_abs( phi_total_i[ip_v] ) + c) * dx[1] + (u_abs( phi_total_i[ip_w] ) + c) * dx[2];
        lambda_conv = u_min( lambda_conv, cells->volume[i] / ds1 );

        double ds2  = dot_n( dx, dx, DIM );
        lambda_visc = u_min( lambda_visc, s_rho * ds2 / cells->volume[i] );
    }

    lambda_conv *= cfl_scale;
    lambda_visc *= dfl_scale * mu_mix * kappa_pr * lambda_visc;

    is_viscous_dt   = (lambda_visc < lambda_conv);

    return u_min( lambda_conv, lambda_visc );
}