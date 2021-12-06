/*******************************************************************************
 * @file navier_stokes.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <math.h>
#include "navier_stokes_module.h"
#include "mesh/mesh_module.h"
#include "fv/fv_module.h"
#include "timedisc/timedisc_module.h"

const double RM = 8.31446261815324;
const double molar_mass_air = 28.96e-3;

double cfl_scale = 1.00;
double dfl_scale = 1.00;
double mu_mix = 18.13e-6;
double R_mix = RM / molar_mass_air;
double Pr = 0.718;
double kappa = 1.4;

double kappa_m1;
double kappa_p1;
double s_kappa;
double s_kappa_m1;
double s_kappa_p1;
double cv;
double cp;
double kappa_pr;
double lambda;

int ic_rho = -1;
int ic_rho_u = -1;
int ic_rho_v = -1;
int ic_rho_w = -1;
int ic_rho_e = -1;
int ip_u = -1;
int ip_v = -1;
int ip_w = -1;
int ip_p = -1;
int ip_T = -1;

/*******************************************************************************
 * @brief Define navier_stokes
 ******************************************************************************/
void navier_stokes_define()
{
    REGISTER_INITIALIZE_ROUTINE(navier_stokes_initialize);
    REGISTER_FINALIZE_ROUTINE(navier_stokes_finalize);

    SET_PARAMETER("Equation/Navier-Stokes/cfl_scale", NumberParameter,
                  &cfl_scale,
                  "The cfl_loc scale factor (0 ... 1)", NULL, 0);
    SET_PARAMETER("Equation/Navier-Stokes/dfl_scale", NumberParameter,
                  &dfl_scale,
                  "The DFL scale factor (0 ... 1)", NULL, 0);
    SET_PARAMETER("Equation/Navier-Stokes/mu_mix", NumberParameter, &mu_mix,
                  "The dynamic viscosity, N s m-2", NULL, 0);
    SET_PARAMETER("Equation/Navier-Stokes/R_mix", NumberParameter, &R_mix,
                  "The specific gas constant, J kg-1 K-1", NULL, 0);
    SET_PARAMETER("Equation/Navier-Stokes/Pr", NumberParameter, &Pr,
                  "The Prandtl number", NULL, 0);
    SET_PARAMETER("Equation/Navier-Stokes/kappa", NumberParameter, &kappa,
                  "The isentropic exponent", NULL, 0);

    boundary_define();
    flux_define();
}

/*******************************************************************************
 * @brief Finalize navier_stokes
 ******************************************************************************/
void navier_stokes_finalize()
{
    free_equation();
}

/*******************************************************************************
 * @brief Initialize navier_stokes
 ******************************************************************************/
void navier_stokes_initialize()
{
    GET_PARAMETER("Equation/Navier-Stokes/cfl_scale", NumberParameter,
                  &cfl_scale);
    GET_PARAMETER("Equation/Navier-Stokes/dfl_scale", NumberParameter,
                  &dfl_scale);
    GET_PARAMETER("Equation/Navier-Stokes/mu_mix", NumberParameter, &mu_mix);
    GET_PARAMETER("Equation/Navier-Stokes/R_mix", NumberParameter, &R_mix);
    GET_PARAMETER("Equation/Navier-Stokes/Pr", NumberParameter, &Pr);
    GET_PARAMETER("Equation/Navier-Stokes/kappa", NumberParameter, &kappa);

    kappa_m1 = kappa - 1.0;
    kappa_p1 = kappa + 1.0;
    s_kappa = 1.0 / kappa;
    s_kappa_m1 = 1.0 / kappa_m1;
    s_kappa_p1 = 1.0 / kappa_p1;
    cv = R_mix * s_kappa_m1;
    cp = kappa * cv;
    kappa_pr = kappa / Pr;
    lambda = mu_mix * cp / Pr;

    init_equation();
    ic_rho = add_sol_variable("rho");
    ic_rho_u = add_sol_variable("rho_u");
    ic_rho_v = add_sol_variable("rho_v");
    ic_rho_w = add_sol_variable("rho_w");
    ic_rho_e = add_sol_variable("rho_e");

    ip_u = add_dep_variable("u");
    ip_v = add_dep_variable("v");
    ip_w = add_dep_variable("w");
    ip_p = add_dep_variable("p");
    ip_T = add_dep_variable("T");

    set_calc_flux(calc_ns_flux);
    set_calc_timestep(calc_ns_timestep);
    set_exact_function(calc_ns_exact_function);
    set_update(ns_update);
    set_update_gradients(ns_update_gradients);
}

/*******************************************************************************
 * @brief Define solution by region ID
 * @param id
 * @param t
 * @param x
 * @param phi
 ******************************************************************************/
void calc_ns_exact_function(int id, double t, double *x, double *phi)
{
#ifdef DEBUG
    UNUSED(t);
    UNUSED(x);
#endif /* DEBUG */
    Regions_t *regions = solver_mesh->regions;
    int n_tot_variables = solver_variables->n_tot_variables;

    switch (id)
    {
    case BoundaryFlow:
        copy_n(&regions->phi_total[regions->flow_region * n_tot_variables], n_tot_variables, phi);
        break;
    default:
        CHECK_EXPRESSION(0);
        break;
    }
}

/*******************************************************************************
 * @brief Calculate the timestep
 * @param time
 * @return double
 ******************************************************************************/
double calc_ns_timestep(double time)
{
#ifdef DEBUG
    UNUSED(time);
#endif /* DEBUG */
    Cells_t *cells = solver_mesh->cells;
    int n_cells = cells->n_domain_cells;
    int n_tot_variables = solver_variables->n_tot_variables;

    double lambda_conv = BDMX;
    double lambda_visc = BDMX;

    for (int i = 0; i < n_cells; ++i)
    {
        double *phi_total_i = &solver_data->phi_total[i * n_tot_variables];
        double *dx = &cells->dx[i * DIM];
        double s_rho = 1.0 / phi_total_i[ic_rho];
        double c = sqrt(kappa * phi_total_i[ip_p] * s_rho);

        double ds1 = (ABS(phi_total_i[ip_u]) + c) * dx[0] +
                     (ABS(phi_total_i[ip_v]) + c) * dx[1] + (ABS(phi_total_i[ip_w]) + c) * dx[2];
        lambda_conv = MIN(lambda_conv, cells->volume[i] / ds1);

        double ds2 = dot_n(dx, dx, DIM);
        lambda_visc = MIN(lambda_visc, s_rho * ds2 / cells->volume[i]);
    }

    lambda_conv *= cfl_scale;
    lambda_visc *= dfl_scale * mu_mix * kappa_pr * lambda_visc;

    set_viscous_dt(lambda_visc < lambda_conv);

    return MIN(lambda_conv, lambda_visc);
}

/*******************************************************************************
 * @brief Update solution for internal data and boundaries
 * @param time
 ******************************************************************************/
void ns_update(double time)
{
    Cells_t *cells = solver_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_tot_variables = solver_variables->n_tot_variables;

    for (int i = 0; i < n_domain_cells; ++i)
        con_to_prim(&solver_data->phi_total[i * n_tot_variables]);

    update_boundaries(time);
}

/*******************************************************************************
 * @brief Update gradeints for internal data and boundaries
 * @param time
 ******************************************************************************/
void ns_update_gradients(double time)
{
    update_gradients_boundaries(time);
}