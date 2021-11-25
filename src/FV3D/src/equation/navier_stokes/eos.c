//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <math.h>
#include "navier_stokes_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"





void prim_to_con(double *phi)
{
    phi[ic_rho_u] = phi[ip_u] * phi[ic_rho]; // rho * u
    phi[ic_rho_v] = phi[ip_v] * phi[ic_rho]; // rho * v
    phi[ic_rho_w] = phi[ip_w] * phi[ic_rho]; // rho * w
    phi[ic_rho_e] = s_kappa_m1 * phi[ip_p] +
                    0.5 * dot_n(&phi[ip_u], &phi[ic_rho_u], DIM); // rho * e
}

void copy_prim_to_con(double *phi_i, double *phi_j)
{
    copy_n(phi_i, all_variables->n_tot_variables, phi_j);
    prim_to_con(phi_j);
}

void con_to_prim(double *phi)
{
    phi[ip_u] = phi[ic_rho_u] / phi[ic_rho]; // u
    phi[ip_v] = phi[ic_rho_v] / phi[ic_rho]; // v
    phi[ip_w] = phi[ic_rho_w] / phi[ic_rho]; // w
    phi[ip_p] = kappa_m1 * (phi[ic_rho_e] -
                            0.5 * dot_n(&phi[ic_rho_u], &phi[ip_u], DIM)); // p
    phi[ip_p] = MAX(1.00E-10, phi[ip_p]);                                // pressure must not be negative
    phi[ip_T] = calc_ig_T(phi[ip_p], phi[ic_rho], R_mix);                  // T
}

void copy_con_to_prim(double *phi_i, double *phi_j)
{
    copy_n(phi_i, all_variables->n_tot_variables, phi_j);
    con_to_prim(phi_j);
}

double calc_ig_p(double rho, double T, double R_mix)
{
    return rho * R_mix * T;
}

double calc_ig_rho(double p, double T, double R_mix)
{
    return p / (R_mix * T);
}

double calc_ig_T(double p, double rho, double R_mix)
{
    return p / (rho * R_mix);
}

double calc_riemann_p(double *phi)
{
    if (phi[ip_u] <= 0.0)
    {
        double c = sqrt(kappa * phi[ip_p] / phi[ic_rho]);
        return phi[ip_p] * pow(
                               MAX(1e-4, 1. + 0.5 * kappa_m1 * phi[ip_u] / c),
                               2 * kappa * s_kappa_m1);
    }
    else
    {
        double ar = 2 * s_kappa_p1 / phi[ip_u];
        double br = kappa_m1 * s_kappa_p1 * phi[ip_p];
        return phi[ip_p] + phi[ip_u] / ar * 0.5 * (phi[ip_u] + sqrt(phi[ip_u] * phi[ip_u] + 4. * ar * (phi[ip_p] + br)));
    }
}