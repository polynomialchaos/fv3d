/*******************************************************************************
 * @file navier_stokes_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef NAVIER_STOKES_MODULE_H
#define NAVIER_STOKES_MODULE_H

#include "fv3d/fv3d_module.h"



extern int navier_stokes_active;

extern const double RM;
extern const double molar_mass_air;

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

extern int ic_rho;
extern int ic_rho_u;
extern int ic_rho_v;
extern int ic_rho_w;
extern int ic_rho_e;
extern int ip_u;
extern int ip_v;
extern int ip_w;
extern int ip_p;
extern int ip_T;

enum BoundaryType
{
    BoundaryFlow,
    BoundaryInflow,
    BoundaryOutflow,
    BoundaryAdiabaticWall,
    BoundaryIsothermalWall,
    BoundarySlipWall,
    BoundarySymmetry,
    BoundaryState,
    BoundaryFunction,
    BoundaryTypeMax
};

void navier_stokes_define();
void boundary_define();
void flux_define();

void update_boundaries(double t);
void update_gradients_boundaries();

void prim_to_con(double *phi);
void copy_prim_to_con(double *phi_i, double *phi_j);
void con_to_prim(double *phi);
void copy_con_to_prim(double *phi_i, double *phi_j);
double calc_ig_p(double rho, double T, double R_mix);
double calc_ig_rho(double p, double T, double R_mix);
double calc_ig_T(double p, double rho, double R_mix);
double calc_riemann_p(double *phi);

void calc_flux();

#endif /* NAVIER_STOKES_MODULE_H */