/*******************************************************************************
 * @file navier_stokes_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
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

/*******************************************************************************
 * @brief Define boundary
 ******************************************************************************/
void boundary_define();

/*******************************************************************************
 * @brief Finalize boundary
 ******************************************************************************/
void boundary_finalize();

/*******************************************************************************
 * @brief Initialize boundary
 ******************************************************************************/
void boundary_initialize();

/*******************************************************************************
 * @brief Define solution by region ID
 * @param id
 * @param t
 * @param x
 * @param phi
 ******************************************************************************/
void calc_ns_exact_function(int id, double t, double *x, double *phi);

/*******************************************************************************
 * @brief Calculate pressure (ideal gas)
 * @param rho
 * @param T
 * @param R_mix
 * @return double
 ******************************************************************************/
double calc_ig_p(double rho, double T, double R_mix);

/*******************************************************************************
 * @brief Calculate density (ideal gas)
 * @param p
 * @param T
 * @param R_mix
 * @return double
 ******************************************************************************/
double calc_ig_rho(double p, double T, double R_mix);

/*******************************************************************************
 * @brief Calculate temperature (ideal gas)
 * @param p
 * @param rho
 * @param R_mix
 * @return double
 ******************************************************************************/
double calc_ig_T(double p, double rho, double R_mix);

/*******************************************************************************
 * @brief Calculate Riemann pressure
 * @param phi
 * @return double
 ******************************************************************************/
double calc_riemann_p(double *phi);

/*******************************************************************************
 * @brief Calculate the flux
 * @param time
 ******************************************************************************/
void calc_ns_flux(double time);

/*******************************************************************************
 * @brief Calculate the timestep
 * @param time
 * @return double
 ******************************************************************************/
double calc_ns_timestep(double time);

/*******************************************************************************
 * @brief Convert conservative to primitive variables
 * @param phi
 ******************************************************************************/
void con_to_prim(double *phi);

/*******************************************************************************
 * @brief Copy and convert conservative to primitive variables
 * @param phi_i
 * @param phi_j
 ******************************************************************************/
void copy_con_to_prim(double *phi_i, double *phi_j);

/*******************************************************************************
 * @brief Copy and convert primitive to conservative variables
 * @param phi_i
 * @param phi_j
 ******************************************************************************/
void copy_prim_to_con(double *phi_i, double *phi_j);

/*******************************************************************************
 * @brief Calculate the Euler 1D convective flux
 ******************************************************************************/
void eval_euler_flux_1d(double *phi, double *f);

/*******************************************************************************
 * @brief Calculate the Euler 1D viscous flux
 ******************************************************************************/
void eval_viscous_flux_1d(double *phi, double *grad_phi_x, double *grad_phi_y, double *grad_phi_z,
                          double *f, double *g, double *h);

/*******************************************************************************
 * @brief Define flux
 ******************************************************************************/
void flux_define();

/*******************************************************************************
 * @brief Finalize flux
 ******************************************************************************/
void flux_finalize();

/*******************************************************************************
 * @brief Initialize flux
 ******************************************************************************/
void flux_initialize();

/*******************************************************************************
 * @brief Define navier_stokes
 ******************************************************************************/
void navier_stokes_define();

/*******************************************************************************
 * @brief Finalize navier_stokes
 ******************************************************************************/
void navier_stokes_finalize();

/*******************************************************************************
 * @brief Initialize navier_stokes
 ******************************************************************************/
void navier_stokes_initialize();

/*******************************************************************************
 * @brief Parse primitive state parameter
 ******************************************************************************/
void parse_primitive_state(cstring_t prefix, double *phi);

/*******************************************************************************
 * @brief Convert primitive to conservative variables
 * @param phi
 ******************************************************************************/
void prim_to_con(double *phi);

/*******************************************************************************
 * @brief Riemann solver from Rusanov
 * @param phi_l
 * @param phi_r
 * @param f
 ******************************************************************************/
void riemann_rusanonv(double *phi_l, double *phi_r, double *f);

/*******************************************************************************
 * @brief Riemann solver from AUSM
 * @param phi_l
 * @param phi_r
 * @param f
 ******************************************************************************/
void riemann_ausm(double *phi_l, double *phi_r, double *f);

/*******************************************************************************
 * @brief Update boundaries
 ******************************************************************************/
void update_boundaries(double t);

/*******************************************************************************
 * @brief Update gradient boundaries
 * @param time
 ******************************************************************************/
void update_gradients_boundaries(double time);

/*******************************************************************************
 * @brief Update solution for internal data and boundaries
 * @param time
 ******************************************************************************/
void ns_update(double time);

/*******************************************************************************
 * @brief Update gradeints for internal data and boundaries
 * @param time
 ******************************************************************************/
void ns_update_gradients(double time);

/*******************************************************************************
 * @brief Viscous flux
 * @param phi_l
 * @param grad_phi_x_l
 * @param grad_phi_y_l
 * @param grad_phi_z_l
 * @param phi_r
 * @param grad_phi_x_r
 * @param grad_phi_y_r
 * @param grad_phi_z_r
 * @param f
 * @param g
 * @param h
 ******************************************************************************/
void viscous_flux(double *phi_l, double *grad_phi_x_l, double *grad_phi_y_l, double *grad_phi_z_l,
                  double *phi_r, double *grad_phi_x_r, double *grad_phi_y_r, double *grad_phi_z_r, double *f, double *g, double *h);

#endif /* NAVIER_STOKES_MODULE_H */