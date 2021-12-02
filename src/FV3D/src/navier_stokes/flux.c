/*******************************************************************************
 * @file flux.c
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

string_t flux_scheme_name = NULL;

typedef void (*void_calc_conv_flux_ft)(double *phi_l, double *phi_r, double *f);
void_calc_conv_flux_ft calc_convective_flux_function_pointer = NULL;

/*******************************************************************************
 * @brief Calculate the flux
 ******************************************************************************/
void calc_flux()
{
    Faces_t *faces = solver_mesh->faces;
    Boundaries_t *boundaries = solver_mesh->boundaries;
    Regions_t *regions = solver_mesh->regions;
    int n_internal_faces = faces->n_internal_faces;
    int n_boundary_faces = faces->n_boundary_faces;
    int n_tot_variables = solver_variables->n_tot_variables;
    int n_sol_variables = solver_variables->n_sol_variables;

    double phi_total_left_r[n_tot_variables];
    double phi_total_right_r[n_tot_variables];
    double flux_c[n_sol_variables];
    double flux_d_x[n_sol_variables], flux_d_y[n_sol_variables], flux_d_z[n_sol_variables];

    for (int ii = 0; ii < n_internal_faces; ++ii)
    {
        int i = faces->internal_faces[ii];
        int *fc = &faces->cells[i * FACE_CELLS];
        double *flux_i = &flux[i * n_sol_variables];
        double *n = &faces->n[i * DIM];
        double *t1 = &faces->t1[i * DIM];
        double *t2 = &faces->t2[i * DIM];

        double *phi_total_left_i = &phi_total_left[i * n_tot_variables];
        double *phi_total_right_i = &phi_total_right[i * n_tot_variables];

        phi_total_left_r[ic_rho] = phi_total_left_i[ic_rho];
        phi_total_left_r[ic_rho_u] = dot_n(&phi_total_left_i[ic_rho_u], n, DIM);
        phi_total_left_r[ic_rho_v] = dot_n(&phi_total_left_i[ic_rho_u], t1, DIM);
        phi_total_left_r[ic_rho_w] = dot_n(&phi_total_left_i[ic_rho_u], t2, DIM);
        phi_total_left_r[ic_rho_e] = phi_total_left_i[ic_rho_e];
        con_to_prim(phi_total_left_r);

        phi_total_right_r[ic_rho] = phi_total_right_i[ic_rho];
        phi_total_right_r[ic_rho_u] = dot_n(&phi_total_right_i[ic_rho_u], n, DIM);
        phi_total_right_r[ic_rho_v] = dot_n(&phi_total_right_i[ic_rho_u], t1, DIM);
        phi_total_right_r[ic_rho_w] = dot_n(&phi_total_right_i[ic_rho_u], t2, DIM);
        phi_total_right_r[ic_rho_e] = phi_total_right_i[ic_rho_e];
        con_to_prim(phi_total_right_r);

        /* calculate convective flux */
        calc_convective_flux_function_pointer(phi_total_left_r, phi_total_right_r, flux_c);

        /* rotate convective flux back to global coordiantes and add to flux */
        flux_i[0] = flux_c[0];
        flux_i[1] = flux_c[1] * n[0] + flux_c[2] * t1[0] + flux_c[3] * t2[0];
        flux_i[2] = flux_c[1] * n[1] + flux_c[2] * t1[1] + flux_c[3] * t2[1];
        flux_i[3] = flux_c[1] * n[2] + flux_c[2] * t1[2] + flux_c[3] * t2[2];
        flux_i[4] = flux_c[4];

        double *grad_phi_total_x_i_0 = &grad_phi_total_x[fc[0] * n_tot_variables];
        double *grad_phi_total_y_i_0 = &grad_phi_total_y[fc[0] * n_tot_variables];
        double *grad_phi_total_z_i_0 = &grad_phi_total_z[fc[0] * n_tot_variables];
        double *grad_phi_total_x_i_1 = &grad_phi_total_x[fc[1] * n_tot_variables];
        double *grad_phi_total_y_i_1 = &grad_phi_total_y[fc[1] * n_tot_variables];
        double *grad_phi_total_z_i_1 = &grad_phi_total_z[fc[1] * n_tot_variables];

        /* calculate diffusive flux */
        viscous_flux(
            phi_total_left_i, grad_phi_total_x_i_0, grad_phi_total_y_i_0, grad_phi_total_z_i_0,
            phi_total_right_i, grad_phi_total_x_i_1, grad_phi_total_y_i_1, grad_phi_total_z_i_1,
            flux_d_x, flux_d_y, flux_d_z);

        /* add diffusive flux to flux */
        for (int j = 0; j < n_sol_variables; ++j)
            flux_i[j] += flux_d_x[j] * n[0] + flux_d_y[j] * n[1] + flux_d_z[j] * n[2];
    }

    for (int ii = 0; ii < n_boundary_faces; ++ii)
    {
        int i = faces->boundary_faces[ii];
        int *fc = &faces->cells[i * FACE_CELLS];
        int id = boundaries->id[faces->boundary[i]];
        double *flux_i = &flux[i * n_sol_variables];
        double *n = &faces->n[i * DIM];
        double *t1 = &faces->t1[i * DIM];
        double *t2 = &faces->t2[i * DIM];

        double *phi_total_left_i = &phi_total_left[i * n_tot_variables];
        double *phi_total_right_i = &phi_total_right[i * n_tot_variables];

        phi_total_left_r[ic_rho] = phi_total_left_i[ic_rho];
        phi_total_left_r[ic_rho_u] = dot_n(&phi_total_left_i[ic_rho_u], n, DIM);
        phi_total_left_r[ic_rho_v] = dot_n(&phi_total_left_i[ic_rho_u], t1, DIM);
        phi_total_left_r[ic_rho_w] = dot_n(&phi_total_left_i[ic_rho_u], t2, DIM);
        phi_total_left_r[ic_rho_e] = phi_total_left_i[ic_rho_e];
        con_to_prim(phi_total_left_r);

        phi_total_right_r[ic_rho] = phi_total_right_i[ic_rho];
        phi_total_right_r[ic_rho_u] = dot_n(&phi_total_right_i[ic_rho_u], n, DIM);
        phi_total_right_r[ic_rho_v] = dot_n(&phi_total_right_i[ic_rho_u], t1, DIM);
        phi_total_right_r[ic_rho_w] = dot_n(&phi_total_right_i[ic_rho_u], t2, DIM);
        phi_total_right_r[ic_rho_e] = phi_total_right_i[ic_rho_e];
        con_to_prim(phi_total_right_r);

        /* calculate convective flux */
        switch (regions->type[id])
        {
        default:
            calc_convective_flux_function_pointer(phi_total_left_r, phi_total_right_r, flux_c);
            break;
        }

        /* rotate convective flux back to global coordiantes and add to flux */
        flux_i[0] = flux_c[0];
        flux_i[1] = flux_c[1] * n[0] + flux_c[2] * t1[0] + flux_c[3] * t2[0];
        flux_i[2] = flux_c[1] * n[1] + flux_c[2] * t1[1] + flux_c[3] * t2[1];
        flux_i[3] = flux_c[1] * n[2] + flux_c[2] * t1[2] + flux_c[3] * t2[2];
        flux_i[4] = flux_c[4];

        double *grad_phi_total_x_i_0 = &grad_phi_total_x[fc[0] * n_tot_variables];
        double *grad_phi_total_y_i_0 = &grad_phi_total_y[fc[0] * n_tot_variables];
        double *grad_phi_total_z_i_0 = &grad_phi_total_z[fc[0] * n_tot_variables];
        double *grad_phi_total_x_i_1 = &grad_phi_total_x[fc[1] * n_tot_variables];
        double *grad_phi_total_y_i_1 = &grad_phi_total_y[fc[1] * n_tot_variables];
        double *grad_phi_total_z_i_1 = &grad_phi_total_z[fc[1] * n_tot_variables];

        /* calculate diffusive flux */
        switch (regions->type[id])
        {
        default:
            viscous_flux(
                phi_total_left_i, grad_phi_total_x_i_0, grad_phi_total_y_i_0, grad_phi_total_z_i_0,
                phi_total_right_i, grad_phi_total_x_i_1, grad_phi_total_y_i_1, grad_phi_total_z_i_1,
                flux_d_x, flux_d_y, flux_d_z);
            break;
        }

        /* add diffusive flux to flux */
        for (int j = 0; j < n_sol_variables; ++j)
            flux_i[j] += flux_d_x[j] * n[0] + flux_d_y[j] * n[1] + flux_d_z[j] * n[2];
    }
}

/*******************************************************************************
 * @brief Calculate the Euler 1D convective flux
 ******************************************************************************/
void eval_euler_flux_1d(double *phi, double *f)
{
    f[0] = phi[ic_rho_u];                           /* rho * u */
    f[1] = phi[ic_rho_u] * phi[ip_u] + phi[ip_p];   /* rho * u * u + p */
    f[2] = phi[ic_rho_u] * phi[ip_v];               /* rho * u * v */
    f[3] = phi[ic_rho_u] * phi[ip_w];               /* rho * u * w */
    f[4] = (phi[ic_rho_e] + phi[ip_p]) * phi[ip_u]; /* (rho * e + p) * u */
}

/*******************************************************************************
 * @brief Calculate the Euler 1D viscous flux
 ******************************************************************************/
void eval_viscous_flux_1d(double *phi, double *grad_phi_x, double *grad_phi_y, double *grad_phi_z,
                          double *f, double *g, double *h)
{
    const double s_23 = 2.0 / 3.0;
    const double s_43 = 4.0 / 3.0;

    double tau_xx = mu_mix * (s_43 * grad_phi_x[1] -
                              s_23 * grad_phi_y[2] - s_23 * grad_phi_z[3]); /* 4/3 * mu * u_x - 2/3 * mu * v_y - 2/3 * mu * w_z */
    double tau_yy = mu_mix * (-s_23 * grad_phi_x[1] +
                              s_43 * grad_phi_y[2] - s_23 * grad_phi_z[3]); /* -2/3 * mu * u_x + 4/3 * mu * v_y - 2/3 * mu * w_z */
    double tau_zz = mu_mix * (-s_23 * grad_phi_x[1] -
                              s_23 * grad_phi_y[2] + s_43 * grad_phi_z[3]); /* -2/3 * mu * u_x - 2/3 * mu * v_y + 4/3 * mu * w_z */

    double tau_xy = mu_mix * (grad_phi_y[1] + grad_phi_x[2]); /* mu * (u_y + v_x) */
    double tau_xz = mu_mix * (grad_phi_z[1] + grad_phi_x[3]); /* mu * (u_z + w_x) */
    double tau_yz = mu_mix * (grad_phi_z[2] + grad_phi_y[3]); /* mu * (y_z + w_y) */

    f[0] = 0.0;
    f[1] = -tau_xx; /* -4/3 * mu * u_x + 2/3 * mu * (v_y + w_z) */
    f[2] = -tau_xy; /* -mu * (u_y + v_x) */
    f[3] = -tau_xz; /* -mu * (u_z + w_x) */
    f[4] = -tau_xx * phi[ip_u] - tau_xy * phi[ip_v] -
           tau_xz * phi[ip_w] - lambda * grad_phi_x[ip_T]; /* -(tau_xx * phi + tau_xy * v + tau_xz * w - q_x) q_x=-lambda * T_x */

    g[0] = 0.0;
    g[1] = -tau_xy; /* -mu * (u_y + v_x) */
    g[2] = -tau_yy; /* -4/3 * mu * v_y + 2/3 * mu * (u_x + w_z) */
    g[3] = -tau_yz; /* -mu * (y_z + w_y) */
    g[4] = -tau_xy * phi[ip_u] - tau_yy * phi[ip_v] -
           tau_yz * phi[ip_w] - lambda * grad_phi_y[ip_T]; /* -(tau_yx * phi + tau_yy * v + tau_yz * w - q_y) q_y=-lambda * T_y */

    h[0] = 0.0;
    h[1] = -tau_xz; /* -mu * (u_z + w_x) */
    h[2] = -tau_yz; /* -mu * (y_z + w_y) */
    h[3] = -tau_zz; /* -4/3 * mu * w_z + 2/3 * mu * (u_x + v_y) */
    h[4] = -tau_xz * phi[ip_u] - tau_yz * phi[ip_v] -
           tau_zz * phi[ip_w] - lambda * grad_phi_z[ip_T]; /* -(tau_zx * phi + tau_zy * v + tau_zz * w - q_z) q_z=-lambda * T_z */
}

/*******************************************************************************
 * @brief Define flux
 ******************************************************************************/
void flux_define()
{
    REGISTER_INITIALIZE_ROUTINE(flux_initialize);
    REGISTER_FINALIZE_ROUTINE(flux_finalize);

    string_t tmp_opt[] = {"AUSM", "Rusanov"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    SET_PARAMETER("Equation/Navier-Stokes/Flux/flux_scheme", StringParameter, &tmp,
                  "The Riemann solver", &tmp_opt, tmp_opt_n);
}

/*******************************************************************************
 * @brief Finalize flux
 ******************************************************************************/
void flux_finalize()
{
    DEALLOCATE(flux_scheme_name);

    calc_convective_flux_function_pointer = NULL;
}

/*******************************************************************************
 * @brief Initialize flux
 ******************************************************************************/
void flux_initialize()
{
    GET_PARAMETER("Equation/Navier-Stokes/Flux/flux_scheme", StringParameter, &flux_scheme_name);

    if (is_equal(flux_scheme_name, "Rusanov"))
    {
        calc_convective_flux_function_pointer = riemann_rusanonv;
    }
    else if (is_equal(flux_scheme_name, "AUSM"))
    {
        calc_convective_flux_function_pointer = riemann_ausm;
    }
    else
    {
        CHECK_EXPRESSION(0);
    }
}

/*******************************************************************************
 * @brief Riemann solver from Rusanov
 * @param phi_l
 * @param phi_r
 * @param f
 ******************************************************************************/
void riemann_rusanonv(double *phi_l, double *phi_r, double *f)
{
    int n_sol_variables = solver_variables->n_sol_variables;

    double c_l = sqrt(kappa * phi_l[ip_p] / phi_l[ic_rho]);
    double c_r = sqrt(kappa * phi_r[ip_p] / phi_r[ic_rho]);
    double eigval = MAX(ABS(phi_l[ip_u]) + c_l, ABS(phi_r[ip_u]) + c_r);

    double f_l[n_sol_variables];
    double f_r[n_sol_variables];
    eval_euler_flux_1d(phi_l, f_l);
    eval_euler_flux_1d(phi_r, f_r);

    for (int j = 0; j < n_sol_variables; ++j)
        f[j] = 0.5 * (f_l[j] + f_r[j]) - 0.5 * eigval * (phi_r[j] - phi_l[j]);
}

/*******************************************************************************
 * @brief Riemann solver from AUSM
 * @param phi_l
 * @param phi_r
 * @param f
 ******************************************************************************/
void riemann_ausm(double *phi_l, double *phi_r, double *f)
{
    double c_l = sqrt(kappa * phi_l[ip_p] / phi_l[ic_rho]);
    double c_r = sqrt(kappa * phi_r[ip_p] / phi_r[ic_rho]);

    double M_l = phi_l[ip_u] / c_l;
    double M_r = phi_r[ip_u] / c_r;

    double H_l = (phi_l[ic_rho_e] + phi_l[ip_p]) / phi_l[ic_rho];
    double H_r = (phi_r[ic_rho_e] + phi_r[ip_p]) / phi_r[ic_rho];

    double M_m, M_p, P_m, P_p;

    /* Positive M and p in the LEFT cell */
    if (M_l <= -1.0)
    {
        M_p = 0.0;
        P_p = 0.0;
    }
    else if (M_l < 1.0)
    {
        M_p = 0.25 * (M_l + 1.0) * (M_l + 1.0);
        P_p = 0.25 * phi_l[ip_p] * (1.0 + M_l) * (1.0 + M_l) * (2.0 - M_l); /* or use P_p = half*(1.0+M_l)*phi_l[ip_v]; */
    }
    else
    {
        M_p = M_l;
        P_p = phi_l[ip_p];
    }

    /* Negative M and p in the RIGHT cell */
    if (M_r <= -1.0)
    {
        M_m = M_r;
        P_m = phi_r[ip_p];
    }
    else if (M_r < 1.0)
    {
        M_m = -0.25 * (M_r - 1.0) * (M_r - 1.0);
        P_m = 0.25 * phi_r[ip_p] * (1.0 - M_r) * (1.0 - M_r) * (2.0 + M_r); /* or use P_m = half*(1.0-M_r)*phi_r[ip_v]; */
    }
    else
    {
        M_m = 0.0;
        P_m = 0.0;
    }

    /* Positive Part of Flux evaluated in the left cell. */
    f[0] = MAX(0.0, M_p + M_m) * c_l * phi_l[ic_rho];
    f[1] = MAX(0.0, M_p + M_m) * c_l * phi_l[ic_rho] * phi_l[ip_u] + P_p;
    f[2] = MAX(0.0, M_p + M_m) * c_l * phi_l[ic_rho] * phi_l[ip_v];
    f[3] = MAX(0.0, M_p + M_m) * c_l * phi_l[ic_rho] * phi_l[ip_w];
    f[4] = MAX(0.0, M_p + M_m) * c_l * phi_l[ic_rho] * H_l;

    /* Negative Part of Flux evaluated in the right cell. */
    f[0] += MIN(0.0, M_p + M_m) * c_r * phi_r[ic_rho];
    f[1] += MIN(0.0, M_p + M_m) * c_r * phi_r[ic_rho] * phi_r[ip_u] + P_m;
    f[2] += MIN(0.0, M_p + M_m) * c_r * phi_r[ic_rho] * phi_r[ip_v];
    f[3] += MIN(0.0, M_p + M_m) * c_r * phi_r[ic_rho] * phi_r[ip_w];
    f[4] += MIN(0.0, M_p + M_m) * c_r * phi_r[ic_rho] * H_r;
}

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
                  double *phi_r, double *grad_phi_x_r, double *grad_phi_y_r, double *grad_phi_z_r, double *f, double *g, double *h)
{
    int n_sol_variables = solver_variables->n_sol_variables;

    double f_l[n_sol_variables], g_l[n_sol_variables], h_l[n_sol_variables];
    double f_r[n_sol_variables], g_r[n_sol_variables], h_r[n_sol_variables];
    eval_viscous_flux_1d(phi_l, grad_phi_x_l, grad_phi_y_l, grad_phi_z_l, f_l, g_l, h_l);
    eval_viscous_flux_1d(phi_r, grad_phi_x_r, grad_phi_y_r, grad_phi_z_r, f_r, g_r, h_r);

    for (int j = 0; j < n_sol_variables; ++j)
    {
        f[j] = 0.5 * (f_l[j] + f_r[j]);
        g[j] = 0.5 * (g_l[j] + g_r[j]);
        h[j] = 0.5 * (h_l[j] + h_r[j]);
    }
}
