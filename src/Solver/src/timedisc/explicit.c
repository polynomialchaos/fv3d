/*******************************************************************************
 * @file explicit.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "timedisc_private.h"

enum
{
    NRKStagesEuler = 1,
    NRKStagesRK33 = 3,
    NRKStagesRK45 = 5
};

double rk_a_euler[NRKStagesEuler] = {0.0};
double rk_b_euler[NRKStagesEuler] = {0.0};
double rk_g_euler[NRKStagesEuler] = {1.0};

double rk_a_rk33[NRKStagesRK33] = {0., -5. / 9., -153. / 128.};
double rk_b_rk33[NRKStagesRK33] = {0., 1. / 3., 3. / 4.};
double rk_g_rk33[NRKStagesRK33] = {1. / 3., 15. / 16., 8. / 15.};

double rk_a_rk45[NRKStagesRK45] = {0., -567301805773.0 / 1357537059087.0,
                                   -2404267990393.0 / 2016746695238.0, -3550918686646.0 / 2091501179385.0,
                                   -1275806237668.0 / 842570457699.0};
double rk_b_rk45[NRKStagesRK45] = {0., 1432997174477.0 / 9575080441755.0,
                                   2526269341429.0 / 6820363962896.0, 2006345519317.0 / 3224310063776.0,
                                   2802321613138.0 / 2924317926251.0};
double rk_g_rk45[NRKStagesRK45] = {1432997174477.0 / 9575080441755.0,
                                   5161836677717.0 / 13612068292357.0, 1720146321549.0 / 2090206949498.0,
                                   3134564353537.0 / 4481467310338.0, 2277821191437.0 / 14882151754819.0};

int n_rk_stages = 0;
double *rk_a = NULL;
double *rk_b = NULL;
double *rk_g = NULL;

/*******************************************************************************
 * @brief Free explicit
 ******************************************************************************/
void free_explicit()
{
    rk_a = NULL;
    rk_b = NULL;
    rk_g = NULL;
}

/*******************************************************************************
 * @brief Initialize explicit
 * @param explicit_scheme
 ******************************************************************************/
void init_explicit(explicit_scheme_t explicit_scheme)
{
    time_step_function_pointer = time_step_lserkw2;

    switch (explicit_scheme)
    {
    case EulerExplicit:
        n_rk_stages = NRKStagesEuler;
        rk_a = rk_a_euler;
        rk_b = rk_b_euler;
        rk_g = rk_g_euler;
        break;
    case RungeKutta33:
        n_rk_stages = NRKStagesRK33;
        rk_a = rk_a_rk33;
        rk_b = rk_b_rk33;
        rk_g = rk_g_rk33;
        break;
    case RungeKutta45:
        n_rk_stages = NRKStagesRK45;
        rk_a = rk_a_rk45;
        rk_b = rk_b_rk45;
        rk_g = rk_g_rk45;
        break;
    default:
        CHECK_EXPRESSION(0);
        break;
    }
}

/*******************************************************************************
 * @brief Explicit time discretizazion routine (LSERKW2)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void time_step_lserkw2(int iter, double t, double dt)
{
#if DEBUG
    UNUSED(iter);
#endif /* DEBUG */

    Cells_t *cells = solver_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_sol_variables = solver_variables->n_sol_variables;
    int n_tot_variables = solver_variables->n_tot_variables;

    double phi_dt_tmp[n_sol_variables * n_domain_cells];

    /* first stage */
    double t_stage = t; /* + dt * rk_b[i_stage] = 0 */
    finite_volume_time_derivative(t_stage);

    for (int i = 0; i < n_domain_cells; ++i)
        for (int j = 0; j < n_sol_variables; ++j)
            phi_dt_tmp[i * n_sol_variables + j] = solver_phi_dt[i * n_sol_variables + j]; /* + phi_dt_tmp * rk_a[i_stage] = 0 */

    for (int i = 0; i < n_domain_cells; ++i)
        for (int j = 0; j < n_sol_variables; ++j)
            solver_phi_total[i * n_tot_variables + j] += phi_dt_tmp[i * n_sol_variables + j] * dt * rk_g[0];

    /* 2nd to n_rk_stages */
    for (int i_stage = 1; i_stage < n_rk_stages; ++i_stage)
    {
        t_stage = t + dt * rk_b[i_stage];
        finite_volume_time_derivative(t_stage);

        for (int i = 0; i < n_domain_cells; ++i)
            for (int j = 0; j < n_sol_variables; ++j)
                phi_dt_tmp[i * n_sol_variables + j] = solver_phi_dt[i * n_sol_variables + j] + phi_dt_tmp[i * n_sol_variables + j] * rk_a[i_stage];

        for (int i = 0; i < n_domain_cells; ++i)
            for (int j = 0; j < n_sol_variables; ++j)
                solver_phi_total[i * n_tot_variables + j] += phi_dt_tmp[i * n_sol_variables + j] * dt * rk_g[i_stage];
    }
}