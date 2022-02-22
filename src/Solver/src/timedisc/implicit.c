
/*******************************************************************************
 * @file implicit.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include <math.h>
#include "timedisc_private.h"

int solver_max_iter_inner = 100;
double solver_tolerance_lsoe = 1e-12;
int solver_max_iter_lsoe = 100;
int solver_max_krylov_dims = 15;
int solver_max_krylov_restarts = 2;

int n_iter_inner = 0;
int n_iter_lsoe = 0;
double **phi_old = NULL;

bool_t is_bicgstab = BTRU;
int n_work_size = 0;
double dt_loc = 0.0;
double tpdt_loc = 0.0;

enum
{
    NBDFStagesEuler = 1,
    NBDFStagesBDF2 = 2
};

double bdf_a_euler[NBDFStagesEuler] = {0.0};
double bdf_b_euler = 1.0;

double bdf_a_bdf2[NBDFStagesBDF2] = {-1. / 3., 1 / 3.};
double bdf_b_bdf2 = 2. / 3.;

int n_bdf_stages = 0;
double *bdf_a = NULL;
double bdf_b = 1.0;

int n_bdf_stages_loc = 0;
double *bdf_a_loc = NULL;
double bdf_b_loc = 1.0;

double *work = NULL;

double *Y_n = NULL;
double *f_Y_n = NULL;
double *dY_n = NULL;
double *dY_dt_n = NULL;
double *jac = NULL;

/*******************************************************************************
 * @brief Numerical jacobian routine
 * @param n_var
 * @param n_cells
 ******************************************************************************/
void calc_jacobian_numerical(int n_var, int n_cells)
{
    int n_tot_variables = solver_variables->n_tot_variables;

    for (int i_var = 0; i_var < n_var; ++i_var)
    {
        double eps_fd = 0.0;
        for (int i = 0; i < n_cells; ++i)
            eps_fd += Y_n[i * n_var + i_var] * Y_n[i * n_var + i_var];

        eps_fd = sqrt(eps_fd) * 1e-4;

        /* positive + eps */
        for (int i = 0; i < n_cells; ++i)
        {
            for (int j = 0; j < n_var; ++j)
                solver_data->phi_total[i * n_tot_variables + j] = Y_n[i * n_var + j];

            solver_data->phi_total[i * n_tot_variables + i_var] += 0.5 * eps_fd;
        }

        finite_volume_time_derivative(tpdt_loc);

        for (int i = 0; i < n_cells; ++i)
        {
            int idx_i = i * n_var * n_var + i_var * n_var;

            for (int j = 0; j < n_var; ++j)
            {
                jac[idx_i + j] = -solver_data->phi_dt[i * n_var + j];
            }
        }

        /* negative + eps */
        for (int i = 0; i < n_cells; ++i)
        {
            for (int j = 0; j < n_var; ++j)
                solver_data->phi_total[i * n_tot_variables + j] = Y_n[i * n_var + j];

            solver_data->phi_total[i * n_tot_variables + i_var] -= 0.5 * eps_fd;
        }

        finite_volume_time_derivative(tpdt_loc);

        for (int i = 0; i < n_cells; ++i)
        {
            int idx_i = i * n_var * n_var + i_var * n_var;

            for (int j = 0; j < n_var; ++j)
            {
                jac[idx_i + j] += solver_data->phi_dt[i * n_var + j];
                jac[idx_i + j] *= bdf_b_loc / (eps_fd + SMALL);
            }

            jac[idx_i + i_var] += 1.0 / dt_loc;
        }
    }
}

/*******************************************************************************
 * @brief Free implicit timedisc
 ******************************************************************************/
void free_implicit()
{
    for (int i = 0; i < n_bdf_stages; ++i)
        DEALLOCATE(phi_old[i]);
    DEALLOCATE(phi_old);

    DEALLOCATE(work);

    DEALLOCATE(Y_n);
    DEALLOCATE(f_Y_n);
    DEALLOCATE(dY_n);
    DEALLOCATE(dY_dt_n);
    DEALLOCATE(jac);
}

/*******************************************************************************
 * @brief Initialize implicit timedisc
 * @param implicit_scheme
 * @param implicit_solver
 * @param max_iter_inner
 * @param tolerance_lsoe
 * @param max_iter_lsoe
 * @param max_krylov_dims
 * @param max_krylov_restarts
 ******************************************************************************/
void init_implicit(implicit_scheme_t implicit_scheme,
                   implicit_solver_t implicit_solver,
                   int max_iter_inner,
                   double tolerance_lsoe, int max_iter_lsoe,
                   int max_krylov_dims, int max_krylov_restarts)
{
    Cells_t *cells = solver_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_sol_variables = solver_variables->n_sol_variables;

    timestep_function_pointer = timestep_newton;

    solver_max_iter_inner = max_iter_inner;
    solver_tolerance_lsoe = tolerance_lsoe;
    solver_max_iter_lsoe = max_iter_lsoe;
    solver_max_krylov_dims = max_krylov_dims;
    solver_max_krylov_restarts = max_krylov_restarts;

    switch (implicit_scheme)
    {
    case EulerImplicit:
        n_bdf_stages = NBDFStagesEuler;
        bdf_a = bdf_a_euler;
        bdf_b = bdf_b_euler;
        break;
    case BDF2:
        n_bdf_stages = NBDFStagesBDF2;
        bdf_a = bdf_a_bdf2;
        bdf_b = bdf_b_bdf2;
        break;
    default:
        CHECK_EXPRESSION(0);
        break;
    }

    switch (implicit_solver)
    {
    case BiCGStab:
        is_bicgstab = BTRU;
        n_work_size = get_bicgstab_n_m_work_size(n_sol_variables, n_domain_cells);
        work = ALLOCATE(sizeof(double) * n_work_size);
        break;
    case GMRes:
        is_bicgstab = BFLS;
        n_work_size = get_gmres_n_m_work_size(n_sol_variables, n_domain_cells, solver_max_krylov_dims);
        work = ALLOCATE(sizeof(double) * n_work_size);
        break;
    default:
        CHECK_EXPRESSION(0);
        break;
    }

    phi_old = ALLOCATE(sizeof(double *) * n_bdf_stages);
    for (int i = 0; i < n_bdf_stages; ++i)
        phi_old[i] = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);

    Y_n = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);
    f_Y_n = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);
    dY_n = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);
    dY_dt_n = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);
    jac = ALLOCATE(sizeof(double) * n_sol_variables * n_sol_variables * n_domain_cells);
}

/*******************************************************************************
 * @brief Implicit time discretizazion routine (Newton)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void timestep_newton(int iter, double t, double dt)
{
    Cells_t *cells = solver_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_sol_variables = solver_variables->n_sol_variables;
    int n_tot_variables = solver_variables->n_tot_variables;

    /* discretization parameters */
    n_bdf_stages_loc = (iter < n_bdf_stages) ? NBDFStagesEuler : n_bdf_stages;
    bdf_a_loc = (iter < n_bdf_stages) ? bdf_a_euler : bdf_a;
    bdf_b_loc = (iter < n_bdf_stages) ? bdf_b_euler : bdf_b;

    /* set local timestep */
    dt_loc = dt;
    tpdt_loc = t + dt_loc;

    /* store the old state before calculating the FVTimeDerivative */
    for (int i = n_bdf_stages_loc - 1; i > 0; --i)
        copy_n(phi_old[i - 1], n_sol_variables * n_domain_cells, phi_old[i]);

    for (int i = 0; i < n_domain_cells; ++i)
        for (int j = 0; j < n_sol_variables; ++j)
            phi_old[0][i * n_sol_variables + j] = solver_data->phi_total[i * n_tot_variables + j];

    /* fill inital values for newton iteration */
    copy_n(phi_old[0], n_sol_variables * n_domain_cells, Y_n);
    copy_n(solver_data->phi_dt, n_sol_variables * n_domain_cells, dY_dt_n);

    /* calculate the inital error for newton abort criterion */
    for (int i = 0; i < n_domain_cells; ++i)
        for (int j = 0; j < n_sol_variables; ++j)
        {
            int idx = i * n_sol_variables + j;
            f_Y_n[idx] = -(Y_n[idx] - phi_old[0][idx]) / dt_loc + bdf_b_loc * dY_dt_n[idx];
        }

    for (int i_stage = 0; i_stage < n_bdf_stages_loc; ++i_stage)
        for (int i = 0; i < n_domain_cells; ++i)
            for (int j = 0; j < n_sol_variables; ++j)
            {
                int idx = i * n_sol_variables + j;
                f_Y_n[idx] -= bdf_a_loc[i_stage] / dt_loc * phi_old[i_stage][idx];
            }

    const double err_f_Y_0 = len_n(f_Y_n, n_sol_variables * n_domain_cells);
    double err_f_Y_old = err_f_Y_0;

    for (n_iter_inner = 1; n_iter_inner <= solver_max_iter_inner; ++n_iter_inner)
    {
        calc_jacobian_numerical(n_sol_variables, n_domain_cells);

        if (err_f_Y_old >= err_f_Y_0)
        {
            /* Jac * dY = fY_n => dY ... Jacobian is determined via finite difference. fY_n = phi - dt * RHS */
            n_iter_lsoe = solver_max_iter_lsoe;
            double residual_lsoe = solver_tolerance_lsoe * err_f_Y_old;
            if (is_bicgstab)
            {
                solve_bicgstab_n_m(n_sol_variables, n_domain_cells, f_Y_n, dY_n,
                                   work, matrix_vector_numerical, &n_iter_lsoe, &residual_lsoe);
            }
            else
            {
                solve_gmres_n_m(n_sol_variables, n_domain_cells, f_Y_n, dY_n,
                                work, matrix_vector_numerical, &n_iter_lsoe, &residual_lsoe, solver_max_krylov_dims, solver_max_krylov_restarts);
            }

            /* Y^(n+1) = Y^(n) + (Y^(n+1)-Y^(n)) */
            for (int i = 0; i < n_domain_cells; ++i)
                for (int j = 0; j < n_sol_variables; ++j)
                {
                    int idx = i * n_sol_variables + j;
                    Y_n[idx] += dY_n[idx];
                    solver_data->phi_total[i * n_tot_variables + j] = Y_n[idx];
                }

            finite_volume_time_derivative(tpdt_loc);
            copy_n(solver_data->phi_dt, n_sol_variables * n_domain_cells, dY_dt_n);

            for (int i = 0; i < n_domain_cells; ++i)
                for (int j = 0; j < n_sol_variables; ++j)
                {
                    int idx = i * n_sol_variables + j;
                    f_Y_n[idx] = -(Y_n[idx] - phi_old[0][idx]) / dt_loc + bdf_b_loc * dY_dt_n[idx];
                }

            for (int i_stage = 0; i_stage < n_bdf_stages_loc; ++i_stage)
                for (int i = 0; i < n_domain_cells; ++i)
                    for (int j = 0; j < n_sol_variables; ++j)
                    {
                        int idx = i * n_sol_variables + j;
                        f_Y_n[idx] -= bdf_a_loc[i_stage] / dt_loc * phi_old[i_stage][idx];
                    }

            err_f_Y_old = len_n(f_Y_n, n_sol_variables * n_domain_cells);
        }
        else
        {
            finite_volume_time_derivative(tpdt_loc);
        }

        int is_error_less = (err_f_Y_old < err_f_Y_0);
        int is_error_less_g;
        MPI_ALL_REDUCE(MPIInt, MPILogAnd, &is_error_less, &is_error_less_g);

        if (is_error_less_g == 1)
            break;
        if (solver_is_transient == BFLS)
            break;
    }

    if (n_iter_inner >= solver_max_iter_inner)
        CHECK_EXPRESSION(0);
}

/*******************************************************************************
 * @brief Matrix vector routine (called by solver)
 * @param x
 * @param b
 * @param n_var
 * @param n_cells
 * @return int
 ******************************************************************************/
int matrix_vector_numerical(double *x, double *b, size_t n_var, size_t n_cells)
{
    for (size_t i = 0; i < n_cells; ++i)
    {
        size_t idx_i = i * n_var * n_var;
        double *b_i = &b[i * n_var];
        double *x_i = &x[i * n_var];

        for (size_t j = 0; j < n_var; ++j)
        {
            b_i[j] = jac[idx_i + j] * x_i[0];
        }

        for (size_t i_var = 1; i_var < n_var; ++i_var)
        {
            for (size_t j = 0; j < n_var; ++j)
            {
                b_i[j] += jac[idx_i + i_var * n_var + j] * x_i[i_var];
            }
        }
    }

    return 0;
}