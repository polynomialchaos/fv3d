/*******************************************************************************
 * @file timedisc.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "analyze/analyze_private.h"
#include "output/output_private.h"
#include "restart/restart_private.h"
#include "timedisc_private.h"

void_timestep_ft timestep_function_pointer = NULL;
double_calc_timestep_ft calc_timestep_function_pointer = NULL;

bool_t solver_is_explicit = BC_TRUE;

bool_t solver_is_viscous_dt = 0;
int solver_max_iter = 1000;
bool_t solver_is_transient = 1;
double solver_abort_residual = 1e-10;
double solver_t_start = 0.0;
double solver_t_end = 0.2;

/*******************************************************************************
 * @brief Free timedisc
 ******************************************************************************/
void free_timedisc()
{
    timestep_function_pointer = NULL;
    calc_timestep_function_pointer = NULL;
}

/*******************************************************************************
 * @brief Initialize timedisc
 * @param timedisc_type
 * @param max_iter
 * @param is_transient
 * @param abort_residual
 * @param t_start
 * @param t_end
 ******************************************************************************/
void init_timedisc(timedisc_type_t timedisc_type, int max_iter,
                   bool_t is_transient, double abort_residual,
                   double t_start, double t_end)
{
    switch (timedisc_type)
    {
    case Explicit:
        solver_is_explicit = BC_TRUE;
        break;
    case Implicit:
        solver_is_explicit = BC_FALSE;
        break;
    default:
        BM_CHECK_EXPRESSION(0);
        break;
    }

    solver_max_iter = max_iter;
    solver_is_transient = is_transient;
    solver_abort_residual = abort_residual;

    solver_t_start = t_start;
    solver_t_end = solver_is_transient == BC_TRUE ? t_end : BC_DMAX;
}

/*******************************************************************************
 * @brief Return is_explicit flag
 * bool_t
 ******************************************************************************/
bool_t is_explicit()
{
    return solver_is_explicit;
}

/*******************************************************************************
 * @brief Print the residual line
 * @param iter
 * @param t
 * @param dt
 * @param do_output
 ******************************************************************************/
void print_residual(int iter, double t, double dt, bool_t do_output)
{
    char output_str = (do_output == BC_TRUE) ? '*' : ' ';
    char viscous_str = (solver_is_viscous_dt == BC_TRUE) ? 'T' : 'F';

    if (solver_is_explicit)
    {
        BM_PRINT("%09d %12.5e %12.5e %c %c:", iter, t, dt, viscous_str, output_str);
    }
    else
    {
        BM_PRINT("%09d %12.5e %12.5e %c %c %6d %6d:", iter, t, dt, viscous_str, output_str, n_iter_inner, n_iter_lsoe);
    }

    for (int i = 0; i < solver_variables->n_sol_variables; ++i)
        BM_PRINT(" %12.5e", solver_residual[i]);

    BM_PRINT("\n");
}

/*******************************************************************************
 * @brief Print the residual line header
 ******************************************************************************/
void print_residual_header()
{
    if (solver_is_explicit)
    {
        BM_PRINT("%9s %12s %12s %1s %1s:", "iter", "time", "dt", "V", "O");
    }
    else
    {
        BM_PRINT("%9s %12s %12s %1s %1s %6s %6s:", "iter", "time", "dt", "V", "O", "inner", "lsoe");
    }

    for (int i = 0; i < solver_variables->n_sol_variables; ++i)
        BM_PRINT(" %12s", solver_variables->sol_variables[i].name);

    BM_PRINT("\n");
}

/*******************************************************************************
 * @brief Time discretizazion routine
 ******************************************************************************/
void timedisc()
{
    /* check for errors */
    check_abort(0);

    bool_t do_output = BC_FALSE;
    bool_t do_finalize = BC_FALSE;

    int iter = 0;
    if (solver_use_restart == BC_TRUE)
        iter = solver_iter_restart;

    double t = solver_t_start;
    if (solver_use_restart == BC_TRUE)
        t = solver_t_restart;

    double dt;

    print_residual_header();
    if (iter == 0)
        finite_volume_time_derivative(t);
    if ((solver_do_output_data == BC_TRUE) && (solver_use_restart == BC_FALSE))
        write_output(iter, t);

    while (1)
    {
        /* check for errors */
        check_abort(0);

        /* calculate time step dt */
        double dt_local = calc_timestep_function_pointer(t);
        BM_MPI_ALL_REDUCE(MPIDouble, MPIMin, &dt_local, &dt);

        if (t + dt > solver_t_end)
        {
            dt = solver_t_end - t;
            do_output = BC_TRUE;
            do_finalize = BC_TRUE;
        }

        /* the timestep to be called */
        timestep_function_pointer(iter, t, dt);
        calc_global_residual(dt);

        /* check for NAN and INF */
        if (is_nan_n(solver_residual, solver_variables->n_sol_variables) ||
            is_inf_n(solver_residual, solver_variables->n_sol_variables))
            BM_CHECK_EXPRESSION(0);

        t = t + dt;
        iter = iter + 1;

        /* steady-state simulation */
        if ((solver_is_transient == BC_FALSE) && (max_n(solver_residual, solver_variables->n_sol_variables) < solver_abort_residual))
        {
            solver_t_end = t;
            do_output = BC_TRUE;
            do_finalize = BC_TRUE;
        }

        if ((solver_i_output_data > 0) && (iter % solver_i_output_data == 0))
        {
            do_output = BC_TRUE;
        }

        /* maximum iteration number reacher */
        if (iter >= solver_iter_restart + solver_max_iter)
        {
            solver_t_end = t;
            do_output = BC_TRUE;
            do_finalize = BC_TRUE;
        }

        print_residual(iter, t, dt, do_output);

        /* output */
        if ((solver_do_output_data == BC_TRUE) && (do_output == BC_TRUE))
        {
            write_output(iter, t);
            do_output = BC_FALSE;
        }

        /* stop if required */
        if (do_finalize)
            break;
    }
}

/*******************************************************************************
 * @brief Set the timestep calculation routine
 * @param fun_ptr
 ******************************************************************************/
void set_calc_timestep(double_calc_timestep_ft fun_ptr)
{
    calc_timestep_function_pointer = fun_ptr;
}

/*******************************************************************************
 * @brief Set the viscous timestep flag
 * @param is_viscous_dt
 ******************************************************************************/
void set_viscous_dt(bool_t is_viscous_dt)
{
    solver_is_viscous_dt = is_viscous_dt;
}