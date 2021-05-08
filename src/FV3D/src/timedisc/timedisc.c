
//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "timedisc_module.h"
#include "explicit/explicit_module.h"
#include "implicit/implicit_module.h"
#include "equation/equation_module.h"
#include "restart/restart_module.h"
#include "analyze/analyze_module.h"
#include "output/output_module.h"
#include "fv/fv_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
void_timestep_fp_t time_step_function_pointer = NULL;
double_calc_timestep_fp_t calc_time_step_function_pointer = NULL;

int is_viscous_dt = 0;

string_t time_step_name = NULL;
int max_iter = 1000;
int is_transient = 1;
double abort_residual = 1e-10;
double t_start = 0.0;
double t_end = 0.2;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void timedisc_initialize();
void timedisc_finalize();

void print_residual_header();
void print_residual(int iter, double t, double dt, int do_output);

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void timedisc_define()
{
    register_initialize_routine(timedisc_initialize);
    register_finalize_routine(timedisc_finalize);

    string_t tmp_opt[] = {"Explicit", "Implicit"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    set_parameter("TimeDisc/time_step", ParameterString, &tmp, "The timestep mehtod", &tmp_opt, tmp_opt_n);
    set_parameter("TimeDisc/max_iter", ParameterDigit, &max_iter, "The maximum number of iterations", NULL, 0);
    set_parameter("TimeDisc/transient", ParameterBool, &is_transient, "The flag wheter to be transient or steady-state", NULL, 0);
    set_parameter("TimeDisc/abort_residual", ParameterNumber, &abort_residual, "The abort residual", NULL, 0);
    set_parameter("TimeDisc/t_start", ParameterNumber, &t_start, "The start time", NULL, 0);
    set_parameter("TimeDisc/t_end", ParameterNumber, &t_end, "The end time", NULL, 0);

    explicit_define();
    implicit_define();
}

void timedisc_initialize()
{
    get_parameter("TimeDisc/time_step", ParameterString, &time_step_name);
    get_parameter("TimeDisc/max_iter", ParameterDigit, &max_iter);
    get_parameter("TimeDisc/transient", ParameterBool, &is_transient);
    get_parameter("TimeDisc/abort_residual", ParameterNumber, &abort_residual);
    get_parameter("TimeDisc/t_start", ParameterNumber, &t_start);
    get_parameter("TimeDisc/t_end", ParameterNumber, &t_end);

    if (is_transient == 0)
        t_end = DOUBLE_MAX;

    if (is_equal(time_step_name, "Explicit"))
    {
        explicit_active = 1;
    }
    else if (is_equal(time_step_name, "Implicit"))
    {
        implicit_active = 1;
    }
    else
    {
        check_error(0);
    }
}

void timedisc_finalize()
{
    time_step_function_pointer = NULL;
    calc_time_step_function_pointer = NULL;

    deallocate(time_step_name);
}

void timedisc()
{
    // check for errors
    check_abort(0);

    int do_output = 0;
    int do_finalize = 0;

    int iter = 0;
    if (use_restart == 1)
        iter = iter_restart;

    double t = t_start;
    if (use_restart == 1)
        t = t_restart;

    double dt;

    print_residual_header();
    if (iter == 0)
        fv_time_derivative(t);
    if ((do_output_data == 1) && (use_restart == 0))
        write_output(iter, t);

    while (1)
    {
        // check for errors
        check_abort(0);

        // calculate time step dt
        double dt_local = calc_time_step_function_pointer();
        mpi_all_reduce(&dt_local, &dt, MPIDouble, MPIMin);

        if (t + dt > t_end)
        {
            dt = t_end - t;
            do_output = 1;
            do_finalize = 1;
        }

        // the timestep to be called
        time_step_function_pointer(iter, t, dt);
        calc_global_residual(dt);

        // check for NAN and INF
        if (is_nan_n(residual, all_variables->n_sol_variables) ||
            is_inf_n(residual, all_variables->n_sol_variables))
            check_error(0);

        t = t + dt;
        iter = iter + 1;

        // steady-state simulation
        if ((is_transient == 0) && (min_n(residual, all_variables->n_sol_variables) < abort_residual))
        {
            t_end = t;
            do_output = 1;
            do_finalize = 1;
        }

        if ((i_output_data > 0) && (iter % i_output_data == 0))
        {
            do_output = 1;
        }

        // maximum iteration number reacher
        if (iter >= iter_restart + max_iter)
        {
            t_end = t;
            do_output = 1;
            do_finalize = 1;
        }

        print_residual(iter, t, dt, do_output);

        // output
        if ((do_output_data == 1) && (do_output == 1))
        {
            write_output(iter, t);
            do_output = 0;
        }

        // stop if required
        if (do_finalize)
            break;
    }
}

void print_residual_header()
{
    if (explicit_active)
    {
        printf_r("%9s %12s %12s %1s %1s:", "iter", "time", "dt", "V", "O");
    }
    else
    {
        printf_r("%9s %12s %12s %1s %1s %6s %6s:", "iter", "time", "dt", "V", "O", "inner", "lsoe");
    }

    for (int i = 0; i < all_variables->n_sol_variables; ++i)
        printf_r(" %12s", (&all_variables->sol_variables[i])->name);

    printf_r("\n");
}

void print_residual(int iter, double t, double dt, int do_output)
{
    char output_str = (do_output == 1) ? '*' : ' ';
    char viscous_str = (is_viscous_dt == 1) ? 'T' : 'F';

    if (explicit_active)
    {
        printf_r("%09d %12.5e %12.5e %c %c:", iter, t, dt, viscous_str, output_str);
    }
    else
    {
        printf_r("%09d %12.5e %12.5e %c %c %6d %6d:", iter, t, dt, viscous_str, output_str, n_iter_inner, n_iter_lsoe);
    }

    for (int i = 0; i < all_variables->n_sol_variables; ++i)
        printf_r(" %12.5e", residual[i]);

    printf_r("\n");
}