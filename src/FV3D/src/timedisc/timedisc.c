
/*******************************************************************************
 * @file timedisc.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "timedisc_module.h"

string_t timestep_name = NULL;
int max_iter = 1000;
bool_t is_transient = BFLS;
double abort_residual = 1e-10;
double t_start = 0.0;
double t_end = 0.2;

/*******************************************************************************
 * @brief Define timedisc
 ******************************************************************************/
void timedisc_define()
{
    REGISTER_INITIALIZE_ROUTINE(timedisc_initialize);
    REGISTER_FINALIZE_ROUTINE(timedisc_finalize);

    string_t tmp_opt[] = {"Explicit", "Implicit"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    SET_PARAMETER("TimeDisc/timestep", StringParameter, &tmp,
                  "The timestep mehtod", &tmp_opt, tmp_opt_n);
    SET_PARAMETER("TimeDisc/max_iter", DigitParameter, &max_iter,
                  "The maximum number of iterations", NULL, 0);
    SET_PARAMETER("TimeDisc/transient", LogicalParameter, &is_transient,
                  "The flag wheter to be transient or steady-state", NULL, 0);
    SET_PARAMETER("TimeDisc/abort_residual", NumberParameter, &abort_residual,
                  "The abort residual", NULL, 0);
    SET_PARAMETER("TimeDisc/t_start", NumberParameter, &t_start,
                  "The start time", NULL, 0);
    SET_PARAMETER("TimeDisc/t_end", NumberParameter, &t_end,
                  "The end time", NULL, 0);

    explicit_define();
    implicit_define();
}

/*******************************************************************************
 * @brief Finalize timedisc
 ******************************************************************************/
void timedisc_finalize()
{
    DEALLOCATE(timestep_name);
    free_timedisc();
}

/*******************************************************************************
 * @brief Initialize timedisc
 ******************************************************************************/
void timedisc_initialize()
{
    GET_PARAMETER("TimeDisc/timestep", StringParameter, &timestep_name);
    GET_PARAMETER("TimeDisc/max_iter", DigitParameter, &max_iter);
    GET_PARAMETER("TimeDisc/transient", LogicalParameter, &is_transient);
    GET_PARAMETER("TimeDisc/abort_residual", NumberParameter, &abort_residual);
    GET_PARAMETER("TimeDisc/t_start", NumberParameter, &t_start);
    GET_PARAMETER("TimeDisc/t_end", NumberParameter, &t_end);

    if (is_equal(timestep_name, "Explicit"))
    {
        init_timedisc(Explicit, max_iter, is_transient,
                      abort_residual, t_start, t_end);
    }
    else if (is_equal(timestep_name, "Implicit"))
    {
        init_timedisc(Implicit, max_iter, is_transient,
                      abort_residual, t_start, t_end);
    }
    else
    {
        CHECK_EXPRESSION(0);
    }
}
