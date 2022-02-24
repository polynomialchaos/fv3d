/*******************************************************************************
 * @file timedisc.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "timedisc_module.h"

string_t timestep_name = NULL;
int max_iter = 1000;
bool_t is_transient = BC_FALSE;
double abort_residual = 1e-10;
double t_start = 0.0;
double t_end = 0.2;

/*******************************************************************************
 * @brief Define timedisc
 ******************************************************************************/
void timedisc_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(timedisc_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(timedisc_finalize);

    string_t tmp_opt[] = {"Explicit", "Implicit"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    BM_SET_PARAMETER("TimeDisc/timestep", StringParameter, &tmp,
                     "The timestep mehtod", &tmp_opt, tmp_opt_n);
    BM_SET_PARAMETER("TimeDisc/max_iter", DigitParameter, &max_iter,
                     "The maximum number of iterations", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/transient", LogicalParameter, &is_transient,
                     "The flag wheter to be transient or steady-state",
                     NULL, 0);
    BM_SET_PARAMETER("TimeDisc/abort_residual", NumberParameter,
                     &abort_residual,
                     "The abort residual", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/t_start", NumberParameter, &t_start,
                     "The start time", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/t_end", NumberParameter, &t_end,
                     "The end time", NULL, 0);

    explicit_define();
    implicit_define();
}

/*******************************************************************************
 * @brief Finalize timedisc
 ******************************************************************************/
void timedisc_finalize()
{
    BM_DEALLOCATE(timestep_name);
    free_timedisc();
}

/*******************************************************************************
 * @brief Initialize timedisc
 ******************************************************************************/
void timedisc_initialize()
{
    BM_GET_PARAMETER("TimeDisc/timestep", StringParameter, &timestep_name);
    BM_GET_PARAMETER("TimeDisc/max_iter", DigitParameter, &max_iter);
    BM_GET_PARAMETER("TimeDisc/transient", LogicalParameter, &is_transient);
    BM_GET_PARAMETER("TimeDisc/abort_residual", NumberParameter,
                     &abort_residual);
    BM_GET_PARAMETER("TimeDisc/t_start", NumberParameter, &t_start);
    BM_GET_PARAMETER("TimeDisc/t_end", NumberParameter, &t_end);

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
        BM_CHECK_EXPRESSION(0);
    }
}
