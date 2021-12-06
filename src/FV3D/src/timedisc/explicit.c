/*******************************************************************************
 * @file explicit.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "timedisc/timedisc_module.h"

string_t explicit_scheme_name = NULL;

/*******************************************************************************
 * @brief Define explicit
 ******************************************************************************/
void explicit_define()
{
    REGISTER_INITIALIZE_ROUTINE(explicit_initialize);
    REGISTER_FINALIZE_ROUTINE(explicit_finalize);

    string_t tmp_opt[] = {"RK3-3", "RK4-5", "Euler"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    SET_PARAMETER("TimeDisc/Explicit/scheme", StringParameter, &tmp, "The explicit timestep scheme", &tmp_opt, tmp_opt_n);
}

/*******************************************************************************
 * @brief Finalize explicit
 ******************************************************************************/
void explicit_finalize()
{
    DEALLOCATE(explicit_scheme_name);
    free_explicit();
}

/*******************************************************************************
 * @brief Initialize explicit
 ******************************************************************************/
void explicit_initialize()
{
    if (is_explicit() == BFLS)
        return;

    GET_PARAMETER("TimeDisc/Explicit/scheme", StringParameter, &explicit_scheme_name);

    if (is_equal(explicit_scheme_name, "Euler"))
    {
        init_explicit(EulerExplicit);
    }
    else if (is_equal(explicit_scheme_name, "RK3-3"))
    {
        init_explicit(RungeKutta33);
    }
    else if (is_equal(explicit_scheme_name, "RK4-5"))
    {
        init_explicit(RungeKutta45);
    }
    else
    {
        CHECK_EXPRESSION(0);
    }
}