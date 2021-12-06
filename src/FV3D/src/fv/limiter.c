/*******************************************************************************
 * @file limiter.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "fv_module.h"

string_t limiter_name = NULL;

/*******************************************************************************
 * @brief Define limiter
 ******************************************************************************/
void limiter_define()
{
    REGISTER_INITIALIZE_ROUTINE(limiter_initialize);
    REGISTER_FINALIZE_ROUTINE(limiter_finalize);

    string_t tmp_opt[] = {"Barth-Jespersenn", "None"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    SET_PARAMETER("FV/Limiter/limiter", StringParameter, &tmp,
                  "The limiter method", &tmp_opt, tmp_opt_n);
}

/*******************************************************************************
 * @brief Finalize limiter
 ******************************************************************************/
void limiter_finalize()
{
    DEALLOCATE(limiter_name);
    free_limiter();
}

/*******************************************************************************
 * @brief Initialize limiter
 ******************************************************************************/
void limiter_initialize()
{
    GET_PARAMETER("FV/Limiter/limiter", StringParameter, &limiter_name);

    if (is_equal(limiter_name, "Barth-Jespersenn"))
    {
        init_limiter(BarthJespersenn);
    }
    else
    {
        init_limiter(NoLimter);
    }
}
