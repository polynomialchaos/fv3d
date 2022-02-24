/*******************************************************************************
 * @file limiter.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "fv_module.h"

string_t limiter_name = NULL;

/*******************************************************************************
 * @brief Define limiter
 ******************************************************************************/
void limiter_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(limiter_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(limiter_finalize);

    string_t tmp_opt[] = {"Barth-Jespersenn", "None"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    BM_SET_PARAMETER("FV/Limiter/limiter", StringParameter, &tmp,
                     "The limiter method", &tmp_opt, tmp_opt_n);
}

/*******************************************************************************
 * @brief Finalize limiter
 ******************************************************************************/
void limiter_finalize()
{
    BM_DEALLOCATE(limiter_name);
    free_limiter();
}

/*******************************************************************************
 * @brief Initialize limiter
 ******************************************************************************/
void limiter_initialize()
{
    BM_GET_PARAMETER("FV/Limiter/limiter", StringParameter, &limiter_name);

    if (is_equal(limiter_name, "Barth-Jespersenn"))
    {
        init_limiter(BarthJespersenn);
    }
    else
    {
        init_limiter(NoLimter);
    }
}
