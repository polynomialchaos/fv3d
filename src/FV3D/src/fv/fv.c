/*******************************************************************************
 * @file fv.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "fv_module.h"

/*******************************************************************************
 * @brief Define fv
 ******************************************************************************/
void fv_define()
{
    REGISTER_INITIALIZE_ROUTINE(fv_initialize);
    REGISTER_FINALIZE_ROUTINE(fv_finalize);

    reconstruction_define();
    limiter_define();
}

/*******************************************************************************
 * @brief Finalize fv
 ******************************************************************************/
void fv_finalize()
{
    free_finite_volume();
}

/*******************************************************************************
 * @brief Initialize fv
 ******************************************************************************/
void fv_initialize()
{
    init_finite_volume();
}
