/*******************************************************************************
 * @file reconstruction.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "fv_module.h"

string_t reconstruction_name = NULL;

/*******************************************************************************
 * @brief Define reconstruction
 ******************************************************************************/
void reconstruction_define()
{
    REGISTER_INITIALIZE_ROUTINE(reconstruction_initialize);
    REGISTER_FINALIZE_ROUTINE(reconstruction_finalize);

    string_t tmp_opt[] = {"Linear", "First-Order"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    SET_PARAMETER("FV/Reconstruction/reconstruction", StringParameter, &tmp,
                  "The reconstruction method", &tmp_opt, tmp_opt_n);
}

/*******************************************************************************
 * @brief Finalize reconstruction
 ******************************************************************************/
void reconstruction_finalize()
{
    DEALLOCATE(reconstruction_name);
    free_reconstruction();
}

/*******************************************************************************
 * @brief Initialize reconstruction
 ******************************************************************************/
void reconstruction_initialize()
{
    GET_PARAMETER("FV/Reconstruction/reconstruction", StringParameter,
                  &reconstruction_name);

    if (is_equal(reconstruction_name, "First-Order"))
    {
        init_reconstruction(FirstOrder);
    }
    else if (is_equal(reconstruction_name, "Linear"))
    {
        init_reconstruction(Linear);
    }
    else
    {
        CHECK_EXPRESSION(0);
    }
}