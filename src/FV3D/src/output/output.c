/*******************************************************************************
 * @file output.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "output_module.h"
// #include "fv3d_private.h"
// #include "mesh/mesh_module.h"
//
// #include "fv/fv_module.h"
// #include "timedisc/implicit/implicit_module.h"



/*******************************************************************************
 * @brief Define output
 ******************************************************************************/
void output_define()
{
    REGISTER_INITIALIZE_ROUTINE(output_initialize);
    REGISTER_FINALIZE_ROUTINE(output_finalize);

    int i_output_data = -1;
    SET_PARAMETER("Output/i_output_data", DigitParameter, &i_output_data,
                  "The output file frequency "
                  "(-1 ... first/solutions/last, 0 ... disable)",
                  NULL, 0);
}

/*******************************************************************************
 * @brief Finalize output
 ******************************************************************************/
void output_finalize()
{
    free_output();
}

/*******************************************************************************
 * @brief Initialize output
 ******************************************************************************/
void output_initialize()
{
    int i_output_data;
    GET_PARAMETER("Output/i_output_data", DigitParameter, &i_output_data);
    init_output(i_output_data);
}
