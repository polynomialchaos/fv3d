/*******************************************************************************
 * @file analyze.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "analyze_module.h"

/*******************************************************************************
 * @brief Define analyze
 ******************************************************************************/
void analyze_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(analyze_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(analyze_finalize);
}

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void analyze_finalize()
{
    free_analyze();
}

/*******************************************************************************
 * @brief Initialize analyze
 ******************************************************************************/
void analyze_initialize()
{
    init_analyze();
}
