/*******************************************************************************
 * @file analyze.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "analyze_module.h"

/*******************************************************************************
 * @brief Define analyze
 ******************************************************************************/
void analyze_define()
{
    REGISTER_INITIALIZE_ROUTINE(analyze_initialize);
    REGISTER_FINALIZE_ROUTINE(analyze_finalize);
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
