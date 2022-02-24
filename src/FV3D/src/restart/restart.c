/*******************************************************************************
 * @file restart.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "restart_module.h"

bool_t use_restart = BC_FALSE;
string_t restart_file = NULL;

/*******************************************************************************
 * @brief Define restart
 ******************************************************************************/
void restart_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(restart_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(restart_finalize);

    string_t tmp = "untitled.h5";

    BM_SET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart,
                     "The flag to start from restart", NULL, 0);
    BM_SET_PARAMETER("Restart/restart_file", StringParameter, &tmp,
                     "The restart file", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize restart
 ******************************************************************************/
void restart_finalize()
{
    free_restart();
    BM_DEALLOCATE(restart_file);
}

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void restart_initialize()
{
    BM_GET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart);
    BM_GET_PARAMETER("Restart/restart_file", StringParameter, &restart_file);
    init_restart(use_restart, restart_file);
}
