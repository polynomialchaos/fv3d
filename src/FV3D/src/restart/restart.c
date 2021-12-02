/*******************************************************************************
 * @file restart.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "restart_module.h"

/*******************************************************************************
 * @brief Define restart
 ******************************************************************************/
void restart_define()
{
    REGISTER_INITIALIZE_ROUTINE(restart_initialize);
    REGISTER_FINALIZE_ROUTINE(restart_finalize);

    int use_restart = 0;
    string_t restart_file = "untitled.h5";
    SET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart,
                  "The flag to start from restart", NULL, 0);
    SET_PARAMETER("Restart/restart_file", StringParameter, &restart_file,
                  "The restart file", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize restart
 ******************************************************************************/
void restart_finalize()
{
    free_restart();
}

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void restart_initialize()
{
    bool_t use_restart;
    string_t restart_file = NULL;
    GET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart);
    GET_PARAMETER("Restart/restart_file", StringParameter, &restart_file);
    init_restart(use_restart, restart_file);
}
