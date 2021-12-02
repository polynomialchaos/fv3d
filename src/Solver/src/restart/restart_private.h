/*******************************************************************************
 * @file restart_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef RESTART_PRIVATE_H
#define RESTART_PRIVATE_H

#include "solver/solver_module.h"

/*******************************************************************************
 * @brief Free restart
 ******************************************************************************/
void free_restart();

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void init_restart(bool_t use_restart, cstring_t restart_file);

/*******************************************************************************
 * @brief Read restart data
 * @param restart_file
 ******************************************************************************/
void read_restart_data(cstring_t restart_file);

#endif /* RESTART_PRIVATE_H */