/*******************************************************************************
 * @file restart_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef RESTART_MODULE_H
#define RESTART_MODULE_H

#include "fv3d/fv3d_module.h"

/*******************************************************************************
 * @brief Define restart
 ******************************************************************************/
void restart_define();

/*******************************************************************************
 * @brief Finalize restart
 ******************************************************************************/
void restart_finalize();

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void restart_initialize();

#endif /* RESTART_MODULE_H */