/*******************************************************************************
 * @file explicit_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef EXPLICIT_MODULE_H
#define EXPLICIT_MODULE_H

#include "fv3d/fv3d_module.h"

extern int explicit_active;

/*******************************************************************************
 * @brief Define explicit
 ******************************************************************************/
void explicit_define();

/*******************************************************************************
 * @brief Finalize explicit
 ******************************************************************************/
void explicit_finalize();

/*******************************************************************************
 * @brief Initialize explicit
 ******************************************************************************/
void explicit_initialize();

/*******************************************************************************
 * @brief Explicit time discretizazion routine (LSERKW2)
 * @param iter
 * @param t
 * @param dt
 ******************************************************************************/
void time_step_lserkw2(int iter, double t, double dt);

#endif /* EXPLICIT_MODULE_H */