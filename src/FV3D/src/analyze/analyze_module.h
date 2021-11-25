/*******************************************************************************
 * @file analyze_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef ANALYZE_MODULE_H
#define ANALYZE_MODULE_H

#include "fv3d/fv3d_module.h"

extern double *residual;

void analyze_define();

void calc_global_residual(double dt);

#endif /* ANALYZE_MODULE_H */