/*******************************************************************************
 * @file implicit_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef IMPLICIT_MODULE_H
#define IMPLICIT_MODULE_H

#include "fv3d/fv3d_module.h"

extern int implicit_active;
extern int n_iter_inner;
extern int n_iter_lsoe;
extern int n_bdf_stages;
extern double **phi_old;

/*******************************************************************************
 * @brief Define implicit
 ******************************************************************************/
void implicit_define();

/*******************************************************************************
 * @brief Finalize implicit
 ******************************************************************************/
void implicit_finalize();

/*******************************************************************************
 * @brief Initialize implicit
 ******************************************************************************/
void implicit_initialize();

#endif /* IMPLICIT_MODULE_H */