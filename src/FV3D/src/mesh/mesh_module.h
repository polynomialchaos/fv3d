/*******************************************************************************
 * @file mesh_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef MESH_MODULE_H
#define MESH_MODULE_H

#include "fv3d/fv3d_module.h"

/*******************************************************************************
 * @brief Define mesh
 ******************************************************************************/
void mesh_define();

/*******************************************************************************
 * @brief Finalize mesh
 ******************************************************************************/
void mesh_finalize();

/*******************************************************************************
 * @brief Initialize mesh
 ******************************************************************************/
void mesh_initialize();

#endif /* MESH_MODULE_H */