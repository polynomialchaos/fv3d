/*******************************************************************************
 * @file fv3d_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef FV3D_PRIVATE_H
#define FV3D_PRIVATE_H

#include "fv3d/fv3d_module.h"

extern string_t title;

/*******************************************************************************
 * @brief Define fv3d
 ******************************************************************************/
void fv3d_define();

/*******************************************************************************
 * @brief Finalize fv3d
 ******************************************************************************/
void fv3d_finalize();

/*******************************************************************************
 * @brief Initialize fv3d
 ******************************************************************************/
void fv3d_initialize();

#endif /* FV3D_PRIVATE_H */