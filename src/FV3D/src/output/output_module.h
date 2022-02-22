/*******************************************************************************
 * @file output_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef OUTPUT_MODULE_H
#define OUTPUT_MODULE_H

#include "fv3d/fv3d_module.h"

/*******************************************************************************
 * @brief Define output
 ******************************************************************************/
void output_define();

/*******************************************************************************
 * @brief Finalize output
 ******************************************************************************/
void output_finalize();

/*******************************************************************************
 * @brief Initialize output
 ******************************************************************************/
void output_initialize();

#endif /* OUTPUT_MODULE_H */