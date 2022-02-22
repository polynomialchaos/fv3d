/*******************************************************************************
 * @file analyze_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef ANALYZE_MODULE_H
#define ANALYZE_MODULE_H

#include "fv3d/fv3d_module.h"

/*******************************************************************************
 * @brief Define analyze
 ******************************************************************************/
void analyze_define();

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void analyze_finalize();

/*******************************************************************************
 * @brief Initialize analyze
 ******************************************************************************/
void analyze_initialize();

#endif /* ANALYZE_MODULE_H */