/*******************************************************************************
 * @file timedisc_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef TIMEDISC_MODULE_H
#define TIMEDISC_MODULE_H

#include "fv3d/fv3d_module.h"

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

/*******************************************************************************
 * @brief Define timedisc
 ******************************************************************************/
void timedisc_define();

/*******************************************************************************
 * @brief Finalize timedisc
 ******************************************************************************/
void timedisc_finalize();

/*******************************************************************************
 * @brief Initialize timedisc
 ******************************************************************************/
void timedisc_initialize();

#endif /* TIMEDISC_MODULE_H */