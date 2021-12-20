/*******************************************************************************
 * @file equation_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef EQUATION_PRIVATE_H
#define EQUATION_PRIVATE_H

#include "solver/solver_module.h"

/*******************************************************************************
 * @brief Allocate the variables
 ******************************************************************************/
Variables_t *allocate_variables();

/*******************************************************************************
 * @brief Deallocate the variables
 * @param variables
 ******************************************************************************/
void deallocate_variables(Variables_t *variables);

/*******************************************************************************
 * @brief Deallocate the variable
 * @param variable
 ******************************************************************************/
void deallocate_variable(Variable_t *variable);

/*******************************************************************************
 * @brief Initialize an allocated variable
 * @param variable
 * @param name
 ******************************************************************************/
void init_variable(Variable_t *variable, cstring_t name);

/*******************************************************************************
 * @brief Set the total variables
 ******************************************************************************/
void set_tot_variables();

#endif /* EQUATION_PRIVATE_H */