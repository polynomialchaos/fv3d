/*******************************************************************************
 * @file equation_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef EQUATION_MODULE_H
#define EQUATION_MODULE_H

#include "fv3d/fv3d_module.h"

typedef struct Variable
{
    string_t name;
} Variable_t;

typedef struct Variables
{
    int n_sol_variables;
    int n_dep_variables;
    int n_tot_variables;

    Variable_t *sol_variables;
    Variable_t *dep_variables;
    Variable_t **tot_variables;
} Variables_t;

extern Variables_t *all_variables;

/*******************************************************************************
 * @brief Allocate the variables
 ******************************************************************************/
Variables_t *allocate_variables();

/*******************************************************************************
 * @brief Add a solution variable
 * @param variables
 * @param name
 * @return int
 ******************************************************************************/
int add_sol_variable(Variables_t *variables, string_t name);

/*******************************************************************************
 * @brief Add a dependent variable
 * @param variables
 * @param name
 * @return int
 ******************************************************************************/
int add_dep_variable(Variables_t *variables, string_t name);

/*******************************************************************************
 * @brief Deallocate the variables
 * @param variables
 ******************************************************************************/
void deallocate_variables(Variables_t *variables);

/*******************************************************************************
 * @brief Define equation
 ******************************************************************************/
void equation_define();

/*******************************************************************************
 * @brief Finalize equation
 ******************************************************************************/
void equation_finalize();

/*******************************************************************************
 * @brief Initialize equation
 ******************************************************************************/
void equation_initialize();

/*******************************************************************************
 * @brief Set the total variables
 * @param variables
 ******************************************************************************/
void set_tot_variables(Variables_t *variables);

#endif /* EQUATION_MODULE_H */