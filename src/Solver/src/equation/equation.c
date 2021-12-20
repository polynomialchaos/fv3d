/*******************************************************************************
 * @file equation.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "equation_private.h"

Variables_t *solver_variables = NULL;

/*******************************************************************************
 * @brief Allocate the variables
 ******************************************************************************/
Variables_t *allocate_variables()
{
    Variables_t *tmp = ALLOCATE(sizeof(Variables_t));

    tmp->n_sol_variables = 0;
    tmp->n_dep_variables = 0;
    tmp->tot_variables = 0;

    tmp->sol_variables = NULL;
    tmp->dep_variables = NULL;
    tmp->tot_variables = NULL;

    return tmp;
}

/*******************************************************************************
 * @brief Add a solution variable and return array index
 * @param name
 * @return int
 ******************************************************************************/
int add_sol_variable(cstring_t name)
{
    CHECK_EXPRESSION(solver_variables->n_dep_variables == 0);

    solver_variables->n_sol_variables += 1;
    solver_variables->sol_variables =
        REALLOCATE(solver_variables->sol_variables,
                   sizeof(Variable_t) * solver_variables->n_sol_variables);

    Variable_t *tmp =
        &solver_variables->sol_variables[solver_variables->n_sol_variables - 1];
    init_variable(tmp, name);

    set_tot_variables(solver_variables);
    return solver_variables->n_tot_variables - 1;
}

/*******************************************************************************
 * @brief Add a dependent variable and return array index
 * @param name
 * @return int
 ******************************************************************************/
int add_dep_variable(cstring_t name)
{
    solver_variables->n_dep_variables += 1;
    solver_variables->dep_variables =
        REALLOCATE(solver_variables->dep_variables,
                   sizeof(Variable_t) * solver_variables->n_dep_variables);

    Variable_t *tmp =
        &solver_variables->dep_variables[solver_variables->n_dep_variables - 1];
    init_variable(tmp, name);

    set_tot_variables(solver_variables);
    return solver_variables->n_tot_variables - 1;
}

/*******************************************************************************
 * @brief Deallocate the variable
 * @param variable
 ******************************************************************************/
void deallocate_variable(Variable_t *variable)
{
    DEALLOCATE(variable->name);
}

/*******************************************************************************
 * @brief Deallocate the variables
 * @param variables
 ******************************************************************************/
void deallocate_variables(Variables_t *variables)
{
    if (variables == NULL)
        return;

    for (int i = 0; i < variables->n_sol_variables; ++i)
        deallocate_variable(&variables->sol_variables[i]);
    DEALLOCATE(variables->sol_variables);

    for (int i = 0; i < variables->n_dep_variables; ++i)
        deallocate_variable(&variables->dep_variables[i]);
    DEALLOCATE(variables->dep_variables);

    DEALLOCATE(variables->tot_variables);
}

/*******************************************************************************
 * @brief Free equation
 ******************************************************************************/
void free_equation()
{
    deallocate_variables(solver_variables);
    DEALLOCATE(solver_variables);
}

/*******************************************************************************
 * @brief Initialize equation
 ******************************************************************************/
void init_equation()
{
    solver_variables = allocate_variables();
}

/*******************************************************************************
 * @brief Initialize an allocated variable
 * @param variable
 * @param name
 ******************************************************************************/
void init_variable(Variable_t *variable, cstring_t name)
{
    variable->name = allocate_strcpy(name);
}

/*******************************************************************************
 * @brief Print the variables
 ******************************************************************************/
void print_variables()
{
    PRINTF("VARIABLES\n");

    PRINTF("n_sol_variables = %d\n", solver_variables->n_sol_variables);
    for (int i = 0; i < solver_variables->n_sol_variables; ++i)
        PRINTF("%d: %s\n", i, (&solver_variables->sol_variables[i])->name);

    PRINTF("n_dep_variables  = %d\n", solver_variables->n_dep_variables);
    for (int i = 0; i < solver_variables->n_dep_variables; ++i)
        PRINTF("%d: %s\n", i, (&solver_variables->dep_variables[i])->name);

    PRINTF("n_tot_variables = %d\n", solver_variables->n_tot_variables);
    for (int i = 0; i < solver_variables->n_tot_variables; ++i)
        PRINTF("%d: %s\n", i, solver_variables->tot_variables[i]->name);
}

/*******************************************************************************
 * @brief Set the total variables
 ******************************************************************************/
void set_tot_variables()
{
    solver_variables->n_tot_variables =
        solver_variables->n_sol_variables + solver_variables->n_dep_variables;
    solver_variables->tot_variables =
        REALLOCATE(solver_variables->tot_variables,
                   sizeof(Variable_t *) * solver_variables->n_tot_variables);

    for (int i = 0; i < solver_variables->n_sol_variables; ++i)
        solver_variables->tot_variables[i] =
            &solver_variables->sol_variables[i];

    for (int i = 0; i < solver_variables->n_dep_variables; ++i)
        solver_variables->tot_variables
            [solver_variables->n_sol_variables + i] =
            &solver_variables->dep_variables[i];
}
