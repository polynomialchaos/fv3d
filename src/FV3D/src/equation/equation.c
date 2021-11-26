/*******************************************************************************
 * @file equation.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "equation_module.h"
#include "navier_stokes/navier_stokes_module.h"

string_t equation_name = NULL;

Variables_t *all_variables = NULL;

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
 * @brief Add a solution variable
 * @param variables
 * @param name
 * @return int
 ******************************************************************************/
int add_sol_variable(Variables_t *variables, string_t name)
{
    variables->n_sol_variables += 1;
    variables->sol_variables = REALLOCATE(variables->sol_variables, sizeof(Variable_t) * variables->n_sol_variables);

    Variable_t *tmp = &variables->sol_variables[variables->n_sol_variables - 1];
    tmp->name = allocate_strcpy(name);

    set_tot_variables(variables);
    return variables->n_tot_variables - 1;
}

/*******************************************************************************
 * @brief Add a dependent variable
 * @param variables
 * @param name
 * @return int
 ******************************************************************************/
int add_dep_variable(Variables_t *variables, string_t name)
{
    variables->n_dep_variables += 1;
    variables->dep_variables = REALLOCATE(variables->dep_variables, sizeof(Variable_t) * variables->n_dep_variables);

    Variable_t *tmp = &variables->dep_variables[variables->n_dep_variables - 1];
    tmp->name = allocate_strcpy(name);

    set_tot_variables(variables);
    return variables->n_tot_variables - 1;
}

/*******************************************************************************
 * @brief Deallocate the variables
 * @param variables
 ******************************************************************************/
void deallocate_variables(Variables_t **variables)
{
    if (*variables == NULL)
        return;

    for (int i = 0; i < (*variables)->n_sol_variables; ++i)
        DEALLOCATE((&(*variables)->sol_variables[i])->name);

    DEALLOCATE((*variables)->sol_variables);

    for (int i = 0; i < (*variables)->n_dep_variables; ++i)
        DEALLOCATE((&(*variables)->dep_variables[i])->name);

    DEALLOCATE((*variables)->dep_variables);

    DEALLOCATE((*variables)->tot_variables);

    DEALLOCATE((*variables));
}

/*******************************************************************************
 * @brief Define equation
 ******************************************************************************/
void equation_define()
{
    REGISTER_INITIALIZE_ROUTINE(equation_initialize);
    REGISTER_FINALIZE_ROUTINE(equation_finalize);

    string_t tmp_opt[] = {"Navier-Stokes"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    SET_PARAMETER("Equation/equation", StringParameter, &tmp,
                  "The eqaution to solve", &tmp_opt, tmp_opt_n);

    navier_stokes_define();
}

/*******************************************************************************
 * @brief Finalize equation
 ******************************************************************************/
void equation_finalize()
{
    deallocate_variables(&all_variables);

    DEALLOCATE(equation_name);
}

/*******************************************************************************
 * @brief Initialize equation
 ******************************************************************************/
void equation_initialize()
{
    GET_PARAMETER("Equation/equation", StringParameter, &equation_name);

    if (is_equal(equation_name, "Navier-Stokes"))
    {
        navier_stokes_active = 1;
    }
    else
    {
        CHECK_EXPRESSION(0);
    }

    all_variables = allocate_variables();
}

/*******************************************************************************
 * @brief Print the variables
 * @param variables
 ******************************************************************************/
void print_variables(Variables_t *variables)
{
    PRINTF("VARIABLES\n");

    PRINTF("n_sol_variables = %d\n", variables->n_sol_variables);
    for (int i = 0; i < variables->n_sol_variables; ++i)
        PRINTF("%d: %s\n", i, (&variables->sol_variables[i])->name);

    PRINTF("n_dep_variables  = %d\n", variables->n_dep_variables);
    for (int i = 0; i < variables->n_dep_variables; ++i)
        PRINTF("%d: %s\n", i, (&variables->dep_variables[i])->name);

    PRINTF("n_tot_variables = %d\n", variables->n_tot_variables);
    for (int i = 0; i < variables->n_tot_variables; ++i)
        PRINTF("%d: %s\n", i, variables->tot_variables[i]->name);
}

/*******************************************************************************
 * @brief Set the total variables
 * @param variables
 ******************************************************************************/
void set_tot_variables(Variables_t *variables)
{
    variables->n_tot_variables = variables->n_sol_variables + variables->n_dep_variables;
    variables->tot_variables = REALLOCATE(variables->tot_variables, sizeof(Variable_t *) * variables->n_tot_variables);

    for (int i = 0; i < variables->n_sol_variables; ++i)
        variables->tot_variables[i] = &variables->sol_variables[i];

    for (int i = 0; i < variables->n_dep_variables; ++i)
        variables->tot_variables[variables->n_sol_variables + i] = &variables->dep_variables[i];
}
