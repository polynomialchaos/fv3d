//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "equation_module.h"
#include "navier_stokes/navier_stokes_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
string_t equation_name = NULL;

Variables_t *all_variables = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_initialize();
void equation_finalize();

Variables_t *allocate_variables();
void set_tot_variables(Variables_t *variables);
void print_variables(Variables_t *variables);
void deallocate_variables(Variables_t **tmp);

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_define()
{
    register_initialize_routine(equation_initialize);
    register_finalize_routine(equation_finalize);

    string_t tmp_opt[] = {"Navier-Stokes"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    set_parameter("Equation/equation", ParameterString, &tmp,
                  "The eqaution to solve", &tmp_opt, tmp_opt_n);

    navier_stokes_define();
}

void equation_initialize()
{
    get_parameter("Equation/equation", ParameterString, &equation_name);

    if (is_equal(equation_name, "Navier-Stokes"))
    {
        navier_stokes_active = 1;
    }
    else
    {
        check_error(0);
    }

    all_variables = allocate_variables();
}

void equation_finalize()
{
    deallocate_variables(&all_variables);

    deallocate(equation_name);
}

int add_sol_variable(Variables_t *variables, string_t name)
{
    variables->n_sol_variables += 1;
    variables->sol_variables = reallocate(variables->sol_variables, sizeof(Variable_t) * variables->n_sol_variables);

    Variable_t *tmp = &variables->sol_variables[variables->n_sol_variables - 1];
    tmp->name = allocate_strcpy(name);

    set_tot_variables(variables);
    return variables->n_tot_variables - 1;
}

int add_dep_variable(Variables_t *variables, string_t name)
{
    variables->n_dep_variables += 1;
    variables->dep_variables = reallocate(variables->dep_variables, sizeof(Variable_t) * variables->n_dep_variables);

    Variable_t *tmp = &variables->dep_variables[variables->n_dep_variables - 1];
    tmp->name = allocate_strcpy(name);

    set_tot_variables(variables);
    return variables->n_tot_variables - 1;
}

Variables_t *allocate_variables()
{
    Variables_t *tmp = allocate(sizeof(Variables_t));

    tmp->n_sol_variables = 0;
    tmp->n_dep_variables = 0;
    tmp->tot_variables = 0;

    tmp->sol_variables = NULL;
    tmp->dep_variables = NULL;
    tmp->tot_variables = NULL;

    return tmp;
}

void set_tot_variables(Variables_t *variables)
{
    variables->n_tot_variables = variables->n_sol_variables + variables->n_dep_variables;
    variables->tot_variables = reallocate(variables->tot_variables, sizeof(Variable_t *) * variables->n_tot_variables);

    for (int i = 0; i < variables->n_sol_variables; i++)
        variables->tot_variables[i] = &variables->sol_variables[i];

    for (int i = 0; i < variables->n_dep_variables; i++)
        variables->tot_variables[variables->n_sol_variables + i] = &variables->dep_variables[i];
}

void print_variables(Variables_t *variables)
{
    printf_r("VARIABLES\n");

    printf_r("n_sol_variables = %d\n", variables->n_sol_variables);
    for (int i = 0; i < variables->n_sol_variables; i++)
        printf_r("%d: %s\n", i, (&variables->sol_variables[i])->name);

    printf_r("n_dep_variables  = %d\n", variables->n_dep_variables);
    for (int i = 0; i < variables->n_dep_variables; i++)
        printf_r("%d: %s\n", i, (&variables->dep_variables[i])->name);

    printf_r("n_tot_variables = %d\n", variables->n_tot_variables);
    for (int i = 0; i < variables->n_tot_variables; i++)
        printf_r("%d: %s\n", i, variables->tot_variables[i]->name);
}

void deallocate_variables(Variables_t **variables)
{
    if (*variables == NULL)
        return;

    for (int i = 0; i < (*variables)->n_sol_variables; i++)
        deallocate((&(*variables)->sol_variables[i])->name);

    deallocate((*variables)->sol_variables);

    for (int i = 0; i < (*variables)->n_dep_variables; i++)
        deallocate((&(*variables)->dep_variables[i])->name);

    deallocate((*variables)->dep_variables);

    deallocate((*variables)->tot_variables);

    deallocate((*variables));
}