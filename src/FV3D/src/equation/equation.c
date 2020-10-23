//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "equation_private.h"
#include "navier_stokes/navier_stokes_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
string_t equation_name      = NULL;

Variables_t *all_variables  = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_initialize();
void equation_finalize();

Variables_t *allocate_variables();
void deallocate_variables( Variables_t **tmp );

int add_sol_variable( Variables_t *variables, string_t name );
int add_dep_variable( Variables_t *variables, string_t name );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_define()
{
    register_initialize_routine( equation_initialize );
    register_finalize_routine( equation_finalize );

    string_t tmp_opt[] = {"Navier-Stokes"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    set_parameter( "Equation/equation", ParameterString, &tmp,
        "The eqaution to solve", &tmp_opt, tmp_opt_n );

    navier_stokes_define();
}

void equation_initialize()
{
    get_parameter( "Equation/equation", ParameterString, &equation_name );

    if (is_equal( equation_name, "Navier-Stokes" ))
    {
        navier_stokes_active = 1;
    }
    else
    {
        check_error( 0 );
    }

    all_variables = allocate_variables();
}

void equation_finalize()
{
    deallocate( equation_name );
    deallocate_variables( &all_variables );

    // nullify( exact_func_routine )
    // nullify( update_routine )
    // nullify( update_gradients_routine )
    // nullify( calc_time_step_routine )
    // nullify( calc_flux_routine )
}

int add_sol_variable( Variables_t *variables, string_t name )
{
    variables->n_sol_variables  += 1;
    variables->sol_variables = reallocate( variables->sol_variables, sizeof( Variable_t ) * variables->n_sol_variables );

    int i_var = variables->n_sol_variables - 1;
    Variable_t *tmp = &variables->sol_variables[i_var];
    tmp->name = allocate_strcpy( name );

    return i_var;
}

int add_dep_variable( Variables_t *variables, string_t name )
{
    variables->n_dep_variables  += 1;
    variables->dep_variables     = reallocate( variables->dep_variables, sizeof( Variable_t ) * variables->n_dep_variables );

    int i_var = variables->n_dep_variables - 1;
    Variable_t *tmp = &variables->dep_variables[i_var];
    tmp->name = allocate_strcpy( name );

    return i_var;
}

Variables_t *allocate_variables()
{
    Variables_t *tmp = allocate( sizeof( Variables_t ) );

    tmp->n_sol_variables    = 0;
    tmp->n_dep_variables    = 0;

    tmp->sol_variables  = NULL;
    tmp->dep_variables  = NULL;

    return tmp;
}

void deallocate_variables( Variables_t **variables )
{
    if (*variables == NULL) return;

    for (int i = 0; i < (*variables)->n_sol_variables; i++ )
        deallocate( (&(*variables)->sol_variables[i])->name );

    deallocate( (*variables)->sol_variables );

    for (int i = 0; i < (*variables)->n_dep_variables; i++ )
        deallocate( (&(*variables)->dep_variables[i])->name );

    deallocate( (*variables)->dep_variables );

    deallocate( (*variables) );
}