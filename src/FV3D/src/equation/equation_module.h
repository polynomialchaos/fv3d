//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef EQUATION_MODULE_H
#define EQUATION_MODULE_H

#include "fv3d_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
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

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_define();

int add_sol_variable( Variables_t *variables, string_t name );
int add_dep_variable( Variables_t *variables, string_t name );

#endif /* EQUATION_MODULE_H */