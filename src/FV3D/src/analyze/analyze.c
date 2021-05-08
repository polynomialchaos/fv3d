//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "analyze_module.h"
#include "equation/equation_module.h"
#include "mesh/mesh_module.h"
#include "fv/fv_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
double *residual = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void analyze_initialize();
void analyze_finalize();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void analyze_define()
{
    register_initialize_routine(analyze_initialize);
    register_finalize_routine(analyze_finalize);
}

void analyze_initialize()
{
    residual = allocate(sizeof(double) * all_variables->n_sol_variables);
    set_value_n(0.0, residual, all_variables->n_sol_variables);
}

void analyze_finalize()
{
    deallocate(residual);
}

void calc_global_residual(double dt)
{
    Cells_t *cells = global_mesh->cells;
    int n_sol_variables = all_variables->n_sol_variables;
    int n_domain_cells = cells->n_domain_cells;

    double tmp[n_sol_variables];
    set_value_n(0.0, tmp, n_sol_variables);

    for (int i = 0; i < n_domain_cells; i++)
    {
        double volume = cells->volume[i];

        for (int j = 0; j < n_sol_variables; j++)
        {
            tmp[j] += u_abs(phi_dt[i * n_sol_variables + j]) * volume;
        }
    }

    for (int i = 0; i < n_sol_variables; i++)
        tmp[i] *= dt / global_mesh->global_volume;

    mpi_all_reduce_n(tmp, residual, n_sol_variables, MPIDouble, MPISum);
}