/*******************************************************************************
 * @file analyze.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "analyze_module.h"
#include "equation/equation_module.h"
#include "mesh/mesh_module.h"
#include "fv/fv_module.h"

double *residual = NULL;

void analyze_initialize();
void analyze_finalize();

void analyze_define()
{
    REGISTER_INITIALIZE_ROUTINE(analyze_initialize);
    REGISTER_FINALIZE_ROUTINE(analyze_finalize);
}

void analyze_initialize()
{
    int n_sol_variables = all_variables->n_sol_variables;

    residual = ALLOCATE(sizeof(double) * n_sol_variables);
    set_value_n(0.0, n_sol_variables, residual);
}

void analyze_finalize()
{
    DEALLOCATE(residual);
}

void calc_global_residual(double dt)
{
    Cells_t *cells = global_mesh->cells;
    int n_sol_variables = all_variables->n_sol_variables;
    int n_domain_cells = cells->n_domain_cells;

    double tmp[n_sol_variables];
    set_value_n(0.0, n_sol_variables, tmp);

    for (int i = 0; i < n_domain_cells; ++i)
    {
        double volume = cells->volume[i];

        for (int j = 0; j < n_sol_variables; ++j)
        {
            tmp[j] += ABS(phi_dt[i * n_sol_variables + j]) * volume;
        }
    }

    for (int i = 0; i < n_sol_variables; ++i)
        tmp[i] *= dt / global_mesh->global_volume;

    mpi_all_reduce_n(MPIDouble, MPISum, tmp, n_sol_variables, residual);
}