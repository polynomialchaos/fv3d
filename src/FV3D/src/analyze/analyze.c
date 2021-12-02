/*******************************************************************************
 * @file analyze.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "analyze_module.h"
#include "mesh/mesh_module.h"
#include "fv/fv_module.h"

double *residual = NULL;

/*******************************************************************************
 * @brief Define analyze
 ******************************************************************************/
void analyze_define()
{
    REGISTER_INITIALIZE_ROUTINE(analyze_initialize);
    REGISTER_FINALIZE_ROUTINE(analyze_finalize);
}

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void analyze_finalize()
{
    DEALLOCATE(residual);
}

/*******************************************************************************
 * @brief Initialize analyze
 ******************************************************************************/
void analyze_initialize()
{
    int n_sol_variables = solver_variables->n_sol_variables;

    residual = ALLOCATE(sizeof(double) * n_sol_variables);
    set_value_n(0.0, n_sol_variables, residual);
}

/*******************************************************************************
 * @brief Calculate the global residual
 * @param dt
 ******************************************************************************/
void calc_global_residual(double dt)
{
    Cells_t *cells = solver_mesh->cells;
    int n_sol_variables = solver_variables->n_sol_variables;
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
        tmp[i] *= dt / solver_mesh->global_volume;

    mpi_all_reduce_n(MPIDouble, MPISum, tmp, n_sol_variables, residual);
}