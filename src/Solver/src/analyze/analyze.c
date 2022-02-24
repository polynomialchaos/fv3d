/*******************************************************************************
 * @file analyze.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "analyze_private.h"

double *solver_residual = NULL;

/*******************************************************************************
 * @brief Free analyze
 ******************************************************************************/
void free_analyze()
{
    BM_DEALLOCATE(solver_residual);
}

/*******************************************************************************
 * @brief Initialize analyze
 ******************************************************************************/
void init_analyze()
{
    int n_sol_variables = solver_variables->n_sol_variables;

    solver_residual = BM_ALLOCATE(sizeof(double) * n_sol_variables);
    set_value_n(0.0, n_sol_variables, solver_residual);
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
            tmp[j] +=
                BM_ABS(solver_data->phi_dt[i * n_sol_variables + j]) * volume;
        }
    }

    for (int i = 0; i < n_sol_variables; ++i)
        tmp[i] *= dt / solver_mesh->global_volume;

    mpi_all_reduce_n(MPIDouble, MPISum, tmp, n_sol_variables, solver_residual);
}