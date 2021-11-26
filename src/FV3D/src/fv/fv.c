/*******************************************************************************
 * @file fv.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "fv_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"

void_update_ft update_function_pointer = NULL;
void_calc_flux_ft calc_flux_function_pointer = NULL;
void_calc_exact_ft calc_exact_function_pointer = NULL;

double *phi_total = NULL;
double *grad_phi_total_x = NULL;
double *grad_phi_total_y = NULL;
double *grad_phi_total_z = NULL;

double *phi_total_left = NULL;
double *phi_total_right = NULL;

double *phi_dt = NULL;
double *flux = NULL;

/*******************************************************************************
 * @brief Define fv
 ******************************************************************************/
void fv_define()
{
    REGISTER_INITIALIZE_ROUTINE(fv_initialize);
    REGISTER_FINALIZE_ROUTINE(fv_finalize);

    reconstruction_define();
    limiter_define();
}

/*******************************************************************************
 * @brief Finalize fv
 ******************************************************************************/
void fv_finalize()
{
    update_function_pointer = NULL;
    calc_flux_function_pointer = NULL;
    calc_exact_function_pointer = NULL;

    DEALLOCATE(phi_total);
    DEALLOCATE(grad_phi_total_x);
    DEALLOCATE(grad_phi_total_y);
    DEALLOCATE(grad_phi_total_z);

    DEALLOCATE(phi_total_left);
    DEALLOCATE(phi_total_right);

    DEALLOCATE(phi_dt);
    DEALLOCATE(flux);
}

/*******************************************************************************
 * @brief Initialize fv
 ******************************************************************************/
void fv_initialize()
{
    Cells_t *cells = global_mesh->cells;
    Boundaries_t *boundaries = global_mesh->boundaries;
    Faces_t *faces = global_mesh->faces;
    int n_local_cells = cells->n_local_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_faces = faces->n_faces;
    int n_tot_variables = all_variables->n_tot_variables;
    int n_sol_variables = all_variables->n_sol_variables;

    phi_total = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    grad_phi_total_x = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    grad_phi_total_y = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    grad_phi_total_z = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));

    phi_total_left = ALLOCATE(sizeof(double) * n_tot_variables * n_faces);
    phi_total_right = ALLOCATE(sizeof(double) * n_tot_variables * n_faces);

    phi_dt = ALLOCATE(sizeof(double) * n_sol_variables * (n_local_cells + n_boundaries));
    flux = ALLOCATE(sizeof(double) * n_sol_variables * n_faces);

    set_solution();
}

/*******************************************************************************
 * @brief The finite volume time derivative
 * @param t
 ******************************************************************************/
void fv_time_derivative(double t)
{
    Cells_t *cells = global_mesh->cells;
    Boundaries_t *boundaries = global_mesh->boundaries;
    Faces_t *faces = global_mesh->faces;
    int n_local_cells = cells->n_local_cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_faces = faces->n_faces;
    int n_sol_variables = all_variables->n_sol_variables;

    update_function_pointer(t);
    reconstruction_function_pointer();

    calc_flux_function_pointer();

    /* the temporal derivative */
    set_value_n(0.0, n_sol_variables * (n_local_cells + n_boundaries), phi_dt);

    for (int i = 0; i < n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];
        double area = faces->area[i];

        for (int j = 0; j < n_sol_variables; ++j)
        {
            phi_dt[fc[0] * n_sol_variables + j] += flux[i * n_sol_variables + j] * area;
            phi_dt[fc[1] * n_sol_variables + j] -= flux[i * n_sol_variables + j] * area;
        }
    }

    for (int i = 0; i < n_domain_cells; ++i)
    {
        double volume = cells->volume[i];

        for (int j = 0; j < n_sol_variables; ++j)
        {
            phi_dt[i * n_sol_variables + j] = -phi_dt[i * n_sol_variables + j] / volume;
        }
    }
}

/*******************************************************************************
 * @brief Set/initialize the solution pointer
 * @param t
 ******************************************************************************/
void set_solution()
{
    Cells_t *cells = global_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_tot_variables = all_variables->n_tot_variables;

    for (int i = 0; i < n_domain_cells; ++i)
        calc_exact_function_pointer(0, 0.0, &cells->x[i * DIM], &phi_total[i * n_tot_variables]);

    update_function_pointer(0.0);
}