/*******************************************************************************
 * @file finite_volume.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "finite_volume_private.h"

void_update_ft update_function_pointer = NULL;
void_calc_flux_ft calc_flux_function_pointer = NULL;
void_calc_exact_ft calc_exact_function_pointer = NULL;

double *solver_phi_total = NULL;
double *solver_grad_phi_total_x = NULL;
double *solver_grad_phi_total_y = NULL;
double *solver_grad_phi_total_z = NULL;

double *solver_phi_total_left = NULL;
double *solver_phi_total_right = NULL;

double *solver_phi_dt = NULL;
double *solver_flux = NULL;

/*******************************************************************************
 * @brief The finite volume time derivative
 * @param t
 ******************************************************************************/
void finite_volume_time_derivative(double t)
{
    Cells_t *cells = solver_mesh->cells;
    Boundaries_t *boundaries = solver_mesh->boundaries;
    Faces_t *faces = solver_mesh->faces;
    int n_local_cells = cells->n_local_cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_faces = faces->n_faces;
    int n_sol_variables = solver_variables->n_sol_variables;

    update_function_pointer(t);
    reconstruction_function_pointer();

    calc_flux_function_pointer();

    /* the temporal derivative */
    set_value_n(0.0, n_sol_variables * (n_local_cells + n_boundaries), solver_phi_dt);

    for (int i = 0; i < n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];
        double area = faces->area[i];

        for (int j = 0; j < n_sol_variables; ++j)
        {
            solver_phi_dt[fc[0] * n_sol_variables + j] += solver_flux[i * n_sol_variables + j] * area;
            solver_phi_dt[fc[1] * n_sol_variables + j] -= solver_flux[i * n_sol_variables + j] * area;
        }
    }

    for (int i = 0; i < n_domain_cells; ++i)
    {
        double volume = cells->volume[i];

        for (int j = 0; j < n_sol_variables; ++j)
        {
            solver_phi_dt[i * n_sol_variables + j] = -solver_phi_dt[i * n_sol_variables + j] / volume;
        }
    }
}

/*******************************************************************************
 * @brief Free finite_volume
 ******************************************************************************/
void free_finite_volume()
{
    update_function_pointer = NULL;
    calc_flux_function_pointer = NULL;
    calc_exact_function_pointer = NULL;

    DEALLOCATE(solver_phi_total);
    DEALLOCATE(solver_grad_phi_total_x);
    DEALLOCATE(solver_grad_phi_total_y);
    DEALLOCATE(solver_grad_phi_total_z);

    DEALLOCATE(solver_phi_total_left);
    DEALLOCATE(solver_phi_total_right);

    DEALLOCATE(solver_phi_dt);
    DEALLOCATE(solver_flux);
}


/*******************************************************************************
 * @brief Initialize finite_volume
 ******************************************************************************/
void init_finite_volume()
{
    Cells_t *cells = solver_mesh->cells;
    Boundaries_t *boundaries = solver_mesh->boundaries;
    Faces_t *faces = solver_mesh->faces;
    int n_local_cells = cells->n_local_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_faces = faces->n_faces;
    int n_tot_variables = solver_variables->n_tot_variables;
    int n_sol_variables = solver_variables->n_sol_variables;

    solver_phi_total = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    solver_grad_phi_total_x = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    solver_grad_phi_total_y = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    solver_grad_phi_total_z = ALLOCATE(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));

    solver_phi_total_left = ALLOCATE(sizeof(double) * n_tot_variables * n_faces);
    solver_phi_total_right = ALLOCATE(sizeof(double) * n_tot_variables * n_faces);

    solver_phi_dt = ALLOCATE(sizeof(double) * n_sol_variables * (n_local_cells + n_boundaries));
    solver_flux = ALLOCATE(sizeof(double) * n_sol_variables * n_faces);

    set_solution();
}



/*******************************************************************************
 * @brief Set/initialize the solution pointer
 * @param t
 ******************************************************************************/
void set_solution()
{
    Cells_t *cells = solver_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_tot_variables = solver_variables->n_tot_variables;

    for (int i = 0; i < n_domain_cells; ++i)
        calc_exact_function_pointer(0, 0.0, &cells->x[i * DIM], &solver_phi_total[i * n_tot_variables]);

    update_function_pointer(0.0);
}