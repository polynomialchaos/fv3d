/*******************************************************************************
 * @file finite_volume.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "equation/equation_private.h"
#include "finite_volume_private.h"

void_calc_flux_ft calc_flux_function_pointer = NULL;
void_calc_exact_ft calc_exact_function_pointer = NULL;
void_update_ft update_function_pointer = NULL;
void_update_gradients_ft update_gradients_function_pointer = NULL;

data_t *solver_data = NULL;

/*******************************************************************************
 * @brief Allocate the solver data
 ******************************************************************************/
data_t *allocate_data()
{
    const Cells_t *cells = solver_mesh->cells;
    const Boundaries_t *boundaries = solver_mesh->boundaries;
    const Faces_t *faces = solver_mesh->faces;

    const int n_local_cells = cells->n_local_cells;
    const int n_boundaries = boundaries->n_boundaries;
    const int n_elements = n_local_cells + n_boundaries;
    const int n_faces = faces->n_faces;

    const int n_tot_variables = solver_variables->n_tot_variables;
    const int n_sol_variables = solver_variables->n_sol_variables;

    data_t *tmp = ALLOCATE(sizeof(data_t));

    tmp->flux = ALLOCATE(sizeof(double) * n_sol_variables * n_faces);
    tmp->phi_total_left = ALLOCATE(sizeof(double) * n_tot_variables * n_faces);
    tmp->phi_total_right = ALLOCATE(sizeof(double) * n_tot_variables * n_faces);

    tmp->grad_phi_total_x =
        ALLOCATE(sizeof(double) * n_tot_variables * n_elements);
    tmp->grad_phi_total_y =
        ALLOCATE(sizeof(double) * n_tot_variables * n_elements);
    tmp->grad_phi_total_z =
        ALLOCATE(sizeof(double) * n_tot_variables * n_elements);
    tmp->phi_dt = ALLOCATE(sizeof(double) * n_sol_variables * n_elements);
    tmp->phi_total = ALLOCATE(sizeof(double) * n_tot_variables * n_elements);

    return tmp;
}

/*******************************************************************************
 * @brief Deallocate the solver data
 ******************************************************************************/
void deallocate_data(data_t *this)
{
    DEALLOCATE(this->flux);
    DEALLOCATE(this->phi_total_left);
    DEALLOCATE(this->phi_total_right);

    DEALLOCATE(this->grad_phi_total_x);
    DEALLOCATE(this->grad_phi_total_y);
    DEALLOCATE(this->grad_phi_total_z);
    DEALLOCATE(this->phi_dt);
    DEALLOCATE(this->phi_total);
}

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
    set_value_n(0.0, n_sol_variables * (n_local_cells + n_boundaries), solver_data->phi_dt);

    for (int i = 0; i < n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];
        double area = faces->area[i];

        for (int j = 0; j < n_sol_variables; ++j)
        {
            solver_data->phi_dt[fc[0] * n_sol_variables + j] += solver_data->flux[i * n_sol_variables + j] * area;
            solver_data->phi_dt[fc[1] * n_sol_variables + j] -= solver_data->flux[i * n_sol_variables + j] * area;
        }
    }

    for (int i = 0; i < n_domain_cells; ++i)
    {
        double volume = cells->volume[i];

        for (int j = 0; j < n_sol_variables; ++j)
        {
            solver_data->phi_dt[i * n_sol_variables + j] = -solver_data->phi_dt[i * n_sol_variables + j] / volume;
        }
    }
}

/*******************************************************************************
 * @brief Free finite_volume
 ******************************************************************************/
void free_finite_volume()
{
    update_function_pointer = NULL;
    update_gradients_function_pointer = NULL;
    calc_flux_function_pointer = NULL;
    calc_exact_function_pointer = NULL;

    deallocate_data(solver_data);
    DEALLOCATE(solver_data);
}

/*******************************************************************************
 * @brief Initialize finite_volume
 ******************************************************************************/
void init_finite_volume()
{
    solver_data = allocate_data();
    set_solution();
}

/*******************************************************************************
 * @brief Set the flux calculation routine
 * @param fun_ptr
 ******************************************************************************/
void set_calc_flux(void_calc_flux_ft fun_ptr)
{
    calc_flux_function_pointer = fun_ptr;
}

/*******************************************************************************
 * @brief Set the exact function routine
 * @param fun_ptr
 ******************************************************************************/
void set_exact_function(void_calc_exact_ft fun_ptr)
{
    calc_exact_function_pointer = fun_ptr;
}

/*******************************************************************************
 * @brief Set the update routine
 * @param fun_ptr
 ******************************************************************************/
void set_update(void_update_ft fun_ptr)
{
    update_function_pointer = fun_ptr;
}

/*******************************************************************************
 * @brief Set the update gradients routine
 * @param fun_ptr
 ******************************************************************************/
void set_update_gradients(void_update_gradients_ft fun_ptr)
{
    update_gradients_function_pointer = fun_ptr;
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
        calc_exact_function_pointer(0, 0.0, &cells->x[i * DIM], &solver_data->phi_total[i * n_tot_variables]);

    update_function_pointer(0.0);
}