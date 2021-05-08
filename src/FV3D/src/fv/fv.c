//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "fv_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
void_update_fp_t update_function_pointer = NULL;
void_calc_flux_fp_t calc_flux_function_pointer = NULL;
void_calc_exact_fp_t calc_exact_function_pointer = NULL;

double *phi_total = NULL;
double *grad_phi_total_x = NULL;
double *grad_phi_total_y = NULL;
double *grad_phi_total_z = NULL;

double *phi_total_left = NULL;
double *phi_total_right = NULL;

double *phi_dt = NULL;
double *flux = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void fv_initialize();
void fv_finalize();

void set_solution();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void fv_define()
{
    register_initialize_routine(fv_initialize);
    register_finalize_routine(fv_finalize);

    reconstruction_define();
    limiter_define();
}

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

    phi_total = allocate(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    grad_phi_total_x = allocate(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    grad_phi_total_y = allocate(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));
    grad_phi_total_z = allocate(sizeof(double) * n_tot_variables * (n_local_cells + n_boundaries));

    phi_total_left = allocate(sizeof(double) * n_tot_variables * n_faces);
    phi_total_right = allocate(sizeof(double) * n_tot_variables * n_faces);

    phi_dt = allocate(sizeof(double) * n_sol_variables * (n_local_cells + n_boundaries));
    flux = allocate(sizeof(double) * n_sol_variables * n_faces);

    set_solution();
}

void fv_finalize()
{
    update_function_pointer = NULL;
    calc_flux_function_pointer = NULL;
    calc_exact_function_pointer = NULL;

    deallocate(phi_total);
    deallocate(grad_phi_total_x);
    deallocate(grad_phi_total_y);
    deallocate(grad_phi_total_z);

    deallocate(phi_total_left);
    deallocate(phi_total_right);

    deallocate(phi_dt);
    deallocate(flux);
}

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

    // the temporal derivative
    set_value_n(0.0, phi_dt, n_sol_variables * (n_local_cells + n_boundaries));

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

void set_solution()
{
    Cells_t *cells = global_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_tot_variables = all_variables->n_tot_variables;

    for (int i = 0; i < n_domain_cells; ++i)
        calc_exact_function_pointer(0, 0.0, &cells->x[i * DIM], &phi_total[i * n_tot_variables]);

    update_function_pointer(0.0);
}