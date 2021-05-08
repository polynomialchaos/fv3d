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
void_reconstruction_fp_t reconstruction_function_pointer = NULL;
void_update_gradients_fp_t update_gradients_function_pointer = NULL;

string_t reconstruction_name = NULL;

double *send_buffer = NULL;
double *receive_buffer = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void reconstruction_initialize();
void reconstruction_finalize();

void reconstruction_first_order();
void reconstruction_linear();
void calc_gradients();
void update_parallel();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void reconstruction_define()
{
    register_initialize_routine(reconstruction_initialize);
    register_finalize_routine(reconstruction_finalize);

    string_t tmp_opt[] = {"Linear", "First-Order"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    set_parameter("FV/Reconstruction/reconstruction", ParameterString, &tmp,
                  "The reconstruction method", &tmp_opt, tmp_opt_n);
}

void reconstruction_initialize()
{
    get_parameter("FV/Reconstruction/reconstruction", ParameterString, &reconstruction_name);

    if (is_equal(reconstruction_name, "First-Order"))
    {
        reconstruction_function_pointer = reconstruction_first_order;
    }
    else if (is_equal(reconstruction_name, "Linear"))
    {
        reconstruction_function_pointer = reconstruction_linear;
    }
    else
    {
        check_error(0);
    }

    if (get_is_parallel())
    {
        Partition_t *partition = global_mesh->partition;
        int n_partition_sends = partition->n_partition_sends;
        int n_partition_receives = partition->n_partition_receives;
        int n_tot_variables = all_variables->n_tot_variables;

        send_buffer = allocate(sizeof(double) * n_tot_variables * n_partition_sends);
        receive_buffer = allocate(sizeof(double) * n_tot_variables * n_partition_receives);
    }
}

void reconstruction_finalize()
{
    reconstruction_function_pointer = NULL;
    update_gradients_function_pointer = NULL;

    deallocate(reconstruction_name);
    deallocate(send_buffer);
    deallocate(receive_buffer);
}

void reconstruction_first_order()
{
    Faces_t *faces = global_mesh->faces;
    int n_faces = faces->n_faces;
    int n_tot_variables = all_variables->n_tot_variables;

    calc_gradients();

    for (int i = 0; i < n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            phi_total_left[i * n_tot_variables + j] = phi_total[fc[0] * n_tot_variables + j];
            phi_total_right[i * n_tot_variables + j] = phi_total[fc[1] * n_tot_variables + j];
        }
    }
}

void reconstruction_linear()
{
    Faces_t *faces = global_mesh->faces;
    int n_internal_faces = faces->n_internal_faces;
    int n_boundary_faces = faces->n_boundary_faces;
    int n_tot_variables = all_variables->n_tot_variables;

    calc_gradients();

    for (int ii = 0; ii < n_internal_faces; ++ii)
    {
        int i = faces->internal_faces[ii];
        int *fc = &faces->cells[i * FACE_CELLS];

        double *r = &faces->dist_cell_1[i * DIM];
        double *grad_phi_total_x_i = &grad_phi_total_x[fc[0] * n_tot_variables];
        double *grad_phi_total_y_i = &grad_phi_total_y[fc[0] * n_tot_variables];
        double *grad_phi_total_z_i = &grad_phi_total_z[fc[0] * n_tot_variables];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double slope = r[0] * grad_phi_total_x_i[j] + r[1] * grad_phi_total_y_i[j] + r[2] * grad_phi_total_z_i[j];
            double lim = limiter_function_pointer(fc[0], j, slope);

            phi_total_left[i * n_tot_variables + j] = phi_total[fc[0] * n_tot_variables + j] + lim * slope;
        }

        r = &faces->dist_cell_2[i * DIM];
        grad_phi_total_x_i = &grad_phi_total_x[fc[1] * n_tot_variables];
        grad_phi_total_y_i = &grad_phi_total_y[fc[1] * n_tot_variables];
        grad_phi_total_z_i = &grad_phi_total_z[fc[1] * n_tot_variables];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double slope = r[0] * grad_phi_total_x_i[j] + r[1] * grad_phi_total_y_i[j] + r[2] * grad_phi_total_z_i[j];
            double lim = limiter_function_pointer(fc[1], j, slope);

            phi_total_right[i * n_tot_variables + j] = phi_total[fc[1] * n_tot_variables + j] + lim * slope;
        }
    }

    for (int ii = 0; ii < n_boundary_faces; ++ii)
    {
        int i = faces->boundary_faces[ii];
        int *fc = &faces->cells[i * FACE_CELLS];

        double *r = &faces->dist_cell_1[i * DIM];
        double *grad_phi_total_x_i = &grad_phi_total_x[fc[0] * n_tot_variables];
        double *grad_phi_total_y_i = &grad_phi_total_y[fc[0] * n_tot_variables];
        double *grad_phi_total_z_i = &grad_phi_total_z[fc[0] * n_tot_variables];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double slope = r[0] * grad_phi_total_x_i[j] + r[1] * grad_phi_total_y_i[j] + r[2] * grad_phi_total_z_i[j];
            double lim = limiter_function_pointer(fc[0], j, slope);

            phi_total_left[i * n_tot_variables + j] = phi_total[fc[0] * n_tot_variables + j] + lim * slope;
        }

        for (int j = 0; j < n_tot_variables; ++j)
        {
            phi_total_right[i * n_tot_variables + j] = phi_total[fc[1] * n_tot_variables + j];
        }
    }
}

void calc_gradients()
{
    Cells_t *cells = global_mesh->cells;
    Boundaries_t *boundaries = global_mesh->boundaries;
    Faces_t *faces = global_mesh->faces;
    int n_local_cells = cells->n_local_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_faces = faces->n_faces;
    int n_tot_variables = all_variables->n_tot_variables;

    if (get_is_parallel())
        update_parallel(phi_total);

    set_value_n(0.0, grad_phi_total_x, n_tot_variables * (n_local_cells + n_boundaries));
    set_value_n(0.0, grad_phi_total_y, n_tot_variables * (n_local_cells + n_boundaries));
    set_value_n(0.0, grad_phi_total_z, n_tot_variables * (n_local_cells + n_boundaries));

    for (int i = 0; i < n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];
        double *n = &faces->n[i * DIM];
        double area = faces->area[i];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double phi_mean = phi_total[fc[0] * n_tot_variables + j] + faces->lambda[i] *
                                                                           (phi_total[fc[1] * n_tot_variables + j] - phi_total[fc[0] * n_tot_variables + j]);

            grad_phi_total_x[fc[0] * n_tot_variables + j] += phi_mean * n[0] * area;
            grad_phi_total_y[fc[0] * n_tot_variables + j] += phi_mean * n[1] * area;
            grad_phi_total_z[fc[0] * n_tot_variables + j] += phi_mean * n[2] * area;

            grad_phi_total_x[fc[1] * n_tot_variables + j] -= phi_mean * n[0] * area;
            grad_phi_total_y[fc[1] * n_tot_variables + j] -= phi_mean * n[1] * area;
            grad_phi_total_z[fc[1] * n_tot_variables + j] -= phi_mean * n[2] * area;
        }
    }

    if (get_is_parallel())
        update_parallel(grad_phi_total_x);
    if (get_is_parallel())
        update_parallel(grad_phi_total_y);
    if (get_is_parallel())
        update_parallel(grad_phi_total_z);

    for (int i = 0; i < n_local_cells; ++i)
    {
        double s_volume = 1.0 / cells->volume[i];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            grad_phi_total_x[i * n_tot_variables + j] *= s_volume;
            grad_phi_total_y[i * n_tot_variables + j] *= s_volume;
            grad_phi_total_z[i * n_tot_variables + j] *= s_volume;
        }
    }

    update_gradients_function_pointer();
}

void update_parallel(double *phi_local)
{
    Partition_t *partition = global_mesh->partition;
    int rank = get_rank();
    int n_partitions = partition->n_partitions;
    int n_partitions_sends = partition->n_partition_sends;
    int *n_partitions_sends_to = partition->n_partition_sends_to;
    int n_partition_receives = partition->n_partition_receives;
    int *n_partition_receives_from = partition->n_partition_receives_from;
    int n_tot_variables = all_variables->n_tot_variables;

    for (int s_rank = 0; s_rank < n_partitions; ++s_rank)
    {
        if (s_rank == rank)
        {
            for (int r_rank = 0; r_rank < n_partitions; ++r_rank)
            {
                if (n_partitions_sends_to[r_rank] == 0)
                    continue;
                int *partition_sends_to = &partition->partition_sends_to[r_rank * n_partitions_sends];

                for (int i = 0; i < n_partitions_sends_to[r_rank]; ++i)
                {
                    double *phi_local_i = &phi_local[partition_sends_to[i] * n_tot_variables];

                    for (int j = 0; j < n_tot_variables; ++j)
                    {
                        send_buffer[i * n_tot_variables + j] = phi_local_i[j];
                    }
                }

                mpi_send(send_buffer, n_partitions_sends_to[r_rank] * n_tot_variables, MPIDouble, r_rank);
            }
        }
        else
        {
            if (n_partition_receives_from[s_rank] == 0)
                continue;
            int *partition_receives_from = &partition->partition_receives_from[s_rank * n_partition_receives];

            mpi_receive(receive_buffer, n_partition_receives_from[s_rank] * n_tot_variables, MPIDouble, s_rank);

            for (int i = 0; i < n_partition_receives_from[s_rank]; ++i)
            {
                double *phi_local_i = &phi_local[partition_receives_from[i] * n_tot_variables];

                for (int j = 0; j < n_tot_variables; ++j)
                {
                    phi_local_i[j] = receive_buffer[i * n_tot_variables + j];
                }
            }
        }
    }
}