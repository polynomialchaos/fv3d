/*******************************************************************************
 * @file reconstruction.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "finite_volume_private.h"

void_reconstruction_ft reconstruction_function_pointer = NULL;
void_update_gradients_ft update_gradients_function_pointer = NULL;

double *send_buffer = NULL;
double *receive_buffer = NULL;

/*******************************************************************************
 * @brief Calculate (reconstruct) gradients
 ******************************************************************************/
void calc_gradients()
{
    Cells_t *cells = solver_mesh->cells;
    Boundaries_t *boundaries = solver_mesh->boundaries;
    Faces_t *faces = solver_mesh->faces;
    int n_local_cells = cells->n_local_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_faces = faces->n_faces;
    int n_tot_variables = solver_variables->n_tot_variables;

    if (is_parallel())
        update_parallel(solver_phi_total);

    set_value_n(0.0, n_tot_variables * (n_local_cells + n_boundaries), solver_grad_phi_total_x);
    set_value_n(0.0, n_tot_variables * (n_local_cells + n_boundaries), solver_grad_phi_total_y);
    set_value_n(0.0, n_tot_variables * (n_local_cells + n_boundaries), solver_grad_phi_total_z);

    for (int i = 0; i < n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];
        double *n = &faces->n[i * DIM];
        double area = faces->area[i];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double phi_mean = solver_phi_total[fc[0] * n_tot_variables + j] + faces->lambda[i] *
                                                                           (solver_phi_total[fc[1] * n_tot_variables + j] - solver_phi_total[fc[0] * n_tot_variables + j]);

            solver_grad_phi_total_x[fc[0] * n_tot_variables + j] += phi_mean * n[0] * area;
            solver_grad_phi_total_y[fc[0] * n_tot_variables + j] += phi_mean * n[1] * area;
            solver_grad_phi_total_z[fc[0] * n_tot_variables + j] += phi_mean * n[2] * area;

            solver_grad_phi_total_x[fc[1] * n_tot_variables + j] -= phi_mean * n[0] * area;
            solver_grad_phi_total_y[fc[1] * n_tot_variables + j] -= phi_mean * n[1] * area;
            solver_grad_phi_total_z[fc[1] * n_tot_variables + j] -= phi_mean * n[2] * area;
        }
    }

    if (is_parallel())
        update_parallel(solver_grad_phi_total_x);
    if (is_parallel())
        update_parallel(solver_grad_phi_total_y);
    if (is_parallel())
        update_parallel(solver_grad_phi_total_z);

    for (int i = 0; i < n_local_cells; ++i)
    {
        double s_volume = 1.0 / cells->volume[i];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            solver_grad_phi_total_x[i * n_tot_variables + j] *= s_volume;
            solver_grad_phi_total_y[i * n_tot_variables + j] *= s_volume;
            solver_grad_phi_total_z[i * n_tot_variables + j] *= s_volume;
        }
    }

    update_gradients_function_pointer();
}

/*******************************************************************************
 * @brief Free reconstruction
 ******************************************************************************/
void free_reconstruction()
{
    reconstruction_function_pointer = NULL;
    update_gradients_function_pointer = NULL;

    DEALLOCATE(send_buffer);
    DEALLOCATE(receive_buffer);
}

/*******************************************************************************
 * @brief Initialize reconstruction
 ******************************************************************************/
void init_reconstruction(reconstruction_type_t reconstruction_type)
{
    switch (reconstruction_type)
    {
    case FirstOrder:
        reconstruction_function_pointer = reconstruction_first_order;
        break;
    case Linear:
        reconstruction_function_pointer = reconstruction_linear;
        break;
    default:
        CHECK_EXPRESSION(0);
        break;
    }

    if (is_parallel())
    {
        Partition_t *partition = solver_mesh->partition;
        int n_partition_sends = partition->n_partition_sends;
        int n_partition_receives = partition->n_partition_receives;
        int n_tot_variables = solver_variables->n_tot_variables;

        send_buffer = ALLOCATE(sizeof(double) * n_tot_variables * n_partition_sends);
        receive_buffer = ALLOCATE(sizeof(double) * n_tot_variables * n_partition_receives);
    }
}

/*******************************************************************************
 * @brief First-order reconstruction
 ******************************************************************************/
void reconstruction_first_order()
{
    Faces_t *faces = solver_mesh->faces;
    int n_faces = faces->n_faces;
    int n_tot_variables = solver_variables->n_tot_variables;

    calc_gradients();

    for (int i = 0; i < n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            solver_phi_total_left[i * n_tot_variables + j] = solver_phi_total[fc[0] * n_tot_variables + j];
            solver_phi_total_right[i * n_tot_variables + j] = solver_phi_total[fc[1] * n_tot_variables + j];
        }
    }
}

/*******************************************************************************
 * @brief Second-order (linear) reconstruction
 ******************************************************************************/
void reconstruction_linear()
{
    Faces_t *faces = solver_mesh->faces;
    int n_internal_faces = faces->n_internal_faces;
    int n_boundary_faces = faces->n_boundary_faces;
    int n_tot_variables = solver_variables->n_tot_variables;

    calc_gradients();

    for (int ii = 0; ii < n_internal_faces; ++ii)
    {
        int i = faces->internal_faces[ii];
        int *fc = &faces->cells[i * FACE_CELLS];

        double *r = &faces->dist_cell_1[i * DIM];
        double *grad_phi_total_x_i = &solver_grad_phi_total_x[fc[0] * n_tot_variables];
        double *grad_phi_total_y_i = &solver_grad_phi_total_y[fc[0] * n_tot_variables];
        double *grad_phi_total_z_i = &solver_grad_phi_total_z[fc[0] * n_tot_variables];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double slope = r[0] * grad_phi_total_x_i[j] + r[1] * grad_phi_total_y_i[j] + r[2] * grad_phi_total_z_i[j];
            double lim = limiter_function_pointer(fc[0], j, slope);

            solver_phi_total_left[i * n_tot_variables + j] = solver_phi_total[fc[0] * n_tot_variables + j] + lim * slope;
        }

        r = &faces->dist_cell_2[i * DIM];
        grad_phi_total_x_i = &solver_grad_phi_total_x[fc[1] * n_tot_variables];
        grad_phi_total_y_i = &solver_grad_phi_total_y[fc[1] * n_tot_variables];
        grad_phi_total_z_i = &solver_grad_phi_total_z[fc[1] * n_tot_variables];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double slope = r[0] * grad_phi_total_x_i[j] + r[1] * grad_phi_total_y_i[j] + r[2] * grad_phi_total_z_i[j];
            double lim = limiter_function_pointer(fc[1], j, slope);

            solver_phi_total_right[i * n_tot_variables + j] = solver_phi_total[fc[1] * n_tot_variables + j] + lim * slope;
        }
    }

    for (int ii = 0; ii < n_boundary_faces; ++ii)
    {
        int i = faces->boundary_faces[ii];
        int *fc = &faces->cells[i * FACE_CELLS];

        double *r = &faces->dist_cell_1[i * DIM];
        double *grad_phi_total_x_i = &solver_grad_phi_total_x[fc[0] * n_tot_variables];
        double *grad_phi_total_y_i = &solver_grad_phi_total_y[fc[0] * n_tot_variables];
        double *grad_phi_total_z_i = &solver_grad_phi_total_z[fc[0] * n_tot_variables];

        for (int j = 0; j < n_tot_variables; ++j)
        {
            double slope = r[0] * grad_phi_total_x_i[j] + r[1] * grad_phi_total_y_i[j] + r[2] * grad_phi_total_z_i[j];
            double lim = limiter_function_pointer(fc[0], j, slope);

            solver_phi_total_left[i * n_tot_variables + j] = solver_phi_total[fc[0] * n_tot_variables + j] + lim * slope;
        }

        for (int j = 0; j < n_tot_variables; ++j)
        {
            solver_phi_total_right[i * n_tot_variables + j] = solver_phi_total[fc[1] * n_tot_variables + j];
        }
    }
}

/*******************************************************************************
 * @brief Update the given variable across all domains
 * @param phi_local
 ******************************************************************************/
void update_parallel(double *phi_local)
{
    Partition_t *partition = solver_mesh->partition;
    int rank = get_rank_number();
    int n_partitions = partition->n_partitions;
    int n_partitions_sends = partition->n_partition_sends;
    int *n_partitions_sends_to = partition->n_partition_sends_to;
    int n_partition_receives = partition->n_partition_receives;
    int *n_partition_receives_from = partition->n_partition_receives_from;
    int n_tot_variables = solver_variables->n_tot_variables;

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

                mpi_send(MPIDouble, r_rank, send_buffer, n_partitions_sends_to[r_rank] * n_tot_variables);
            }
        }
        else
        {
            if (n_partition_receives_from[s_rank] == 0)
                continue;
            int *partition_receives_from = &partition->partition_receives_from[s_rank * n_partition_receives];

            mpi_receive(MPIDouble, s_rank, n_partition_receives_from[s_rank] * n_tot_variables, receive_buffer);

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