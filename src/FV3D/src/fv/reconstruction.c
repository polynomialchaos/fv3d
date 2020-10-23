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
    register_initialize_routine( reconstruction_initialize );
    register_finalize_routine( reconstruction_finalize );

    string_t tmp_opt[] = {"Linear", "First-Order"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    set_parameter( "FV/Reconstruction/reconstruction", ParameterString, &tmp,
        "The reconstruction method", &tmp_opt, tmp_opt_n );
}

void reconstruction_initialize()
{
    get_parameter( "FV/Reconstruction/reconstruction", ParameterString, &reconstruction_name );

    if (is_equal( reconstruction_name, "First-Order" ))
    {
        reconstruction_function_pointer = reconstruction_first_order;
    }
    else if (is_equal( reconstruction_name, "Linear" ))
    {
        reconstruction_function_pointer = reconstruction_linear;
    }
    else
    {
        check_error( 0 );
    }

    int n_tot_variables         = all_variables->n_tot_variables;
    int n_partition_sends       = global_mesh->partition->n_partition_sends;
    int n_partition_receives    = global_mesh->partition->n_partition_receives;

    send_buffer     = allocate( sizeof( double ) * n_tot_variables * n_partition_sends );
    receive_buffer  = allocate( sizeof( double ) * n_tot_variables * n_partition_receives );
}

void reconstruction_finalize()
{
    deallocate( reconstruction_name );
    reconstruction_function_pointer = NULL;

    deallocate( send_buffer );
    deallocate( receive_buffer );
}

void reconstruction_first_order()
{
    int n_local_faces   = global_mesh->faces->n_local_faces;
    int n_tot_variables = all_variables->n_tot_variables;

    calc_gradients();

    for ( int i = 0; i < n_local_faces; i++ )
    {
        int *fc = &global_mesh->faces->cells[i*FACE_CELLS];

        for ( int j = 0; j < n_tot_variables; j++ )
        {
            phi_total_left[i*n_tot_variables+j]     = phi_total[fc[0]*n_tot_variables+j];
            phi_total_right[i*n_tot_variables+j]    = phi_total[fc[1]*n_tot_variables+j];
        }
    }
}

void reconstruction_linear()
{
    int n_internal_faces    = global_mesh->faces->n_internal_faces;
    int n_boundary_faces    = global_mesh->faces->n_boundary_faces;
    int n_tot_variables     = all_variables->n_tot_variables;

    calc_gradients();

    for ( int ii = 0; ii < n_internal_faces; ii++ )
    {
        int i = global_mesh->faces->internal_faces[ii];
        int *fc = &global_mesh->faces->cells[i*FACE_CELLS];

        double *r                  = &global_mesh->faces->dist_cell_1[i*DIM];
        double *grad_phi_total_x_i = &grad_phi_total_x[n_tot_variables*fc[0]];
        double *grad_phi_total_y_i = &grad_phi_total_y[n_tot_variables*fc[0]];
        double *grad_phi_total_z_i = &grad_phi_total_z[n_tot_variables*fc[0]];

        for ( int j = 0; j < n_tot_variables; j++ )
        {
            double slope    = r[0] * grad_phi_total_x_i[i] + r[1] * grad_phi_total_y_i[i] + r[2] * grad_phi_total_z_i[i];
            double lim      = limiter_function_pointer( fc[0], j, slope );

            phi_total_left[i*n_tot_variables+j] = phi_total[fc[0]*n_tot_variables+j] + lim * slope;
        }

        r                   = &global_mesh->faces->dist_cell_2[i*DIM];
        grad_phi_total_x_i  = &grad_phi_total_x[n_tot_variables*fc[1]];
        grad_phi_total_y_i  = &grad_phi_total_y[n_tot_variables*fc[1]];
        grad_phi_total_z_i  = &grad_phi_total_z[n_tot_variables*fc[1]];

        for ( int j = 0; j < n_tot_variables; j++ )
        {
            double slope    = r[0] * grad_phi_total_x_i[i] + r[1] * grad_phi_total_y_i[i] + r[2] * grad_phi_total_z_i[i];
            double lim      = limiter_function_pointer( fc[1], j, slope );

            phi_total_right[i*n_tot_variables+j] = phi_total[fc[1]*n_tot_variables+j] + lim * slope;
        }
    }

    for ( int ii = 0; ii < n_boundary_faces; ii++ )
    {
        int i = global_mesh->faces->boundary_faces[ii];
        int *fc = &global_mesh->faces->cells[i*FACE_CELLS];

        double *r                  = &global_mesh->faces->dist_cell_1[i*DIM];
        double *grad_phi_total_x_i = &grad_phi_total_x[n_tot_variables*fc[0]];
        double *grad_phi_total_y_i = &grad_phi_total_y[n_tot_variables*fc[0]];
        double *grad_phi_total_z_i = &grad_phi_total_z[n_tot_variables*fc[0]];

        for ( int j = 0; j < n_tot_variables; j++ )
        {
            double slope    = r[0] * grad_phi_total_x_i[i] + r[1] * grad_phi_total_y_i[i] + r[2] * grad_phi_total_z_i[i];
            double lim      = limiter_function_pointer( fc[0], j, slope );

            phi_total_left[i*n_tot_variables+j] = phi_total[fc[0]*n_tot_variables+j] + lim * slope;
        }

        for ( int j = 0; j < n_tot_variables; j++ )
        {
            phi_total_right[i*n_tot_variables+j] = phi_total[fc[1]*n_tot_variables+j];
        }
    }
}

void calc_gradients()
{
    if (get_is_parallel()) update_parallel( phi_total );

// !         grad_phi_total_x = 0.0
// !         grad_phi_total_y = 0.0
// !         grad_phi_total_z = 0.0

// !         do i = 1, n_faces
// !             fc = faces(i)%cells

// !             phi_mean = phi_total(:,fc(1)) + faces(i)%lambda * (phi_total(:,fc(2)) - phi_total(:,fc(1)))

// !             grad_phi_total_x(:,fc(1))   = grad_phi_total_x(:,fc(1)) + phi_mean(:) * faces(i)%n(1) * faces(i)%area
// !             grad_phi_total_y(:,fc(1))   = grad_phi_total_y(:,fc(1)) + phi_mean(:) * faces(i)%n(2) * faces(i)%area
// !             grad_phi_total_z(:,fc(1))   = grad_phi_total_z(:,fc(1)) + phi_mean(:) * faces(i)%n(3) * faces(i)%area

// !             grad_phi_total_x(:,fc(2))   = grad_phi_total_x(:,fc(2)) - phi_mean(:) * faces(i)%n(1) * faces(i)%area
// !             grad_phi_total_y(:,fc(2))   = grad_phi_total_y(:,fc(2)) - phi_mean(:) * faces(i)%n(2) * faces(i)%area
// !             grad_phi_total_z(:,fc(2))   = grad_phi_total_z(:,fc(2)) - phi_mean(:) * faces(i)%n(3) * faces(i)%area
// !         end do

    if (get_is_parallel()) update_parallel( grad_phi_total_x );
    if (get_is_parallel()) update_parallel( grad_phi_total_y );
    if (get_is_parallel()) update_parallel( grad_phi_total_z );

// !         do i = 1, n_cells + n_partition_receives
// !             s_volume = 1.0 / cells(i)%volume

// !             grad_phi_total_x(:,i)   = grad_phi_total_x(:,i) * s_volume
// !             grad_phi_total_y(:,i)   = grad_phi_total_y(:,i) * s_volume
// !             grad_phi_total_z(:,i)   = grad_phi_total_z(:,i) * s_volume
// !         end do

    if (update_gradients_function_pointer != NULL) update_gradients_function_pointer();
}

void update_parallel( double *phi_local )
{
    int rank                        = get_rank();
    int n_partitions                = global_mesh->partition->n_partitions;
    int n_partitions_sends          = global_mesh->partition->n_partition_sends;
    int *n_partitions_sends_to      = global_mesh->partition->n_partition_sends_to;
    int n_partition_receives        = global_mesh->partition->n_partition_receives;
    int *n_partition_receives_from  = global_mesh->partition->n_partition_receives_from;
    int n_tot_variables             = all_variables->n_tot_variables;

    for ( int s_rank = 0; s_rank < n_partitions; s_rank++ )
    {

        if (s_rank == rank)
        {
            for ( int r_rank = 0; r_rank < n_partitions; r_rank++ )
            {
                if (n_partitions_sends_to[r_rank] == 0) continue;
                int *partition_sends_to = &global_mesh->partition->partition_sends_to[r_rank*n_partitions_sends];

                for ( int i = 0; i < n_partitions_sends_to[r_rank]; i++ )
                {
                    double *phi_local_i = &phi_local[partition_sends_to[i]*n_tot_variables];

                    for ( int j = 0; j < n_tot_variables; j++ )
                    {
                        send_buffer[i*n_tot_variables+j] = phi_local_i[j];
                    }
                }

                mpi_send( send_buffer, n_partitions_sends_to[r_rank] * n_tot_variables, MPIDouble, r_rank );
            }
        }
        else
        {
            if (n_partition_receives_from[s_rank] == 0) continue;
            int *partition_receives_from = &global_mesh->partition->partition_receives_from[s_rank*n_partition_receives];

            mpi_receive( receive_buffer, n_partition_receives_from[s_rank] * n_tot_variables, MPIDouble, s_rank );

            for ( int i = 0; i < n_partition_receives_from[s_rank]; i++ )
            {
                double *phi_local_i = &phi_local[partition_receives_from[i]*n_tot_variables];

                for ( int j = 0; j < n_tot_variables; j++ )
                {
                    phi_local_i[j] = receive_buffer[i*n_tot_variables+j];
                }
            }
        }
    }
}