//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "reconstruction_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
void_reconstruction_fp_t reconstruction_routine = NULL;
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
        reconstruction_routine = reconstruction_first_order;
    }
    else if (is_equal( reconstruction_name, "Linear" ))
    {
        reconstruction_routine = reconstruction_linear;
    }
    else
    {
        check_error( 0 );
    }

    // ! allocate( send_buffer(n_tot_variables * n_partition_sends) )
    // ! allocate( receive_buffer(n_tot_variables * n_partition_receives) )
}

void reconstruction_finalize()
{
    deallocate( reconstruction_name );
    reconstruction_routine = NULL;

    send_buffer = NULL;
    receive_buffer = NULL;
}

void reconstruction_first_order()
{
//     !         call calc_gradients()

// !         do i = 1, n_faces
// !             fc = faces(i)%cells

// !             ! reconstruction of variables
// !             phi_total_left(:,i)     = phi_total(:,fc(1))
// !             phi_total_right(:,i)    = phi_total(:,fc(2))
// !         end do
}

void reconstruction_linear()
{
//     !         call calc_gradients()

// !         do ii = 1, n_internal_faces
// !             i  = internal_faces(ii)
// !             fc = faces(i)%cells

// !             r       => dist_cell_1(:,i)
// !             slope   = r(1) * grad_phi_total_x(:,fc(1)) + r(2) * grad_phi_total_y(:,fc(1)) + r(3) * grad_phi_total_z(:,fc(1))
// !             call limiter_routine( fc(1), slope, lim )

// !             phi_total_left(:,i)     = phi_total(:,fc(1)) + lim * slope

// !             r       => dist_cell_2(:,i)
// !             slope   = r(1) * grad_phi_total_x(:,fc(2)) + r(2) * grad_phi_total_y(:,fc(2)) + r(3) * grad_phi_total_z(:,fc(2))
// !             call limiter_routine( fc(2), slope, lim )

// !             phi_total_right(:,i)    = phi_total(:,fc(2)) + lim * slope
// !         end do

// !         do ii = 1, n_boundary_faces
// !             i  = boundary_faces(ii)
// !             fc = faces(i)%cells

// !             r       => dist_cell_1(:,i)
// !             slope   = r(1) * grad_phi_total_x(:,fc(1)) + r(2) * grad_phi_total_y(:,fc(1)) + r(3) * grad_phi_total_z(:,fc(1))
// !             call limiter_routine( fc(1), slope, lim )

// !             phi_total_left(:,i)     = phi_total(:,fc(1)) + lim * slope

// !             phi_total_right(:,i)    = phi_total(:,fc(2))
// !         end do
}

void calc_gradients()
{
// !         if( get_is_parallel() ) call update_parallel( phi_total )

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

// !         if( get_is_parallel() ) call update_parallel( grad_phi_total_x )
// !         if( get_is_parallel() ) call update_parallel( grad_phi_total_y )
// !         if( get_is_parallel() ) call update_parallel( grad_phi_total_z )

// !         do i = 1, n_cells + n_partition_receives
// !             s_volume = 1.0 / cells(i)%volume

// !             grad_phi_total_x(:,i)   = grad_phi_total_x(:,i) * s_volume
// !             grad_phi_total_y(:,i)   = grad_phi_total_y(:,i) * s_volume
// !             grad_phi_total_z(:,i)   = grad_phi_total_z(:,i) * s_volume
// !         end do

// !         call update_gradients_routine()
}

void update_parallel()
{
// !         t_rank = get_i_rank()

// !         do s_rank = 1, n_partitions
// !             if( t_rank .eq. s_rank-1 ) then
// !                 do r_rank = 1, n_partitions
// !                     if( n_partition_sends_to(r_rank) .eq. 0 ) cycle

// !                     do i = 1, n_partition_sends_to(r_rank)
// !                         j = (i - 1) * n_tot_variables
// !                         send_buffer(j+1:j+n_tot_variables) = phi_local(:,partition_sends_to(i,r_rank))
// !                     end do

// !                     call send_mpi( send_buffer(:j+n_tot_variables), r_rank-1 )
// !                 end do
// !             else if( n_partition_receives_from(s_rank) .gt. 0 ) then
// !                 k = (n_partition_receives_from(s_rank) - 1) * n_tot_variables
// !                 call receive_mpi( receive_buffer(:k+n_tot_variables), s_rank-1 )

// !                 do i = 1, n_partition_receives_from(s_rank)
// !                     j = (i - 1) * n_tot_variables
// !                     phi_local(:,partition_receives_from(i,s_rank)) = receive_buffer(j+1:j+n_tot_variables)
// !                 end do
// !             end if
// !         end do
}