//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "fv_private.h"
#include "reconstruction_private.h"
#include "limiter_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------

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
    register_initialize_routine( fv_initialize );
    register_finalize_routine( fv_finalize );

    reconstruction_define();
    limiter_define();
}

void fv_initialize()
{
        //     ! generate total arrays
        // allocate( phi_total(n_tot_variables,n_cells+n_partition_receives+n_boundaries) ); phi_total = 0.0

        // allocate( grad_phi_total_x(n_tot_variables,n_cells+n_partition_receives+n_boundaries) ); grad_phi_total_x = 0.0
        // allocate( grad_phi_total_y(n_tot_variables,n_cells+n_partition_receives+n_boundaries) ); grad_phi_total_y = 0.0
        // allocate( grad_phi_total_z(n_tot_variables,n_cells+n_partition_receives+n_boundaries) ); grad_phi_total_z = 0.0

        // allocate( phi_dt(n_variables,n_cells+n_partition_receives+n_boundaries) ); phi_dt = 0.0

        // phi => phi_total(:n_variables,:)

        // ! generate solver arrays
        // allocate( phi_total_left(n_tot_variables,n_faces) ); phi_total_left = 0.0
        // allocate( phi_total_right(n_tot_variables,n_faces) ); phi_total_right = 0.0
        // allocate( flux(n_variables,n_faces) ); flux = 0.0

        // call initialize_solution()
}

void fv_finalize()
{
//             _DEALLOCATE( phi_total )

//         _DEALLOCATE( grad_phi_total_x )
//         _DEALLOCATE( grad_phi_total_y )
//         _DEALLOCATE( grad_phi_total_z )

//         _DEALLOCATE( phi_total_left )
//         _DEALLOCATE( phi_total_left )

//         nullify( phi )
//         _DEALLOCATE( phi_dt )
//         _DEALLOCATE( flux )
}

void fv_time_derivative( double t )
{
        // call update_routine( t )
        // call reconstruction_routine()

        // call calc_flux_routine()

        // ! the temporal derivative
        // phi_dt = 0.0
        // do i = 1, n_faces
        //     fc = faces(i)%cells

        //     phi_dt(:,fc(1)) = phi_dt(:,fc(1)) + flux(:,i) * faces(i)%area
        //     phi_dt(:,fc(2)) = phi_dt(:,fc(2)) - flux(:,i) * faces(i)%area
        // end do

        // do i = 1, n_cells
        //     phi_dt(:,i) = -phi_dt(:,i) / cells(i)%volume
        // end do
}

void set_solution()
{
        // do i = 1, n_cells
        //     call exact_func_routine( 0, 0.0, cells(i)%x, phi_total(:,i) )
        // end do

        // call update_routine( 0.0 )
}