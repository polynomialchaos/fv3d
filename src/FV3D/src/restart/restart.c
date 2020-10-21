//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "restart_private.h"
#include "output/output_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int use_restart = 0;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void restart_initialize();
void restart_finalize();
void read_restart_data();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void restart_define()
{
    register_initialize_routine( restart_initialize );
    register_finalize_routine( restart_finalize );

    string_t tmp = "untitled_mesh.h5";
    set_parameter( "Restart/use_restart", ParameterBool, &use_restart, "The flag to start from restart", NULL, 0 );
}

void restart_initialize()
{
    get_parameter( "Restart/use_restart", ParameterBool, &use_restart );

    if (use_restart == 1)
        read_restart_data();
    else
        create_file_header();
}

void restart_finalize()
{
        //     _DEALLOCATE( phi_restart )
        // _DEALLOCATE( phi_dt_restart )
        // _DEALLOCATE( phi_old_restart )
}

void read_restart_data()
{
        // file_id = open_hdf5_file( output_file )
        // last_id = open_hdf5_group( file_id, 'SOLUTION' )

        //     call get_hdf5_attribute( last_id, 'iter', iter_restart )
        //     call get_hdf5_attribute( last_id, 'time', time_restart )

        //     allocate( phi_restart(n_variables,n_cells+n_partition_receives+n_boundaries) )
        //     allocate( phi_total_restart(n_tot_variables,n_cells+n_partition_receives+n_boundaries) )
        //     allocate( phi_dt_restart(n_variables,n_cells+n_partition_receives+n_boundaries) )

        //     if( get_is_parallel() ) then
        //         call get_hdf5_dataset( last_id, 'phi', phi_restart(:,:n_cells), stride=partition_cells )
        //         call get_hdf5_dataset( last_id, 'phi_total', phi_total_restart(:,:n_cells), stride=partition_cells )
        //         call get_hdf5_dataset( last_id, 'phi_dt', phi_dt_restart(:,:n_cells), stride=partition_cells )

        //         if( exists_hdf5_attribute( last_id, 'n_stages' ) ) then
        //             call get_hdf5_attribute( last_id, 'n_stages', n_old_stages )

        //             allocate( phi_old_restart(n_variables,n_cells,n_old_stages) )
        //             do i_stage = 1, n_old_stages
        //                 call get_hdf5_dataset( last_id, 'phi_old:' // set_string( i_stage ), &
        //                     phi_old_restart(:,:,i_stage), stride=partition_cells )
        //             end do
        //         end if
        //     else
        //         call get_hdf5_dataset( last_id, 'phi', phi_restart(:,:n_cells) )
        //         call get_hdf5_dataset( last_id, 'phi_total', phi_total_restart(:,:n_cells) )
        //         call get_hdf5_dataset( last_id, 'phi_dt', phi_dt_restart(:,:n_cells) )

        //         if( exists_hdf5_attribute( last_id, 'n_stages' ) ) then
        //             call get_hdf5_attribute( last_id, 'n_stages', n_old_stages )

        //             allocate( phi_old_restart(n_variables,n_cells,n_old_stages) )
        //             do i_stage = 1, n_old_stages
        //                 call get_hdf5_dataset( last_id, 'phi_old:' // set_string( i_stage ), phi_old_restart(:,:,i_stage) )
        //             end do
        //         end if
        //     end if

        // call close_hdf5_group( last_id, 'SOLUTION' )
        // call close_hdf5_file( file_id, output_file )

        // phi         = phi_restart
        // phi_total   = phi_total_restart
        // phi_dt      = phi_dt_restart

        // if( n_old_stages .lt. n_stages ) &
        //     call add_error( __LINE__, __FILE__, &
        //         'Restart file does not provide enough old stage solutions (req:' // set_string( n_stages ) // &
        //         ',got:' // set_string( n_old_stages ) // ')!' )

        // do i_stage = 1, n_stages
        //     phi_old(:,:,i_stage) = phi_old_restart(:,:,i_stage)
        // end do

        // _DEALLOCATE( phi_restart )
        // _DEALLOCATE( phi_total_restart )
        // _DEALLOCATE( phi_dt_restart )
        // _DEALLOCATE( phi_old_restart )
}