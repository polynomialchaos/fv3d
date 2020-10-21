//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "output_private.h"
#include "main_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int i_output_data = -1;
int do_output_data = 0;
string_t output_file = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void output_initialize();
void output_finalize();
void write_output_solution( int iter, double time, string_t iter_string );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void output_define()
{
    register_initialize_routine( output_initialize );
    register_finalize_routine( output_finalize );

    set_parameter( "Output/i_output_data", ParameterDigit, &i_output_data,
        "The output file frequency  (-1 ... first/solutions/last, 0 ... disable)", NULL, 0 );
}

void output_initialize()
{
    get_parameter( "Output/i_output_data", ParameterDigit, &i_output_data );
    do_output_data = (i_output_data != 0);

    const_string_t suffix = ".h5";
    output_file = allocate( sizeof( string_t ) * (strlen( title ) + strlen( suffix ) + 1) );
    strcpy( title, output_file );
    strcat( output_file, suffix );
}

void output_finalize()
{
    deallocate( output_file );
}

void create_file_header()
{
        // file_id = create_hdf5_file( output_file )

        //     call set_hdf5_attribute( file_id, 'n_variables', n_variables )
        //     call set_hdf5_attribute( file_id, 'n_tot_variables', n_tot_variables )

        //     do i = 1, n_tot_variables; symbols(i) = variables(i)%name; end do
        //     call set_hdf5_dataset( file_id, 'variables', symbols )

        //     group_id = create_hdf5_group( file_id, 'SOLUTIONS' )
        //     call close_hdf5_group( group_id, 'SOLUTIONS' )

        // call close_hdf5_file( file_id, output_file )

        // if( get_is_root() ) call write_base_xdmf()
}

void write_output( int iter, double time )
{
        // if( .not. valid_hdf5_file( output_file ) ) call create_file_header()

        // tmp = set_string( iter, '(I0.9)' )
        // call write_output_solution( iter, time, tmp )
}

void write_output_solution( int iter, double time, string_t iter_string )
{
        // file_id = open_hdf5_file( output_file )
        // group_id = open_hdf5_group( file_id, 'SOLUTIONS' )

        // solution_id = create_hdf5_group( group_id, iter_string )
        //     call set_hdf5_attribute( solution_id, 'iter', iter )
        //     call set_hdf5_attribute( solution_id, 'time', time )

        //     if( get_is_parallel() ) then
        //         call set_hdf5_dataset( solution_id, 'phi', phi(:,:n_cells), &
        //             glob_rank=2, glob_dim=[n_variables,n_global_cells], stride=partition_cells )
        //         call set_hdf5_dataset( solution_id, 'phi_total', phi_total(:,:n_cells), &
        //             glob_rank=2, glob_dim=[n_tot_variables,n_global_cells], stride=partition_cells )
        //         call set_hdf5_dataset( solution_id, 'phi_dt', phi_dt(:,:n_cells), &
        //             glob_rank=2, glob_dim=[n_variables,n_global_cells], stride=partition_cells )

        //         if( n_stages .gt. 0 ) then
        //             call set_hdf5_attribute( solution_id, 'n_stages', n_stages )
        //             do i_stage = 1, n_stages
        //                 call set_hdf5_dataset( solution_id, 'phi_old:' // set_string( i_stage ), phi_old(:,:,i_stage), &
        //                     glob_rank=2, glob_dim=[n_variables,n_global_cells], stride=partition_cells )
        //             end do
        //         end if
        //     else
        //         call set_hdf5_dataset( solution_id, 'phi', phi(:,:n_cells) )
        //         call set_hdf5_dataset( solution_id, 'phi_total', phi_total(:,:n_cells) )
        //         call set_hdf5_dataset( solution_id, 'phi_dt', phi_dt(:,:n_cells) )

        //         if( n_stages .gt. 0 ) then
        //             call set_hdf5_attribute( solution_id, 'n_stages', n_stages )
        //             do i_stage = 1, n_stages
        //                 call set_hdf5_dataset( solution_id, 'phi_old:' // set_string( i_stage ), phi_old(:,:,i_stage) )
        //             end do
        //         end if
        //     end if

        // call close_hdf5_group( solution_id, iter_string )
        // call close_hdf5_group( group_id, 'SOLUTIONS' )

        // if( exists_hdf5_link( file_id, 'SOLUTION' ) ) call delete_hdf5_link( file_id, 'SOLUTION' )
        // call set_hdf5_soft_link( file_id, 'SOLUTIONS/' // iter_string, 'SOLUTION' )

        // call close_hdf5_file( file_id, output_file )
}