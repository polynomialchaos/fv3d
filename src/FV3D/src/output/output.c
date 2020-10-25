//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "main_private.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"
#include "output_module.h"
#include "fv/fv_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int i_output_data   = -1;
int do_output_data  = 0;

string_t output_file = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void output_initialize();
void output_finalize();

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

    output_file = allocate_strcat( title, ".h5" );
}

void output_finalize()
{
    deallocate( output_file );
}

void create_file_header()
{
    hid_t file_id = create_hdf5_file( output_file );

        set_hdf5_attribute( file_id, "n_variables", HDF5Int, &all_variables->n_tot_variables );

        size_t max_len = 0;
        for ( int i = 0; i < all_variables->n_tot_variables; i++ )
            max_len = u_max( max_len, strlen( all_variables->tot_variables[i]->name ) );

        string_t *tmp = allocate_hdf5_string_buffer( all_variables->n_tot_variables, max_len+1 );
        for ( int i = 0; i < all_variables->n_tot_variables; i++ )
            strcpy( tmp[i], all_variables->tot_variables[i]->name );

        hsize_t dims[2] = {all_variables->n_tot_variables, max_len+1};
            set_hdf5_dataset_n( file_id, "variables", HDF5String, tmp[0], dims );

        deallocate_hdf5_string_buffer( &tmp );

        hid_t group_id = create_hdf5_group( file_id, "SOLUTIONS" );
        close_hdf5_group( group_id );

    close_hdf5_file( file_id );
}

void write_output( int iter, double time )
{
    if (is_valid_hdf5_file( output_file ) == 0) create_file_header();
    int n_domain_cells  = global_mesh->cells->n_domain_cells;
    int n_tot_variables = all_variables->n_tot_variables;

    char iter_string[10];
    sprintf( iter_string, "%09d", iter );

    hid_t file_id = open_hdf5_file( output_file );

        hid_t group_id = open_hdf5_group( file_id, "SOLUTIONS" );

            hid_t solution_id = create_hdf5_group( group_id, iter_string );

                set_hdf5_attribute( solution_id, "iter", HDF5Int, &iter );
                set_hdf5_attribute( solution_id, "time", HDF5Double, &time );

                if (get_is_parallel())
                {
                    //  call set_hdf5_dataset( solution_id, 'phi', phi(:,:n_cells), &
                    //      glob_rank=2, glob_dim=[n_variables,n_global_cells], stride=partition_cells )
                    //  call set_hdf5_dataset( solution_id, 'phi_total', phi_total(:,:n_cells), &
                    //      glob_rank=2, glob_dim=[n_tot_variables,n_global_cells], stride=partition_cells )
                    //  call set_hdf5_dataset( solution_id, 'phi_dt', phi_dt(:,:n_cells), &
                    //      glob_rank=2, glob_dim=[n_variables,n_global_cells], stride=partition_cells )

                    //  if( n_stages .gt. 0 ) then
                    //      call set_hdf5_attribute( solution_id, 'n_stages', n_stages )
                    //      do i_stage = 1, n_stages
                    //          call set_hdf5_dataset( solution_id, 'phi_old:' // set_string( i_stage ), phi_old(:,:,i_stage), &
                    //              glob_rank=2, glob_dim=[n_variables,n_global_cells], stride=partition_cells )
                    //      end do
                    //  end if
                }
                else
                {
                    hsize_t dims[2] = {n_domain_cells, n_tot_variables};
                    set_hdf5_dataset_n_m( solution_id, "phi_total", HDF5Double, phi_total, 2, dims );

                    // call set_hdf5_dataset( solution_id, 'phi', phi(:,:n_cells) )
                    // call set_hdf5_dataset( solution_id, 'phi_total', phi_total(:,:n_cells) )
                    // call set_hdf5_dataset( solution_id, 'phi_dt', phi_dt(:,:n_cells) )
                    // if( n_stages .gt. 0 ) then
                    //  call set_hdf5_attribute( solution_id, 'n_stages', n_stages )
                    //      do i_stage = 1, n_stages
                    //          call set_hdf5_dataset( solution_id, 'phi_old:' // set_string( i_stage ), phi_old(:,:,i_stage) )
                    //      end do
                    //  end if
                }

            close_hdf5_group( solution_id );

        close_hdf5_group( group_id );

        if (exists_hdf5_link( file_id, "SOLUTION" ) ) delete_hdf5_link( file_id, "SOLUTION" );
        string_t tmp = allocate_strcat( "SOLUTIONS/", iter_string );
        create_hdf5_soft_link( file_id, "SOLUTION", tmp );
        deallocate( tmp );

    close_hdf5_file( file_id );
}