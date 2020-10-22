//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "mesh_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
string_t mesh_file = NULL;

int dimension = 0;
int is_partitioned = 0;

int n_partitions = 0;
int n_partition_cells = 0;
int n_partition_boundaries = 0;
int n_partition_faces = 0;
int n_partition_sends = 0;
int n_partition_receives = 0;

int *partition_cells = NULL;
int *partition_boundaries = NULL;
int *partition_faces = NULL;
int *partition_sends = NULL;
int *partition_sends_pid = NULL;
int *partition_receives = NULL;
int *partition_receives_pid = NULL;

int *n_partition_sends_to = NULL;
int *partition_sends_to = NULL;
int *n_partition_receives_from = NULL;
int *partition_receives_from = NULL;

int n_vertices = 0;
double *vertices = NULL;

int n_global_cells = 0;
int n_local_cells = 0;
int n_cells = 0;
Cell_t *cells = NULL;

int n_global_boundaries = 0;
int n_local_boundaries = 0;
int n_boundaries = 0;
Boundary_t *boundaries = NULL;

int n_global_faces = 0;
int n_local_faces = 0;
int n_faces = 0;
Face_t *faces = NULL;

int n_regions = 0;
Region_t *regions = NULL;

double *dist_cell_1 = NULL;
double *dist_cell_2 = NULL;

double total_volume = 0.0;

int flow_region = 0;

int n_internal_faces = 0;
int *internal_faces = NULL;
int n_boundary_faces = 0;
int *boundary_faces = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_initialize();
void mesh_finalize();
void read_mesh_file();
void remap_local_mesh();
void calc_mesh_metrics();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_define()
{
    register_initialize_routine( mesh_initialize );
    register_finalize_routine( mesh_finalize );

    string_t tmp = "untitled.mesh.h5";
    set_parameter( "Mesh/mesh_file", ParameterString, &tmp, "The mesh file", NULL, 0 );
}

void mesh_initialize()
{
    get_parameter( "Mesh/mesh_file", ParameterString, &mesh_file );

    read_mesh_file();
    if (get_is_parallel()) remap_local_mesh();
    calc_mesh_metrics();
}

void mesh_finalize()
{
    deallocate( mesh_file );

    deallocate( partition_cells );
    deallocate( partition_boundaries );
    deallocate( partition_faces );
    deallocate( partition_sends );
    deallocate( partition_sends_pid );
    deallocate( partition_receives );
    deallocate( partition_receives_pid );

    deallocate( n_partition_sends_to );
    deallocate( partition_sends_to );
    deallocate( n_partition_receives_from );
    deallocate( partition_receives_from );

    deallocate( vertices );

    deallocate( cells );
    deallocate( boundaries );
    deallocate( faces );
    deallocate( regions );

    deallocate( dist_cell_1 );
    deallocate( dist_cell_2 );

    deallocate( internal_faces );
    deallocate( boundary_faces );
}

void read_mesh_file()
{
    int is_parallel = get_is_parallel();
    int rank = get_rank();

    hid_t file_id = open_hdf5_file( mesh_file );
        get_hdf5_attribute( file_id, "dimension", HDF5Int, &dimension );
        get_hdf5_attribute( file_id, "is_partitioned", HDF5Int, &is_partitioned );

        // the partition
        if (is_parallel)
        {
            check_error( (is_partitioned == 1) );

            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "PARTITIONS" );
                get_hdf5_attribute( group_id, "n_partitions", HDF5Int, &n_partitions );
                check_error( (n_partitions == get_n_procs()) );

                dims[0] = 1; offset[0] = rank; count[0] = 1;
                get_hdf5_dataset_chunk_n( group_id, "n_partition_cells", HDF5Int, &n_partition_cells,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_boundaries", HDF5Int, &n_partition_boundaries,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_faces", HDF5Int, &n_partition_faces,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_sends", HDF5Int, &n_partition_sends,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_receives", HDF5Int, &n_partition_receives,
                    1, dims, NULL, offset, count, NULL, NULL );

                partition_cells         = allocate( sizeof( int ) * n_partition_cells );
                partition_boundaries    = allocate( sizeof( int ) * n_partition_boundaries );
                partition_faces         = allocate( sizeof( int ) * n_partition_faces );
                partition_sends         = allocate( sizeof( int ) * n_partition_sends );
                partition_sends_pid     = allocate( sizeof( int ) * n_partition_sends );
                partition_receives      = allocate( sizeof( int ) * n_partition_receives );
                partition_receives_pid  = allocate( sizeof( int ) * n_partition_receives );

                dims[0] = n_partition_cells; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_cells;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_cells", HDF5Int, partition_cells,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_boundaries; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_boundaries;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_boundaries", HDF5Int, partition_boundaries,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_faces; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_faces;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_faces", HDF5Int, partition_faces,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_sends; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_sends;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends", HDF5Int, partition_sends,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends_pid", HDF5Int, partition_sends_pid,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_receives; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_receives;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives", HDF5Int, partition_receives,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives_pid", HDF5Int, partition_receives_pid,
                    1, dims, NULL, offset, count, NULL, NULL );
            close_hdf5_group( group_id );
        }

        // the vertices
        {
            hsize_t dims[2];

            hid_t group_id = open_hdf5_group( file_id, "VERTICES" );
                get_hdf5_attribute( group_id, "n_vertices", HDF5Int, &n_vertices );

                vertices = allocate( sizeof( double ) * 3 * n_vertices );
                dims[0] = 3;
                dims[1] = n_vertices;
                get_hdf5_dataset_chunk_n_m( group_id, "x", HDF5Int, vertices,
                    2, dims, NULL, NULL, NULL, NULL, NULL );
            close_hdf5_group( group_id );
        }


        // the cells
        {
            int max_cell_vertices = 0;
            int max_cell_faces = 0;

            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "CELLS" );
                get_hdf5_attribute( group_id, "n_cells", HDF5Int, &n_cells );
                get_hdf5_attribute( group_id, "max_cell_vertices", HDF5Int, &max_cell_vertices );
                get_hdf5_attribute( group_id, "max_cell_faces", HDF5Int, &max_cell_faces );

                n_global_cells = n_cells;
                if (is_parallel) n_cells = n_partition_cells;
                n_local_cells = n_cells + n_partition_receives;

                cells = allocate( sizeof( Cell_t ) * n_local_cells );

                hsize_t *stride = allocate( sizeof( hsize_t ) * n_local_cells );
                if (is_parallel)
                {
                    for ( int i = 0; i < n_partition_cells; i++)
                        stride[i] = partition_cells[i];
                    for ( int i = n_partition_cells; i < n_local_cells; i++)
                        stride[i] = partition_receives[i-n_partition_cells];
                }
                else
                {
                    for ( int i = 0; i < n_local_cells; i++)
                        stride[i] = i;
                }

                dims[0] = n_local_cells;
                int *tmp_i = allocate( sizeof(int) * n_local_cells );
                get_hdf5_dataset_select_n( group_id, "id", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                    cells[i].id = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                    cells[i].type = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                    cells[i].n_vertices = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "n_faces", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                    cells[i].n_faces = tmp_i[i];
                deallocate( tmp_i );

                dims[0] = n_local_cells; offset[0] = 0; count[0] = n_local_cells;
                dims[1] = max_cell_vertices; offset[1] = 0; count[1] = max_cell_vertices;
                tmp_i = allocate( sizeof(int) * max_cell_vertices * n_local_cells );
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, tmp_i,
                    2, dims, NULL, offset, count, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                {
                    cells[i].vertices = allocate( sizeof( int ) * cells[i].n_vertices );
                    for ( int j = 0; j < cells[i].n_vertices; j++)
                        cells[i].vertices[j] = tmp_i[i*max_cell_vertices+j];
                }
                deallocate( tmp_i );

                dims[0] = n_local_cells; offset[0] = 0; count[0] = n_local_cells;
                dims[1] = max_cell_faces; offset[1] = 0; count[1] = max_cell_faces;
                tmp_i = allocate( sizeof(int) * max_cell_faces * n_local_cells );
                get_hdf5_dataset_select_n_m( group_id, "faces", HDF5Int, tmp_i,
                    2, dims, NULL, offset, count, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                {
                    cells[i].faces = allocate( sizeof( int ) * cells[i].n_faces );
                    for ( int j = 0; j < cells[i].n_faces; j++)
                        cells[i].faces[j] = tmp_i[i*max_cell_faces+j];
                }
                deallocate( tmp_i );

                dims[0] = n_local_cells;
                double *tmp_d = allocate( sizeof( double ) * n_local_cells );
                get_hdf5_dataset_select_n( group_id, "volume", HDF5Double, tmp_d,
                    1, dims, NULL, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                    cells[i].volume = tmp_d[i];
                deallocate( tmp_d );

                dims[0] = n_local_cells; offset[0] = 0; count[0] = n_local_cells;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                tmp_d = allocate( sizeof( double ) * 3 * n_local_cells );
                get_hdf5_dataset_select_n_m( group_id, "x", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_cells );
                for ( int i = 0; i < n_local_cells; i++)
                    for ( int j = 0; j < 3; j++)
                        cells[i].x[j] = tmp_d[i*3+j];
                deallocate( tmp_d );
                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the boundaries
        {
            int max_boundary_vertices = 0;

            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "BOUNDARIES" );
                get_hdf5_attribute( group_id, "n_boundaries", HDF5Int, &n_boundaries );
                get_hdf5_attribute( group_id, "max_boundary_vertices", HDF5Int, &max_boundary_vertices );

                n_global_boundaries = n_boundaries;
                if (is_parallel) n_boundaries = n_partition_boundaries;
                n_local_boundaries = n_boundaries;

                boundaries = allocate( sizeof( Boundary_t ) * n_local_boundaries );

                hsize_t *stride = allocate( sizeof( hsize_t ) * n_local_boundaries );
                if (is_parallel)
                {
                    for ( int i = 0; i < n_partition_boundaries; i++)
                        stride[i] = partition_boundaries[i];
                }
                else
                {
                    for ( int i = 0; i < n_local_boundaries; i++)
                        stride[i] = i;
                }

                dims[0] = n_local_boundaries;
                int *tmp_i = allocate( sizeof(int) * n_local_boundaries );
                get_hdf5_dataset_select_n( group_id, "id", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    boundaries[i].id = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    boundaries[i].type = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    boundaries[i].n_vertices = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "face", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    boundaries[i].face = tmp_i[i];
                deallocate( tmp_i );

                dims[0] = n_local_boundaries; offset[0] = 0; count[0] = n_local_boundaries;
                dims[1] = max_boundary_vertices; offset[1] = 0; count[1] = max_boundary_vertices;
                tmp_i = allocate( sizeof(int) * max_boundary_vertices * n_local_boundaries );
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, tmp_i,
                    2, dims, NULL, offset, count, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                {
                    boundaries[i].vertices = allocate( sizeof( int ) * boundaries[i].n_vertices );
                    for ( int j = 0; j < boundaries[i].n_vertices; j++)
                        boundaries[i].vertices[j] = tmp_i[i*max_boundary_vertices+j];
                }
                deallocate( tmp_i );

                dims[0] = n_local_boundaries;
                double *tmp_d = allocate( sizeof( double ) * n_local_boundaries );
                get_hdf5_dataset_select_n( group_id, "distance", HDF5Double, tmp_d,
                    1, dims, NULL, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    boundaries[i].distance = tmp_d[i];
                deallocate( tmp_d );

                dims[0] = n_local_boundaries; offset[0] = 0; count[0] = n_local_boundaries;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                tmp_d = allocate( sizeof( double ) * 3 * n_local_boundaries );
                get_hdf5_dataset_select_n_m( group_id, "n", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    for ( int j = 0; j < 3; j++)
                        boundaries[i].n[j] = tmp_d[i*3+j];

                get_hdf5_dataset_select_n_m( group_id, "t1", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    for ( int j = 0; j < 3; j++)
                        boundaries[i].t1[j] = tmp_d[i*3+j];

                get_hdf5_dataset_select_n_m( group_id, "t2", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_boundaries );
                for ( int i = 0; i < n_local_boundaries; i++)
                    for ( int j = 0; j < 3; j++)
                        boundaries[i].t2[j] = tmp_d[i*3+j];
                deallocate( tmp_d );
                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the faces
        {
            int max_face_vertices = 0;

            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "FACES" );
                get_hdf5_attribute( group_id, "n_faces", HDF5Int, &n_faces );
                get_hdf5_attribute( group_id, "max_face_vertices", HDF5Int, &max_face_vertices );

                n_global_faces = n_faces;
                if (is_parallel) n_faces = n_partition_faces;
                n_local_faces = n_faces;

                faces = allocate( sizeof( Face_t ) * n_local_faces );

                hsize_t *stride = allocate( sizeof( hsize_t ) * n_local_faces );
                if (is_parallel)
                {
                    for ( int i = 0; i < n_partition_faces; i++)
                        stride[i] = partition_faces[i];
                }
                else
                {
                    for ( int i = 0; i < n_local_faces; i++)
                        stride[i] = i;
                }

                dims[0] = n_local_faces;
                int *tmp_i = allocate( sizeof(int) * n_local_faces );
                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    faces[i].type = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    faces[i].n_vertices = tmp_i[i];

                get_hdf5_dataset_select_n( group_id, "boundary", HDF5Int, tmp_i,
                    1, dims, NULL, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    faces[i].boundary = tmp_i[i];
                deallocate( tmp_i );

                dims[0] = n_local_faces; offset[0] = 0; count[0] = n_local_faces;
                dims[1] = max_face_vertices; offset[1] = 0; count[1] = max_face_vertices;
                tmp_i = allocate( sizeof(int) * max_face_vertices * n_local_faces );
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, tmp_i,
                    2, dims, NULL, offset, count, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                {
                    faces[i].vertices = allocate( sizeof( int ) * faces[i].n_vertices );
                    for ( int j = 0; j < faces[i].n_vertices; j++)
                        faces[i].vertices[j] = tmp_i[i*max_face_vertices+j];
                }
                deallocate( tmp_i );

                dims[0] = n_local_faces; offset[0] = 0; count[0] = n_local_faces;
                dims[1] = 2; offset[1] = 0; count[1] = 2;
                tmp_i = allocate( sizeof(int) * 2 * n_local_faces );
                get_hdf5_dataset_select_n_m( group_id, "cells", HDF5Int, tmp_i,
                    2, dims, NULL, offset, count, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                {
                    for ( int j = 0; j < 2; j++)
                        faces[i].cells[j] = tmp_i[i*2+j];
                }
                deallocate( tmp_i );

                dims[0] = n_local_faces;
                double *tmp_d = allocate( sizeof( double ) * n_local_faces );
                get_hdf5_dataset_select_n( group_id, "area", HDF5Double, tmp_d,
                    1, dims, NULL, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    faces[i].area = tmp_d[i];

                get_hdf5_dataset_select_n( group_id, "lambda", HDF5Double, tmp_d,
                    1, dims, NULL, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    faces[i].lambda = tmp_d[i];
                deallocate( tmp_d );

                dims[0] = n_local_faces; offset[0] = 0; count[0] = n_local_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                tmp_d = allocate( sizeof( double ) * 3 * n_local_faces );
                get_hdf5_dataset_select_n_m( group_id, "x", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    for ( int j = 0; j < 3; j++)
                        faces[i].x[j] = tmp_d[i*3+j];

                get_hdf5_dataset_select_n_m( group_id, "n", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    for ( int j = 0; j < 3; j++)
                        faces[i].n[j] = tmp_d[i*3+j];

                get_hdf5_dataset_select_n_m( group_id, "t1", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    for ( int j = 0; j < 3; j++)
                        faces[i].t1[j] = tmp_d[i*3+j];

                get_hdf5_dataset_select_n_m( group_id, "t2", HDF5Double, tmp_d,
                    2, dims, NULL, offset, count, stride, n_local_faces );
                for ( int i = 0; i < n_local_faces; i++)
                    for ( int j = 0; j < 3; j++)
                        faces[i].t2[j] = tmp_d[i*3+j];
                deallocate( tmp_d );
                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the regions
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "REGIONS" );
                get_hdf5_attribute( group_id, "n_regions", HDF5Int, &n_regions );

                regions = allocate( sizeof( Region_t ) * n_regions );

                dims[0] = n_regions;
                string_t *tmp_s = allocate_hdf5_string_buffer( n_regions, get_hdf5_dataset_size( group_id, "name" ) );
                get_hdf5_dataset_chunk_n( group_id, "name", HDF5String, tmp_s[0],
                    1, dims, NULL, NULL, NULL, NULL, NULL );
                for ( int i = 0; i < n_regions; i++ )
                {
                    regions[i].name = allocate( sizeof( char ) * (len_trim( tmp_s[i] ) + 1) );
                    strcpy( tmp_s[i], regions[i].name );
                }

                deallocate_hdf5_string_buffer( &tmp_s );
            close_hdf5_group( group_id );
        }
    close_hdf5_file( file_id );
}

void remap_local_mesh()
{
        // global_cells        = 0
        // global_boundaries   = 0
        // global_faces        = 0

        // ! define the global index arrays
        // do i = 1, n_cells
        //     global_cells(partition_cells(i)) = i
        // end do

        // do i = 1, n_partition_receives
        //     global_cells(partition_receives(i)) = n_cells + i
        // end do

        // do i = 1, n_partition_sends
        //     partition_sends(i) = global_cells(partition_sends(i))
        // end do

        // do i = 1, n_partition_receives
        //     partition_receives(i) = global_cells(partition_receives(i))
        // end do

        // do i = 1, n_boundaries
        //     global_boundaries(partition_boundaries(i)) = i
        // end do

        // do i = 1, n_faces
        //     global_faces(partition_faces(i)) = i
        // end do

        // ! set the local cell indices
        // do i = 1, n_local_cells
        //     do j = 1, cells(i)%n_faces
        //         cells(i)%faces(j) = global_faces(cells(i)%faces(j))
        //     end do
        // end do

        // ! set the local boundary indices
        // do i = 1, n_boundaries
        //     boundaries(i)%face = global_faces(boundaries(i)%face)
        // end do

        // ! set the local face cell indices
        // do i = 1, n_faces
        //     faces(i)%cells(1) = global_cells(faces(i)%cells(1))

        //     if( faces(i)%cells(2) .gt. 0 ) then
        //         faces(i)%cells(2) = global_cells(faces(i)%cells(2))
        //     else
        //         faces(i)%boundary = global_boundaries(faces(i)%boundary)
        //     end if
        // end do

        // ! set the sends to and receives from list
        // allocate( n_partition_sends_to(n_partitions) ); n_partition_sends_to = 0
        // allocate( partition_sends_to(n_partition_sends,n_partitions) ); partition_sends_to = 0
        // allocate( n_partition_receives_from(n_partitions) ); n_partition_receives_from = 0
        // allocate( partition_receives_from(n_partition_receives,n_partitions) ); partition_receives_from = 0

        // do i = 1, n_partitions
        //     if( i-1 .eq. get_i_rank() ) cycle

        //     do j = 1, n_partition_sends
        //         if( partition_sends_pid(j) .ne. i ) cycle
        //         n_partition_sends_to(i) = n_partition_sends_to(i) + 1
        //         partition_sends_to(n_partition_sends_to(i),i) = partition_sends(j)
        //     end do

        //     do j = 1, n_partition_receives
        //         if( partition_receives_pid(j) .ne. i ) cycle
        //         n_partition_receives_from(i) = n_partition_receives_from(i) + 1
        //         partition_receives_from(n_partition_receives_from(i),i) = partition_receives(j)
        //     end do
        // end do
}

void calc_mesh_metrics()
{
        // total_volume_loc = sum( cells(:)%volume )
        // call allreduce_mpi( total_volume_loc, total_volume, MPI_SUM )

        // allocate( dist_cell_1(3,n_faces) ); dist_cell_1 = 0.0
        // allocate( dist_cell_2(3,n_faces) ); dist_cell_2 = 0.0

        // n_internal_faces = count( faces(:)%boundary .eq. 0 )
        // n_boundary_faces = n_faces - n_internal_faces
        // allocate(internal_faces(n_internal_faces)); internal_faces = 0
        // allocate(boundary_faces(n_boundary_faces)); boundary_faces = 0

        // ! define the cell metrics
        // do i = 1, n_cells
        //     do j = 1, cells(i)%n_faces
        //         cells(i)%ds = cells(i)%ds + &
        //             abs( faces(cells(i)%faces(j))%n ) * faces(cells(i)%faces(j))%area
        //     end do
        // end do

        // ! communication faces need to be updated with valid faces
        // do i = 1, n_faces
        //     fc = faces(i)%cells

        //     if( fc(1) .gt. n_cells ) then
        //         cells(fc(1))%n_faces = count( cells(fc(1))%faces .gt. 0 )
        //         cells(fc(1))%faces(1:cells(fc(1))%n_faces-1) = &
        //             pack( cells(fc(1))%faces, cells(fc(1))%faces .gt. 0 )
        //     endif

        //     if( fc(2) .gt. n_cells ) then
        //         cells(fc(2))%n_faces = count( cells(fc(2))%faces .gt. 0 )
        //         cells(fc(2))%faces(1:cells(fc(1))%n_faces-1) = &
        //             pack( cells(fc(2))%faces, cells(fc(2))%faces .gt. 0 )
        //     endif
        // end do

        // iif = 0; ibf = 0
        // do i = 1, n_faces
        //     fc = faces(i)%cells

        //     dist_cell_1(:,i) = faces(i)%x - cells(fc(1))%x

        //     if( fc(2) .gt. 0 ) then
        //         dist_cell_2(:,i)    = faces(i)%x - cells(fc(2))%x
        //         iif                 = iif + 1
        //         internal_faces(iif) = i
        //     else
        //         dist_cell_2(:,i)    = 0.0
        //         faces(i)%cells(2)   = n_local_cells + faces(i)%boundary
        //         ibf                 = ibf + 1
        //         boundary_faces(ibf) = i
        //     end if
        // end do

        // do i = 1, n_boundaries
        //     regions(boundaries(i)%id)%is_boundary = .true.
        // end do

        // do i = 1, n_regions
        //     if( regions(i)%is_boundary ) then
        //         flow_region = i
        //         exit
        //     end if
        // end do
}