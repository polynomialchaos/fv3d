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
Mesh_t *global_mesh = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_initialize();
void mesh_finalize();
void read_mesh_file( Mesh_t *mesh );
void remap_local_mesh( Mesh_t *mesh );
void calc_mesh_metrics( Mesh_t *mesh );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_define()
{
    register_initialize_routine( mesh_initialize );
    register_finalize_routine( mesh_finalize );

    string_t tmp = "untitled.mesh->h5";
    set_parameter( "Mesh/mesh_file", ParameterString, &tmp, "The mesh file", NULL, 0 );
}

void mesh_initialize()
{
    get_parameter( "Mesh/mesh_file", ParameterString, &mesh_file );

    global_mesh = allocate_mesh();

    read_mesh_file( global_mesh );
    if (get_is_parallel()) remap_local_mesh( global_mesh );
    calc_mesh_metrics( global_mesh );
}

void mesh_finalize()
{
    deallocate( mesh_file );
    deallocate_mesh( &global_mesh );
}

void read_mesh_file( Mesh_t *mesh )
{
    int is_parallel = get_is_parallel();
    int rank = get_rank();

    hid_t file_id = open_hdf5_file( mesh_file );
        get_hdf5_attribute( file_id, "dimension", HDF5Int, &mesh->dimension );
        get_hdf5_attribute( file_id, "is_partitioned", HDF5Int, &mesh->is_partitioned );

        // the partition
        if (is_parallel)
        {
            check_error( (mesh->is_partitioned == 1) );

            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "PARTITIONS" );
                get_hdf5_attribute( group_id, "n_partitions", HDF5Int, &mesh->n_partitions );
                check_error( (mesh->n_partitions == get_n_procs()) );

                dims[0] = 1; offset[0] = rank; count[0] = 1;
                get_hdf5_dataset_chunk_n( group_id, "n_partition_cells", HDF5Int, &mesh->n_partition_cells,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_boundaries", HDF5Int, &mesh->n_partition_boundaries,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_faces", HDF5Int, &mesh->n_partition_faces,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_sends", HDF5Int, &mesh->n_partition_sends,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_receives", HDF5Int, &mesh->n_partition_receives,
                    1, dims, NULL, offset, count, NULL, NULL );

                allocate_partition( mesh->partition, mesh->n_partitions, mesh->n_partition_cells, mesh->n_partition_boundaries,
                    mesh->n_partition_faces, mesh->n_partition_sends, mesh->n_partition_receives );

                dims[0] = mesh->n_partition_cells; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = mesh->n_partition_cells;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_cells", HDF5Int, mesh->partition->partition_cells,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = mesh->n_partition_boundaries; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = mesh->n_partition_boundaries;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_boundaries", HDF5Int, mesh->partition->partition_boundaries,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = mesh->n_partition_faces; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = mesh->n_partition_faces;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_faces", HDF5Int, mesh->partition->partition_faces,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = mesh->n_partition_sends; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = mesh->n_partition_sends;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends", HDF5Int, mesh->partition->partition_sends,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends_pid", HDF5Int, mesh->partition->partition_sends_pid,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = mesh->n_partition_receives; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = mesh->n_partition_receives;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives", HDF5Int, mesh->partition->partition_receives,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives_pid", HDF5Int, mesh->partition->partition_receives_pid,
                    1, dims, NULL, offset, count, NULL, NULL );
            close_hdf5_group( group_id );
        }

        // the vertices
        {
            hsize_t dims[2];

            hid_t group_id = open_hdf5_group( file_id, "VERTICES" );
                get_hdf5_attribute( group_id, "n_vertices", HDF5Int, &mesh->n_vertices );
                allocate_vertices( mesh->vertices, mesh->n_vertices );

                dims[0] = 3;
                dims[1] = mesh->n_vertices;
                get_hdf5_dataset_chunk_n_m( group_id, "x", HDF5Int, mesh->vertices->x,
                    2, dims, NULL, NULL, NULL, NULL, NULL );
            close_hdf5_group( group_id );
        }

        // the cells
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "CELLS" );
                get_hdf5_attribute( group_id, "n_cells", HDF5Int, &mesh->n_cells );
                get_hdf5_attribute( group_id, "max_cell_vertices", HDF5Int, &mesh->max_cell_vertices );
                get_hdf5_attribute( group_id, "max_cell_faces", HDF5Int, &mesh->max_cell_faces );

                mesh->n_global_cells = mesh->n_cells;
                if (is_parallel) mesh->n_cells = mesh->n_partition_cells;
                mesh->n_local_cells = mesh->n_cells + mesh->n_partition_receives;

                allocate_cells( mesh->cells, mesh->n_local_cells, mesh->max_cell_vertices, mesh->max_cell_faces );

                hsize_t *stride = allocate( sizeof( hsize_t ) * mesh->n_local_cells );
                if (is_parallel)
                {
                    for ( int i = 0; i < mesh->n_partition_cells; i++)
                        stride[i] = mesh->partition->partition_cells[i];
                    for ( int i = mesh->n_partition_cells; i < mesh->n_local_cells; i++)
                        stride[i] = mesh->partition->partition_receives[i-mesh->n_partition_cells];
                }
                else
                {
                    for ( int i = 0; i < mesh->n_local_cells; i++)
                        stride[i] = i;
                }

                dims[0] = mesh->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "id", HDF5Int, mesh->cells->id,
                    1, dims, NULL, stride, mesh->n_local_cells );

                dims[0] = mesh->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, mesh->cells->type,
                    1, dims, NULL, stride, mesh->n_local_cells );

                dims[0] = mesh->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, mesh->cells->n_vertices,
                    1, dims, NULL, stride, mesh->n_local_cells );

                dims[0] = mesh->n_local_cells; offset[0] = 0; count[0] = mesh->n_local_cells;
                dims[1] = mesh->max_cell_vertices; offset[1] = 0; count[1] = mesh->max_cell_vertices;
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, mesh->cells->vertices,
                    2, dims, NULL, offset, count, stride, mesh->n_local_cells );

                dims[0] = mesh->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "n_faces", HDF5Int, mesh->cells->n_faces,
                    1, dims, NULL, stride, mesh->n_local_cells );

                dims[0] = mesh->n_local_cells; offset[0] = 0; count[0] = mesh->n_local_cells;
                dims[1] = mesh->max_cell_faces; offset[1] = 0; count[1] = mesh->max_cell_faces;
                get_hdf5_dataset_select_n_m( group_id, "faces", HDF5Int, mesh->cells->faces,
                    2, dims, NULL, offset, count, stride, mesh->n_local_cells );

                dims[0] = mesh->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "volume", HDF5Double, mesh->cells->volume,
                    1, dims, NULL, stride, mesh->n_local_cells );

                dims[0] = mesh->n_local_cells; offset[0] = 0; count[0] = mesh->n_local_cells;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "x", HDF5Double, mesh->cells->x,
                    2, dims, NULL, offset, count, stride, mesh->n_local_cells );

                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the boundaries
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "BOUNDARIES" );
                get_hdf5_attribute( group_id, "n_boundaries", HDF5Int, &mesh->n_boundaries );
                get_hdf5_attribute( group_id, "max_boundary_vertices", HDF5Int, &mesh->max_boundary_vertices );

                mesh->n_global_boundaries = mesh->n_boundaries;
                if (is_parallel) mesh->n_boundaries = mesh->n_partition_boundaries;
                mesh->n_local_boundaries = mesh->n_boundaries;

                allocate_boundaries( mesh->boundaries, mesh->n_local_boundaries, mesh->max_boundary_vertices );

                hsize_t *stride = allocate( sizeof( hsize_t ) * mesh->n_local_boundaries );
                if (is_parallel)
                {
                    for ( int i = 0; i < mesh->n_partition_boundaries; i++)
                        stride[i] = mesh->partition->partition_boundaries[i];
                }
                else
                {
                    for ( int i = 0; i < mesh->n_local_boundaries; i++)
                        stride[i] = i;
                }

                dims[0] = mesh->n_local_boundaries;
                get_hdf5_dataset_select_n( group_id, "id", HDF5Int, mesh->boundaries->id,
                    1, dims, NULL, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries;
                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, mesh->boundaries->type,
                    1, dims, NULL, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries;
                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, mesh->boundaries->n_vertices,
                    1, dims, NULL, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries; offset[0] = 0; count[0] = mesh->n_local_boundaries;
                dims[1] = mesh->max_boundary_vertices; offset[1] = 0; count[1] = mesh->max_boundary_vertices;
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, mesh->boundaries->vertices,
                    2, dims, NULL, offset, count, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries;
                get_hdf5_dataset_select_n( group_id, "face", HDF5Int, mesh->boundaries->face,
                    1, dims, NULL, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries;
                get_hdf5_dataset_select_n( group_id, "distance", HDF5Double, mesh->boundaries->distance,
                    1, dims, NULL, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries; offset[0] = 0; count[0] = mesh->n_local_boundaries;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "n", HDF5Double, mesh->boundaries->n,
                    2, dims, NULL, offset, count, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries; offset[0] = 0; count[0] = mesh->n_local_boundaries;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t1", HDF5Double, mesh->boundaries->t1,
                    2, dims, NULL, offset, count, stride, mesh->n_local_boundaries );

                dims[0] = mesh->n_local_boundaries; offset[0] = 0; count[0] = mesh->n_local_boundaries;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t2", HDF5Double, mesh->boundaries->t2,
                    2, dims, NULL, offset, count, stride, mesh->n_local_boundaries );

                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the faces
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "FACES" );
                get_hdf5_attribute( group_id, "n_faces", HDF5Int, &mesh->n_faces );
                get_hdf5_attribute( group_id, "max_face_vertices", HDF5Int, &mesh->max_face_vertices );

                mesh->n_global_faces = mesh->n_faces;
                if (is_parallel) mesh->n_faces = mesh->n_partition_faces;
                mesh->n_local_faces = mesh->n_faces;

                allocate_faces( mesh->faces, mesh->n_faces, mesh->max_face_vertices );

                hsize_t *stride = allocate( sizeof( hsize_t ) * mesh->n_local_faces );
                if (is_parallel)
                {
                    for ( int i = 0; i < mesh->n_partition_faces; i++)
                        stride[i] = mesh->partition->partition_faces[i];
                }
                else
                {
                    for ( int i = 0; i < mesh->n_local_faces; i++)
                        stride[i] = i;
                }

                dims[0] = mesh->n_local_faces;
                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, mesh->faces->type,
                    1, dims, NULL, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces;
                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, mesh->faces->n_vertices,
                    1, dims, NULL, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces; offset[0] = 0; count[0] = mesh->n_local_faces;
                dims[1] = mesh->max_face_vertices; offset[1] = 0; count[1] = mesh->max_face_vertices;
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, mesh->faces->vertices,
                    2, dims, NULL, offset, count, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces; offset[0] = 0; count[0] = mesh->n_local_faces;
                dims[1] = 2; offset[1] = 0; count[1] = 2;
                get_hdf5_dataset_select_n_m( group_id, "cells", HDF5Int, mesh->faces->cells,
                    2, dims, NULL, offset, count, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces;
                get_hdf5_dataset_select_n( group_id, "boundary", HDF5Int, mesh->faces->boundary,
                    1, dims, NULL, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces;
                get_hdf5_dataset_select_n( group_id, "area", HDF5Double, mesh->faces->area,
                    1, dims, NULL, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces;
                get_hdf5_dataset_select_n( group_id, "lambda", HDF5Double, mesh->faces->lambda,
                    1, dims, NULL, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces; offset[0] = 0; count[0] = mesh->n_local_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "x", HDF5Double, mesh->faces->x,
                    2, dims, NULL, offset, count, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces; offset[0] = 0; count[0] = mesh->n_local_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "n", HDF5Double, mesh->faces->n,
                    2, dims, NULL, offset, count, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces; offset[0] = 0; count[0] = mesh->n_local_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t1", HDF5Double, mesh->faces->t1,
                    2, dims, NULL, offset, count, stride, mesh->n_local_faces );

                dims[0] = mesh->n_local_faces; offset[0] = 0; count[0] = mesh->n_local_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t2", HDF5Double, mesh->faces->t2,
                    2, dims, NULL, offset, count, stride, mesh->n_local_faces );

                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the regions
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            hid_t group_id = open_hdf5_group( file_id, "REGIONS" );
                get_hdf5_attribute( group_id, "n_regions", HDF5Int, &mesh->n_regions );

                int max_length = get_hdf5_dataset_size( group_id, "name" );
                allocate_regions( mesh->regions, mesh->n_regions, max_length );

                dims[0] = mesh->n_regions;
                get_hdf5_dataset_chunk_n( group_id, "name", HDF5String, mesh->regions->name[0],
                    1, dims, NULL, NULL, NULL, NULL, NULL );
            close_hdf5_group( group_id );
        }
    close_hdf5_file( file_id );
}

void remap_local_mesh( Mesh_t *mesh )
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

void calc_mesh_metrics( Mesh_t *mesh )
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