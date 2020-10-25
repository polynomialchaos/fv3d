//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "mesh_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
Mesh_t *global_mesh = NULL;

string_t mesh_file = NULL;

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

    string_t tmp = "untitled.mesh.h5";
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
    deallocate_mesh( &global_mesh );

    deallocate( mesh_file );
}

void read_mesh_file( Mesh_t *mesh )
{
    Partition_t     *partition  = mesh->partition;
    Vertices_t      *vertices   = mesh->vertices;
    Cells_t         *cells      = mesh->cells;
    Boundaries_t    *boundaries = mesh->boundaries;
    Faces_t         *faces      = mesh->faces;
    Regions_t       *regions    = mesh->regions;

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

            int n_partitions            = 0;
            int n_partition_cells       = 0;
            int n_partition_boundaries  = 0;
            int n_partition_faces       = 0;
            int n_partition_sends       = 0;
            int n_partition_receives    = 0;

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

                partition = allocate_partition( mesh, n_partitions, n_partition_cells, n_partition_boundaries, n_partition_faces,
                    n_partition_sends, n_partition_receives );

                dims[0] = partition->n_partition_cells; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = partition->n_partition_cells;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_cells", HDF5Int, partition->partition_cells,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = partition->n_partition_boundaries; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = partition->n_partition_boundaries;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_boundaries", HDF5Int, partition->partition_boundaries,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = partition->n_partition_faces; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = partition->n_partition_faces;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_faces", HDF5Int, partition->partition_faces,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = partition->n_partition_sends; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = partition->n_partition_sends;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends", HDF5Int, partition->partition_sends,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends_pid", HDF5Int, partition->partition_sends_pid,
                    1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = partition->n_partition_receives; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = partition->n_partition_receives;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives", HDF5Int, partition->partition_receives,
                    1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives_pid", HDF5Int, partition->partition_receives_pid,
                    1, dims, NULL, offset, count, NULL, NULL );
            close_hdf5_group( group_id );
        }

        // the vertices
        {
            hsize_t dims[2];

            int n_vertices = 0;

            hid_t group_id = open_hdf5_group( file_id, "VERTICES" );
                get_hdf5_attribute( group_id, "n_vertices", HDF5Int, &n_vertices );

                vertices = allocate_vertices( mesh, n_vertices );

                dims[0] = vertices->n_vertices;
                dims[1] = 3;
                get_hdf5_dataset_n_m( group_id, "x", HDF5Double, vertices->x, 2, dims );
            close_hdf5_group( group_id );
        }

        // the cells
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            int n_global_cells      = 0;
            int max_cell_vertices   = 0;
            int max_cell_faces      = 0;

            hid_t group_id = open_hdf5_group( file_id, "CELLS" );
                get_hdf5_attribute( group_id, "n_cells", HDF5Int, &n_global_cells );
                get_hdf5_attribute( group_id, "max_cell_vertices", HDF5Int, &max_cell_vertices );
                get_hdf5_attribute( group_id, "max_cell_faces", HDF5Int, &max_cell_faces );

                cells = allocate_cells( mesh, n_global_cells, max_cell_vertices, max_cell_faces );

                hsize_t *stride = allocate( sizeof( hsize_t ) * cells->n_local_cells );
                if (is_parallel)
                {
                    for ( int i = 0; i < partition->n_partition_cells; i++)
                        stride[i] = partition->partition_cells[i];
                    for ( int i = partition->n_partition_cells; i < cells->n_local_cells; i++)
                        stride[i] = partition->partition_receives[i-partition->n_partition_cells];
                }
                else
                {
                    for ( int i = 0; i < cells->n_local_cells; i++)
                        stride[i] = i;
                }

                dims[0] = cells->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "id", HDF5Int, cells->id,
                    1, dims, NULL, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, cells->type,
                    1, dims, NULL, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, cells->n_vertices,
                    1, dims, NULL, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells; offset[0] = 0; count[0] = cells->n_local_cells;
                dims[1] = cells->max_cell_vertices; offset[1] = 0; count[1] = cells->max_cell_vertices;
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, cells->vertices,
                    2, dims, NULL, offset, count, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "n_faces", HDF5Int, cells->n_faces,
                    1, dims, NULL, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells; offset[0] = 0; count[0] = cells->n_local_cells;
                dims[1] = cells->max_cell_faces; offset[1] = 0; count[1] = cells->max_cell_faces;
                get_hdf5_dataset_select_n_m( group_id, "faces", HDF5Int, cells->faces,
                    2, dims, NULL, offset, count, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells;
                get_hdf5_dataset_select_n( group_id, "volume", HDF5Double, cells->volume,
                    1, dims, NULL, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells; offset[0] = 0; count[0] = cells->n_local_cells;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "x", HDF5Double, cells->x,
                    2, dims, NULL, offset, count, stride, cells->n_local_cells );

                dims[0] = cells->n_local_cells; offset[0] = 0; count[0] = cells->n_local_cells;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "dx", HDF5Double, cells->dx,
                    2, dims, NULL, offset, count, stride, cells->n_local_cells );

                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the boundaries
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            int n_global_boundaries     = 0;
            int max_boundary_vertices   = 0;

            hid_t group_id = open_hdf5_group( file_id, "BOUNDARIES" );
                get_hdf5_attribute( group_id, "n_boundaries", HDF5Int, &n_global_boundaries );
                get_hdf5_attribute( group_id, "max_boundary_vertices", HDF5Int, &max_boundary_vertices );

                boundaries = allocate_boundaries( mesh, n_global_boundaries, max_boundary_vertices );

                hsize_t *stride = allocate( sizeof( hsize_t ) * boundaries->n_boundaries );
                if (is_parallel)
                {
                    for ( int i = 0; i < partition->n_partition_boundaries; i++)
                        stride[i] = partition->partition_boundaries[i];
                }
                else
                {
                    for ( int i = 0; i < boundaries->n_boundaries; i++)
                        stride[i] = i;
                }

                dims[0] = boundaries->n_boundaries;
                get_hdf5_dataset_select_n( group_id, "id", HDF5Int, boundaries->id,
                    1, dims, NULL, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries;
                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, boundaries->type,
                    1, dims, NULL, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries;
                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, boundaries->n_vertices,
                    1, dims, NULL, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries; offset[0] = 0; count[0] = boundaries->n_boundaries;
                dims[1] = boundaries->max_boundary_vertices; offset[1] = 0; count[1] = boundaries->max_boundary_vertices;
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, boundaries->vertices,
                    2, dims, NULL, offset, count, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries;
                get_hdf5_dataset_select_n( group_id, "face", HDF5Int, boundaries->face,
                    1, dims, NULL, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries;
                get_hdf5_dataset_select_n( group_id, "distance", HDF5Double, boundaries->distance,
                    1, dims, NULL, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries; offset[0] = 0; count[0] = boundaries->n_boundaries;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "n", HDF5Double, boundaries->n,
                    2, dims, NULL, offset, count, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries; offset[0] = 0; count[0] = boundaries->n_boundaries;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t1", HDF5Double, boundaries->t1,
                    2, dims, NULL, offset, count, stride, boundaries->n_boundaries );

                dims[0] = boundaries->n_boundaries; offset[0] = 0; count[0] = boundaries->n_boundaries;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t2", HDF5Double, boundaries->t2,
                    2, dims, NULL, offset, count, stride, boundaries->n_boundaries );

                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the faces
        {
            hsize_t dims[2];
            hsize_t offset[2];
            hsize_t count[2];

            int n_global_faces      = 0;
            int max_face_vertices   = 0;

            hid_t group_id = open_hdf5_group( file_id, "FACES" );
                get_hdf5_attribute( group_id, "n_faces", HDF5Int, &n_global_faces );
                get_hdf5_attribute( group_id, "max_face_vertices", HDF5Int, &max_face_vertices );

                faces = allocate_faces( mesh, n_global_faces, max_face_vertices );

                hsize_t *stride = allocate( sizeof( hsize_t ) * faces->n_faces );
                if (is_parallel)
                {
                    for ( int i = 0; i < partition->n_partition_faces; i++)
                        stride[i] = partition->partition_faces[i];
                }
                else
                {
                    for ( int i = 0; i < faces->n_faces; i++)
                        stride[i] = i;
                }

                dims[0] = faces->n_faces;
                get_hdf5_dataset_select_n( group_id, "type", HDF5Int, faces->type,
                    1, dims, NULL, stride, faces->n_faces );

                dims[0] = faces->n_faces;
                get_hdf5_dataset_select_n( group_id, "n_vertices", HDF5Int, faces->n_vertices,
                    1, dims, NULL, stride, faces->n_faces );

                dims[0] = faces->n_faces; offset[0] = 0; count[0] = faces->n_faces;
                dims[1] = faces->max_face_vertices; offset[1] = 0; count[1] = faces->max_face_vertices;
                get_hdf5_dataset_select_n_m( group_id, "vertices", HDF5Int, faces->vertices,
                    2, dims, NULL, offset, count, stride, faces->n_faces );

                dims[0] = faces->n_faces; offset[0] = 0; count[0] = faces->n_faces;
                dims[1] = 2; offset[1] = 0; count[1] = 2;
                get_hdf5_dataset_select_n_m( group_id, "cells", HDF5Int, faces->cells,
                    2, dims, NULL, offset, count, stride, faces->n_faces );

                dims[0] = faces->n_faces;
                get_hdf5_dataset_select_n( group_id, "boundary", HDF5Int, faces->boundary,
                    1, dims, NULL, stride, faces->n_faces );

                dims[0] = faces->n_faces;
                get_hdf5_dataset_select_n( group_id, "area", HDF5Double, faces->area,
                    1, dims, NULL, stride, faces->n_faces );

                dims[0] = faces->n_faces;
                get_hdf5_dataset_select_n( group_id, "lambda", HDF5Double, faces->lambda,
                    1, dims, NULL, stride, faces->n_faces );

                dims[0] = faces->n_faces; offset[0] = 0; count[0] = faces->n_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "x", HDF5Double, faces->x,
                    2, dims, NULL, offset, count, stride, faces->n_faces );

                dims[0] = faces->n_faces; offset[0] = 0; count[0] = faces->n_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "n", HDF5Double, faces->n,
                    2, dims, NULL, offset, count, stride, faces->n_faces );

                dims[0] = faces->n_faces; offset[0] = 0; count[0] = faces->n_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t1", HDF5Double, faces->t1,
                    2, dims, NULL, offset, count, stride, faces->n_faces );

                dims[0] = faces->n_faces; offset[0] = 0; count[0] = faces->n_faces;
                dims[1] = 3; offset[1] = 0; count[1] = 3;
                get_hdf5_dataset_select_n_m( group_id, "t2", HDF5Double, faces->t2,
                    2, dims, NULL, offset, count, stride, faces->n_faces );

                deallocate( stride );
            close_hdf5_group( group_id );
        }

        // the regions
        {
            hsize_t dims[2];

            int n_regions = 0;

            hid_t group_id = open_hdf5_group( file_id, "REGIONS" );
                get_hdf5_attribute( group_id, "n_regions", HDF5Int, &n_regions );

                int max_length = get_hdf5_dataset_size( group_id, "name" );
                regions = allocate_regions( mesh, n_regions, max_length );

                dims[0] = regions->n_regions;
                get_hdf5_dataset_n( group_id, "name", HDF5String, regions->name[0], 1, dims );

                dims[0] = regions->n_regions;
                get_hdf5_dataset_n( group_id, "is_boundary", HDF5Int, regions->is_boundary, 1, dims );
            close_hdf5_group( group_id );
        }
    close_hdf5_file( file_id );
}

void remap_local_mesh( Mesh_t *mesh )
{
    Partition_t     *partition  = mesh->partition;
    Cells_t         *cells      = mesh->cells;
    Boundaries_t    *boundaries = mesh->boundaries;
    Faces_t         *faces      = mesh->faces;

    int *global_cells       = allocate( sizeof( int ) * cells->n_global_cells );
    int *global_boundaries  = allocate( sizeof( int ) * boundaries->n_global_boundaries );
    int *global_faces       = allocate( sizeof( int ) * faces->n_global_faces );

    set_value_int_n( -1, global_cells, cells->n_global_cells );
    set_value_int_n( -1, global_boundaries, boundaries->n_global_boundaries );
    set_value_int_n( -1, global_faces, faces->n_global_faces );

    for ( int i = 0; i < cells->n_domain_cells; i++ )
        global_cells[partition->partition_cells[i]] = i;

    for ( int i = 0; i < partition->n_partition_sends; i++ )
        check_error( (global_cells[partition->partition_sends[i]] != -1) );

    for ( int i = 0; i < partition->n_partition_receives; i++ )
        global_cells[partition->partition_receives[i]] = cells->n_domain_cells + i;

    for ( int i = 0; i < partition->n_partition_sends; i++ )
        partition->partition_sends[i] = global_cells[partition->partition_sends[i]];

    for ( int i = 0; i < partition->n_partition_receives; i++ )
        partition->partition_receives[i] = global_cells[partition->partition_receives[i]];

    for ( int i = 0; i < boundaries->n_boundaries; i++ )
        global_boundaries[partition->partition_boundaries[i]] = i;

    for ( int i = 0; i < faces->n_faces; i++ )
        global_faces[partition->partition_faces[i]] = i;

    for ( int i = 0; i < cells->n_local_cells; i++ )
    {
        int *cell_faces = &cells->faces[i*cells->max_cell_faces];

        for ( int j = 0; j < cells->n_faces[i]; j++ )
            cell_faces[j] = global_faces[cell_faces[j]];
    }

    for ( int i = 0; i < boundaries->n_boundaries; i++ )
    {
        int *boundary_face  = &boundaries->face[i];

        boundary_face[0] = global_faces[boundary_face[0]];
    }

    for ( int i = 0; i < faces->n_faces; i++ )
    {
        int *fc     = &faces->cells[i*FACE_CELLS];
        int *face_boundary  = &faces->boundary[i];

        fc[0] = global_cells[fc[0]];

        if (fc[1] >= 0)
        {
            fc[1] = global_cells[fc[1]];
        }
        else
        {
            face_boundary[0] = global_boundaries[face_boundary[0]];
        }

        if (fc[0] >= cells->n_domain_cells)
        {
            int n = cells->n_faces[fc[0]];
            int *cell_faces = &cells->faces[fc[0]*cells->max_cell_faces];

            cells->n_faces[fc[0]] = 0;
            for ( int j = 0; j < n; j++ )
            {
                if (cell_faces[j] < 0) continue;
                cell_faces[cells->n_faces[fc[0]]] = cell_faces[j];
                cells->n_faces[fc[0]] += 1;
            }
        }

        if (fc[1] >= cells->n_domain_cells)
        {
            int n = cells->n_faces[fc[1]];
            int *cell_faces = &cells->faces[fc[1]*cells->max_cell_faces];

            cells->n_faces[fc[1]] = 0;
            for ( int j = 0; j < n; j++ )
            {
                if (cell_faces[j] < 0) continue;
                cell_faces[cells->n_faces[fc[1]]] = cell_faces[j];
                cells->n_faces[fc[1]] += 1;
            }
        }
    }

    for ( int i = 0; i < partition->n_partitions; i++ )
    {
        if (i == get_rank()) continue;

        int *partition_sends_to = &partition->partition_sends_to[i*partition->n_partition_sends];
        for ( int j = 0; j < partition->n_partition_sends; j++ )
        {
            if (partition->partition_sends_pid[j] != i) continue;
            partition_sends_to[partition->n_partition_sends_to[i]] =
                partition->partition_sends[j];
            partition->n_partition_sends_to[i] += 1;
        }

        int *partition_receives_from = &partition->partition_receives_from[i*partition->n_partition_receives];
        for ( int j = 0; j < partition->n_partition_receives; j++ )
        {
            if (partition->partition_receives_pid[j] != i) continue;
            partition_receives_from[partition->n_partition_receives_from[i]] =
                partition->partition_receives[j];
            partition->n_partition_receives_from[i] += 1;
        }
    }

    deallocate( global_cells );
    deallocate( global_boundaries );
    deallocate( global_faces );
}

void calc_mesh_metrics( Mesh_t *mesh )
{
    Cells_t     *cells      = mesh->cells;
    Faces_t     *faces      = mesh->faces;
    Regions_t   *regions    = mesh->regions;

    mesh->local_volume = sum_n( &cells->volume[0], cells->n_domain_cells );
    mpi_all_reduce( &mesh->local_volume, &mesh->global_volume, MPIDouble, MPISum );

    faces->n_internal_faces = 0;
    for ( int i = 0; i < faces->n_faces; i++)
        if (faces->boundary[i] < 0)
            faces->n_internal_faces += 1;

    faces->n_boundary_faces = faces->n_faces - faces->n_internal_faces;

    faces->dist_cell_1    = allocate( sizeof( double ) * faces->n_faces * DIM );
    faces->dist_cell_2    = allocate( sizeof( double ) * faces->n_faces * DIM );
    faces->internal_faces = allocate( sizeof( double ) * faces->n_internal_faces );
    faces->boundary_faces = allocate( sizeof( double ) * faces->n_boundary_faces );

    int k = 0; int l = 0;
    for ( int i = 0; i < faces->n_faces; i++)
    {
        int *fc = &faces->cells[i*FACE_CELLS];

        for ( int j = 0; j < DIM; j++)
            faces->dist_cell_1[i*DIM+j] = faces->x[i*DIM+j] -
                cells->x[fc[0]*DIM+j];

        if (fc[1] >= 0)
        {
            for ( int j = 0; j < DIM; j++)
                faces->dist_cell_2[i*DIM+j] = faces->x[i*DIM+j] -
                    cells->x[fc[1]*DIM+j];

            faces->internal_faces[k] = i;
            k += 1;
        }
        else
        {
            for ( int j = 0; j < DIM; j++)
                faces->dist_cell_2[i*DIM+j] = 0.0;

            fc[1] = cells->n_local_cells + faces->boundary[i];

            faces->boundary_faces[l] = i;
            l += 1;
        }
    }

    for ( int i = 0; i < regions->n_regions; i++)
    {
        if (regions->is_boundary[i] == 0)
        {
            regions->flow_region = i;
            break;
        }
    }
}