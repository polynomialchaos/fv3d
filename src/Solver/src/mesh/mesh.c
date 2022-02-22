/*******************************************************************************
 * @file mesh.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include <string.h>
#include "mesh_private.h"

Mesh_t *solver_mesh = NULL;

/*******************************************************************************
 * @brief Calculate the mesh metrics
 * @param mesh
 ******************************************************************************/
void calc_mesh_metrics(Mesh_t *mesh)
{
    Cells_t *cells = mesh->cells;
    Faces_t *faces = mesh->faces;
    Regions_t *regions = mesh->regions;

    mesh->local_volume = sum_n(&cells->volume[0], cells->n_domain_cells);
    MPI_ALL_REDUCE(MPIDouble, MPISum, &mesh->local_volume, &mesh->global_volume);

    faces->n_internal_faces = 0;
    for (int i = 0; i < faces->n_faces; ++i)
        if (faces->boundary[i] < 0)
            faces->n_internal_faces += 1;

    faces->n_boundary_faces = faces->n_faces - faces->n_internal_faces;

    faces->dist_cell_1 = ALLOCATE(sizeof(double) * faces->n_faces * DIM);
    faces->dist_cell_2 = ALLOCATE(sizeof(double) * faces->n_faces * DIM);
    faces->internal_faces = ALLOCATE(sizeof(int) * faces->n_internal_faces);
    faces->boundary_faces = ALLOCATE(sizeof(int) * faces->n_boundary_faces);

    int k = 0;
    int l = 0;
    for (int i = 0; i < faces->n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];

        for (int j = 0; j < DIM; ++j)
            faces->dist_cell_1[i * DIM + j] = faces->x[i * DIM + j] - cells->x[fc[0] * DIM + j];

        if (fc[1] >= 0)
        {
            for (int j = 0; j < DIM; ++j)
                faces->dist_cell_2[i * DIM + j] = faces->x[i * DIM + j] - cells->x[fc[1] * DIM + j];

            faces->internal_faces[k] = i;
            k += 1;
        }
        else
        {
            for (int j = 0; j < DIM; ++j)
                faces->dist_cell_2[i * DIM + j] = 0.0;

            fc[1] = cells->n_local_cells + faces->boundary[i];

            faces->boundary_faces[l] = i;
            l += 1;
        }
    }

    CHECK_EXPRESSION(k == faces->n_internal_faces);
    CHECK_EXPRESSION(l == faces->n_boundary_faces);

    for (int i = 0; i < regions->n_regions; ++i)
    {
        if (regions->is_boundary[i] == 0)
        {
            regions->flow_region = i;
            break;
        }
    }
}

/*******************************************************************************
 * @brief Free mesh
 ******************************************************************************/
void free_mesh()
{
    deallocate_mesh(solver_mesh);
    DEALLOCATE(solver_mesh);
}

/*******************************************************************************
 * @brief Initialize mesh
 * @param mesh_file
 ******************************************************************************/
void init_mesh(cstring_t mesh_file)
{
    solver_mesh = read_mesh_file(mesh_file);

    if (is_parallel())
        remap_local_mesh(solver_mesh);

    calc_mesh_metrics(solver_mesh);
}

/*******************************************************************************
 * @brief Print mesh information
 * @param mesh
 ******************************************************************************/
void print_mesh_info()
{
    print_mesh(solver_mesh);
}

/*******************************************************************************
 * @brief Read the mesh file
 * @param mesh
 ******************************************************************************/
Mesh_t *read_mesh_file(cstring_t mesh_file)
{
    Mesh_t *mesh = allocate_mesh();
    Partition_t *partition = mesh->partition;
    Vertices_t *vertices = mesh->vertices;
    Cells_t *cells = mesh->cells;
    Boundaries_t *boundaries = mesh->boundaries;
    Faces_t *faces = mesh->faces;
    Regions_t *regions = mesh->regions;

    int rank = get_rank_number();

    hid_t file_id = open_hdf5_file_read_only(mesh_file);
    GET_HDF5_ATTRIBUTE(file_id, "dimension", HDF5Int,
                       &mesh->dimension);
    GET_HDF5_ATTRIBUTE(file_id, "is_partitioned", HDF5Int,
                       &mesh->is_partitioned);

    /* the partition */
    if (is_parallel())
    {
        CHECK_EXPRESSION(mesh->is_partitioned == 1);

        hsize_t dims[2];
        hsize_t offset[2];
        hsize_t count[2];

        int n_partitions = 0;
        int n_partition_cells = 0;
        int n_partition_boundaries = 0;
        int n_partition_faces = 0;
        int n_partition_sends = 0;
        int n_partition_receives = 0;

        hid_t group_id = open_hdf5_group(file_id, "PARTITIONS");
        GET_HDF5_ATTRIBUTE(group_id, "n_partitions", HDF5Int,
                           &n_partitions);
        CHECK_EXPRESSION(n_partitions == get_number_of_procs());

        dims[0] = n_partitions;
        offset[0] = rank;
        count[0] = 1;
        GET_HDF5_DATASET_CHUNK_N(group_id, "n_partition_cells", HDF5Int,
                                 count[0], 1, dims,
                                 offset, count, NULL, NULL,
                                 &n_partition_cells);
        GET_HDF5_DATASET_CHUNK_N(group_id, "n_partition_boundaries", HDF5Int,
                                 count[0], 1, dims,
                                 offset, count, NULL, NULL,
                                 &n_partition_boundaries);
        GET_HDF5_DATASET_CHUNK_N(group_id, "n_partition_faces", HDF5Int,
                                 count[0], 1, dims,
                                 offset, count, NULL, NULL,
                                 &n_partition_faces);
        GET_HDF5_DATASET_CHUNK_N(group_id, "n_partition_sends", HDF5Int,
                                 count[0], 1, dims,
                                 offset, count, NULL, NULL,
                                 &n_partition_sends);
        GET_HDF5_DATASET_CHUNK_N(group_id, "n_partition_receives", HDF5Int,
                                 count[0], 1, dims,
                                 offset, count, NULL, NULL,
                                 &n_partition_receives);

        partition = allocate_partition(mesh,
                                       n_partitions,
                                       n_partition_cells,
                                       n_partition_boundaries,
                                       n_partition_faces,
                                       n_partition_sends,
                                       n_partition_receives);

        dims[0] = n_partitions;
        dims[1] = partition->n_partition_cells;
        offset[0] = rank;
        offset[1] = 0;
        count[0] = 1;
        count[1] = partition->n_partition_cells;
        GET_HDF5_DATASET_CHUNK_N_M(group_id, "partition_cells", HDF5Int,
                                   count, 2, dims,
                                   offset, count, NULL, NULL,
                                   partition->partition_cells);

        dims[0] = n_partitions;
        dims[1] = partition->n_partition_boundaries;
        offset[0] = rank;
        offset[1] = 0;
        count[0] = 1;
        count[1] = partition->n_partition_boundaries;
        GET_HDF5_DATASET_CHUNK_N_M(group_id, "partition_boundaries", HDF5Int,
                                   count, 2, dims,
                                   offset, count, NULL, NULL,
                                   partition->partition_boundaries);

        dims[0] = n_partitions;
        dims[1] = partition->n_partition_faces;
        offset[0] = rank;
        offset[1] = 0;
        count[0] = 1;
        count[1] = partition->n_partition_faces;
        GET_HDF5_DATASET_CHUNK_N_M(group_id, "partition_faces", HDF5Int,
                                   count, 2, dims,
                                   offset, count, NULL, NULL,
                                   partition->partition_faces);

        dims[0] = n_partitions;
        dims[1] = partition->n_partition_sends;
        offset[0] = rank;
        offset[1] = 0;
        count[0] = 1;
        count[1] = partition->n_partition_sends;
        GET_HDF5_DATASET_CHUNK_N_M(group_id, "partition_sends", HDF5Int,
                                   count, 2, dims,
                                   offset, count, NULL, NULL,
                                   partition->partition_sends);
        GET_HDF5_DATASET_CHUNK_N_M(group_id, "partition_sends_pid", HDF5Int,
                                   count, 2, dims,
                                   offset, count, NULL, NULL,
                                   partition->partition_sends_pid);

        dims[0] = n_partitions;
        dims[1] = partition->n_partition_receives;
        offset[0] = rank;
        offset[1] = 0;
        count[0] = 1;
        count[1] = partition->n_partition_receives;
        GET_HDF5_DATASET_CHUNK_N_M(group_id, "partition_receives", HDF5Int,
                                   count, 2, dims,
                                   offset, count, NULL, NULL,
                                   partition->partition_receives);
        GET_HDF5_DATASET_CHUNK_N_M(group_id, "partition_receives_pid", HDF5Int,
                                   count, 2, dims,
                                   offset, count, NULL, NULL,
                                   partition->partition_receives_pid);

        close_hdf5_group(group_id);
    }

    /* the vertices */
    {
        hsize_t dims[2];

        int n_vertices = 0;

        hid_t group_id = open_hdf5_group(file_id, "VERTICES");
        GET_HDF5_ATTRIBUTE(group_id, "n_vertices", HDF5Int, &n_vertices);

        vertices = allocate_vertices(mesh, n_vertices);

        dims[0] = vertices->n_vertices;
        dims[1] = DIM;
        GET_HDF5_DATASET_N_M(group_id, "x", HDF5Double, dims, vertices->x);

        close_hdf5_group(group_id);
    }

    /* the cells */
    {
        hsize_t dims[2];
        hsize_t offset[2];
        hsize_t count[2];

        int n_global_cells = 0;
        int max_cell_vertices = 0;
        int max_cell_faces = 0;

        hid_t group_id = open_hdf5_group(file_id, "CELLS");
        GET_HDF5_ATTRIBUTE(group_id, "n_cells", HDF5Int, &n_global_cells);
        GET_HDF5_ATTRIBUTE(group_id, "max_cell_vertices", HDF5Int,
                           &max_cell_vertices);
        GET_HDF5_ATTRIBUTE(group_id, "max_cell_faces", HDF5Int,
                           &max_cell_faces);

        cells = allocate_cells(mesh, n_global_cells,
                               max_cell_vertices, max_cell_faces);

        dims[0] = cells->n_global_cells;
        count[0] = cells->n_local_cells;
        GET_HDF5_DATASET_SELECT_N(group_id, "id", HDF5Int,
                                  count[0], dims[0],
                                  cells->stride, cells->n_local_cells,
                                  cells->id);

        dims[0] = cells->n_global_cells;
        count[0] = cells->n_local_cells;
        GET_HDF5_DATASET_SELECT_N(group_id, "type", HDF5Int,
                                  count[0], dims[0],
                                  cells->stride, cells->n_local_cells,
                                  cells->type);

        dims[0] = cells->n_global_cells;
        count[0] = cells->n_local_cells;
        GET_HDF5_DATASET_SELECT_N(group_id, "n_vertices", HDF5Int,
                                  count[0], dims[0],
                                  cells->stride, cells->n_local_cells,
                                  cells->n_vertices);

        dims[0] = cells->n_global_cells;
        dims[1] = cells->max_cell_vertices;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = cells->n_local_cells;
        count[1] = cells->max_cell_vertices;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "vertices", HDF5Int,
                                    count, dims, offset[0],
                                    cells->stride, cells->n_local_cells,
                                    cells->vertices);

        dims[0] = cells->n_local_cells;
        count[0] = cells->n_local_cells;
        GET_HDF5_DATASET_SELECT_N(group_id, "n_faces", HDF5Int,
                                  count[0], dims[0],
                                  cells->stride, cells->n_local_cells,
                                  cells->n_faces);

        dims[0] = cells->n_global_cells;
        dims[1] = cells->max_cell_faces;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = cells->n_local_cells;
        count[1] = cells->max_cell_faces;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "faces", HDF5Int,
                                    count, dims, offset[0],
                                    cells->stride, cells->n_local_cells,
                                    cells->faces);

        dims[0] = cells->n_global_cells;
        count[0] = cells->n_local_cells;
        GET_HDF5_DATASET_SELECT_N(group_id, "volume", HDF5Double,
                                  count[0], dims[0],
                                  cells->stride, cells->n_local_cells,
                                  cells->volume);

        dims[0] = cells->n_global_cells;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = cells->n_local_cells;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "x", HDF5Double,
                                    count, dims, offset[0],
                                    cells->stride, cells->n_local_cells,
                                    cells->x);

        dims[0] = cells->n_global_cells;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = cells->n_local_cells;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "dx", HDF5Double,
                                    count, dims, offset[0],
                                    cells->stride, cells->n_local_cells,
                                    cells->dx);

        close_hdf5_group(group_id);
    }

    /* the boundaries */
    {
        hsize_t dims[2];
        hsize_t offset[2];
        hsize_t count[2];

        int n_global_boundaries = 0;
        int max_boundary_vertices = 0;

        hid_t group_id = open_hdf5_group(file_id, "BOUNDARIES");
        GET_HDF5_ATTRIBUTE(group_id, "n_boundaries", HDF5Int,
                           &n_global_boundaries);
        GET_HDF5_ATTRIBUTE(group_id, "max_boundary_vertices", HDF5Int,
                           &max_boundary_vertices);

        boundaries = allocate_boundaries(mesh, n_global_boundaries,
                                         max_boundary_vertices);

        dims[0] = boundaries->n_global_boundaries;
        count[0] = boundaries->n_boundaries;
        GET_HDF5_DATASET_SELECT_N(group_id, "id", HDF5Int,
                                  count[0], dims[0],
                                  boundaries->stride, boundaries->n_boundaries,
                                  boundaries->id);

        dims[0] = boundaries->n_global_boundaries;
        count[0] = boundaries->n_boundaries;
        GET_HDF5_DATASET_SELECT_N(group_id, "type", HDF5Int,
                                  count[0], dims[0],
                                  boundaries->stride, boundaries->n_boundaries,
                                  boundaries->type);

        dims[0] = boundaries->n_global_boundaries;
        count[0] = boundaries->n_boundaries;
        GET_HDF5_DATASET_SELECT_N(group_id, "n_vertices", HDF5Int,
                                  count[0], dims[0],
                                  boundaries->stride, boundaries->n_boundaries,
                                  boundaries->n_vertices);

        dims[0] = boundaries->n_global_boundaries;
        dims[1] = boundaries->max_boundary_vertices;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = boundaries->n_boundaries;
        count[1] = boundaries->max_boundary_vertices;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "vertices", HDF5Int,
                                    count, dims, offset[0],
                                    boundaries->stride,
                                    boundaries->n_boundaries,
                                    boundaries->vertices);

        dims[0] = boundaries->n_global_boundaries;
        count[0] = boundaries->n_boundaries;
        GET_HDF5_DATASET_SELECT_N(group_id, "face", HDF5Int,
                                  count[0], dims[0],
                                  boundaries->stride, boundaries->n_boundaries,
                                  boundaries->face);

        dims[0] = boundaries->n_global_boundaries;
        count[0] = boundaries->n_boundaries;
        GET_HDF5_DATASET_SELECT_N(group_id, "distance", HDF5Double,
                                  count[0], dims[0],
                                  boundaries->stride, boundaries->n_boundaries,
                                  boundaries->distance);

        dims[0] = boundaries->n_global_boundaries;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = boundaries->n_boundaries;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "n", HDF5Double,
                                    count, dims, offset[0],
                                    boundaries->stride,
                                    boundaries->n_boundaries,
                                    boundaries->n);

        dims[0] = boundaries->n_global_boundaries;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = boundaries->n_boundaries;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "t1", HDF5Double,
                                    count, dims, offset[0],
                                    boundaries->stride,
                                    boundaries->n_boundaries,
                                    boundaries->t1);

        dims[0] = boundaries->n_global_boundaries;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = boundaries->n_boundaries;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "t2", HDF5Double,
                                    count, dims, offset[0],
                                    boundaries->stride,
                                    boundaries->n_boundaries,
                                    boundaries->t2);

        close_hdf5_group(group_id);
    }

    /* the faces */
    {
        hsize_t dims[2];
        hsize_t offset[2];
        hsize_t count[2];

        int n_global_faces = 0;
        int max_face_vertices = 0;

        hid_t group_id = open_hdf5_group(file_id, "FACES");
        GET_HDF5_ATTRIBUTE(group_id, "n_faces", HDF5Int, &n_global_faces);
        GET_HDF5_ATTRIBUTE(group_id, "max_face_vertices", HDF5Int,
                           &max_face_vertices);

        faces = allocate_faces(mesh, n_global_faces, max_face_vertices);

        dims[0] = faces->n_global_faces;
        count[0] = faces->n_faces;
        GET_HDF5_DATASET_SELECT_N(group_id, "type", HDF5Int,
                                  count[0], dims[0],
                                  faces->stride, faces->n_faces,
                                  faces->type);

        dims[0] = faces->n_global_faces;
        count[0] = faces->n_faces;
        GET_HDF5_DATASET_SELECT_N(group_id, "n_vertices", HDF5Int,
                                  count[0], dims[0],
                                  faces->stride, faces->n_faces,
                                  faces->n_vertices);

        dims[0] = faces->n_global_faces;
        dims[1] = faces->max_face_vertices;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = faces->n_faces;
        count[1] = faces->max_face_vertices;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "vertices", HDF5Int,
                                    count, dims, offset[0],
                                    faces->stride, faces->n_faces,
                                    faces->vertices);

        dims[0] = faces->n_global_faces;
        dims[1] = FACE_CELLS;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = faces->n_faces;
        count[1] = FACE_CELLS;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "cells", HDF5Int,
                                    count, dims, offset[0],
                                    faces->stride, faces->n_faces,
                                    faces->cells);

        dims[0] = faces->n_global_faces;
        count[0] = faces->n_faces;
        GET_HDF5_DATASET_SELECT_N(group_id, "boundary", HDF5Int,
                                  count[0], dims[0],
                                  faces->stride, faces->n_faces,
                                  faces->boundary);

        dims[0] = faces->n_global_faces;
        count[0] = faces->n_faces;
        GET_HDF5_DATASET_SELECT_N(group_id, "area", HDF5Double,
                                  count[0], dims[0],
                                  faces->stride, faces->n_faces, faces->area);

        dims[0] = faces->n_global_faces;
        count[0] = faces->n_faces;
        GET_HDF5_DATASET_SELECT_N(group_id, "lambda", HDF5Double,
                                  count[0], dims[0],
                                  faces->stride, faces->n_faces, faces->lambda);

        dims[0] = faces->n_global_faces;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = faces->n_faces;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "x", HDF5Double,
                                    count, dims, offset[0],
                                    faces->stride, faces->n_faces, faces->x);

        dims[0] = faces->n_global_faces;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = faces->n_faces;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "n", HDF5Double,
                                    count, dims, offset[0],
                                    faces->stride, faces->n_faces, faces->n);

        dims[0] = faces->n_global_faces;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = faces->n_faces;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "t1", HDF5Double,
                                    count, dims, offset[0],
                                    faces->stride, faces->n_faces, faces->t1);

        dims[0] = faces->n_global_faces;
        dims[1] = DIM;
        offset[0] = 0;
        offset[1] = 0;
        count[0] = faces->n_faces;
        count[1] = DIM;
        GET_HDF5_DATASET_SELECT_N_M(group_id, "t2", HDF5Double,
                                    count, dims, offset[0],
                                    faces->stride, faces->n_faces, faces->t2);

        close_hdf5_group(group_id);
    }

    /* the regions */
    {
        hsize_t dims[2];

        int n_regions = 0;

        hid_t group_id = open_hdf5_group(file_id, "REGIONS");
        GET_HDF5_ATTRIBUTE(group_id, "n_regions", HDF5Int, &n_regions);

        int max_length = get_hdf5_dataset_size(group_id, "name");
        regions = allocate_regions(mesh, n_regions, max_length);

        dims[0] = regions->n_regions;
        GET_HDF5_DATASET_N(group_id, "name", HDF5String,
                           dims[0], regions->name);

        dims[0] = regions->n_regions;
        GET_HDF5_DATASET_N(group_id, "is_boundary", HDF5Int,
                           dims[0], regions->is_boundary);
        close_hdf5_group(group_id);
    }
    close_hdf5_file(file_id);

    return mesh;
}

/*******************************************************************************
 * @brief Remap the mesh (global/local)
 * @param mesh
 ******************************************************************************/
void remap_local_mesh(Mesh_t *mesh)
{
    Partition_t *partition = mesh->partition;
    Cells_t *cells = mesh->cells;
    Boundaries_t *boundaries = mesh->boundaries;
    Faces_t *faces = mesh->faces;

    int *global_cells = ALLOCATE(sizeof(int) * cells->n_global_cells);
    int *global_boundaries =
        ALLOCATE(sizeof(int) * boundaries->n_global_boundaries);
    int *global_faces = ALLOCATE(sizeof(int) * faces->n_global_faces);

    set_value_int_n(-1, cells->n_global_cells, global_cells);
    set_value_int_n(-1, boundaries->n_global_boundaries, global_boundaries);
    set_value_int_n(-1, faces->n_global_faces, global_faces);

    for (int i = 0; i < cells->n_domain_cells; ++i)
        global_cells[partition->partition_cells[i]] = i;

    for (int i = 0; i < partition->n_partition_sends; ++i)
        CHECK_EXPRESSION(global_cells[partition->partition_sends[i]] != -1);

    for (int i = 0; i < partition->n_partition_receives; ++i)
        global_cells[partition->partition_receives[i]] =
            cells->n_domain_cells + i;

    for (int i = 0; i < partition->n_partition_sends; ++i)
        partition->partition_sends[i] =
            global_cells[partition->partition_sends[i]];

    for (int i = 0; i < partition->n_partition_receives; ++i)
        partition->partition_receives[i] =
            global_cells[partition->partition_receives[i]];

    for (int i = 0; i < boundaries->n_boundaries; ++i)
        global_boundaries[partition->partition_boundaries[i]] = i;

    for (int i = 0; i < faces->n_faces; ++i)
        global_faces[partition->partition_faces[i]] = i;

    for (int i = 0; i < cells->n_local_cells; ++i)
    {
        int *cell_faces = &cells->faces[i * cells->max_cell_faces];

        for (int j = 0; j < cells->n_faces[i]; ++j)
            cell_faces[j] = global_faces[cell_faces[j]];
    }

    for (int i = 0; i < boundaries->n_boundaries; ++i)
    {
        int *boundary_face = &boundaries->face[i];

        boundary_face[0] = global_faces[boundary_face[0]];
    }

    for (int i = 0; i < faces->n_faces; ++i)
    {
        int *fc = &faces->cells[i * FACE_CELLS];
        int *face_boundary = &faces->boundary[i];

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
            int *cell_faces = &cells->faces[fc[0] * cells->max_cell_faces];

            cells->n_faces[fc[0]] = 0;
            for (int j = 0; j < n; ++j)
            {
                if (cell_faces[j] < 0)
                    continue;
                cell_faces[cells->n_faces[fc[0]]] = cell_faces[j];
                cells->n_faces[fc[0]] += 1;
            }
        }

        if (fc[1] >= cells->n_domain_cells)
        {
            int n = cells->n_faces[fc[1]];
            int *cell_faces = &cells->faces[fc[1] * cells->max_cell_faces];

            cells->n_faces[fc[1]] = 0;
            for (int j = 0; j < n; ++j)
            {
                if (cell_faces[j] < 0)
                    continue;
                cell_faces[cells->n_faces[fc[1]]] = cell_faces[j];
                cells->n_faces[fc[1]] += 1;
            }
        }
    }

    for (int i = 0; i < partition->n_partitions; ++i)
    {
        if (i == get_rank_number())
            continue;

        int *partition_sends_to =
            &partition->partition_sends_to[i * partition->n_partition_sends];
        for (int j = 0; j < partition->n_partition_sends; ++j)
        {
            if (partition->partition_sends_pid[j] != i)
                continue;
            partition_sends_to[partition->n_partition_sends_to[i]] =
                partition->partition_sends[j];
            partition->n_partition_sends_to[i] += 1;
        }

        int *partition_receives_from =
            &partition->partition_receives_from
                 [i * partition->n_partition_receives];
        for (int j = 0; j < partition->n_partition_receives; ++j)
        {
            if (partition->partition_receives_pid[j] != i)
                continue;
            partition_receives_from[partition->n_partition_receives_from[i]] =
                partition->partition_receives[j];
            partition->n_partition_receives_from[i] += 1;
        }
    }

    DEALLOCATE(global_cells);
    DEALLOCATE(global_boundaries);
    DEALLOCATE(global_faces);
}
