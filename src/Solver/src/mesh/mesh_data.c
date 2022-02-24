/*******************************************************************************
 * @file mesh_data.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include <string.h>
#include "mesh_private.h"

/*******************************************************************************
 * @brief Allocate the boundaries
 * @param mesh
 * @param n_boundaries
 * @param max_boundary_vertices
 * @return Boundaries_t*
 ******************************************************************************/
Boundaries_t *allocate_boundaries(Mesh_t *mesh, int n_boundaries,
                                  int max_boundary_vertices)
{
    mesh->boundaries = BM_ALLOCATE(sizeof(Boundaries_t));
    Boundaries_t *boundaries = mesh->boundaries;

    boundaries->n_global_boundaries = n_boundaries;
    boundaries->n_boundaries = n_boundaries;
    boundaries->max_boundary_vertices = max_boundary_vertices;

    if (is_parallel())
    {
        boundaries->n_boundaries = mesh->partition->n_partition_boundaries;
        n_boundaries = boundaries->n_boundaries;
    }

    boundaries->id = BM_ALLOCATE(sizeof(int) * n_boundaries);
    boundaries->type = BM_ALLOCATE(sizeof(int) * n_boundaries);
    boundaries->n_vertices = BM_ALLOCATE(sizeof(int) * n_boundaries);
    boundaries->vertices =
        BM_ALLOCATE(sizeof(int) * max_boundary_vertices * n_boundaries);
    boundaries->face = BM_ALLOCATE(sizeof(int) * n_boundaries);
    boundaries->distance = BM_ALLOCATE(sizeof(double) * n_boundaries);
    boundaries->n = BM_ALLOCATE(sizeof(double) * DIM * n_boundaries);
    boundaries->t1 = BM_ALLOCATE(sizeof(double) * DIM * n_boundaries);
    boundaries->t2 = BM_ALLOCATE(sizeof(double) * DIM * n_boundaries);

    set_value_int_n(0, n_boundaries, boundaries->n_vertices);

    boundaries->stride = BM_ALLOCATE(sizeof(hsize_t) * n_boundaries);
    if (is_parallel())
    {
        Partition_t *partition = mesh->partition;
        int *partition_boundaries = partition->partition_boundaries;
        int n_partition_boundaries = partition->n_partition_boundaries;

        for (int i = 0; i < n_partition_boundaries; ++i)
            boundaries->stride[i] = partition_boundaries[i];
    }
    else
    {
        for (int i = 0; i < n_boundaries; ++i)
            boundaries->stride[i] = i;
    }

    return boundaries;
}

/*******************************************************************************
 * @brief Allocate the cells
 * @param mesh
 * @param n_cells
 * @param max_cell_vertices
 * @param max_cell_faces
 * @return Cells_t*
 ******************************************************************************/
Cells_t *allocate_cells(Mesh_t *mesh, int n_cells,
                        int max_cell_vertices, int max_cell_faces)
{
    mesh->cells = BM_ALLOCATE(sizeof(Cells_t));
    Cells_t *cells = mesh->cells;

    cells->n_global_cells = n_cells;
    cells->n_local_cells = n_cells;
    cells->n_domain_cells = n_cells;
    cells->max_cell_vertices = max_cell_vertices;
    cells->max_cell_faces = max_cell_faces;

    if (is_parallel())
    {
        cells->n_domain_cells = mesh->partition->n_partition_cells;
        cells->n_local_cells = mesh->partition->n_partition_cells +
                               mesh->partition->n_partition_receives;
        n_cells = cells->n_local_cells;
    }

    cells->id = BM_ALLOCATE(sizeof(int) * n_cells);
    cells->type = BM_ALLOCATE(sizeof(int) * n_cells);
    cells->n_vertices = BM_ALLOCATE(sizeof(int) * n_cells);
    cells->vertices = BM_ALLOCATE(sizeof(int) * max_cell_vertices * n_cells);
    cells->n_faces = BM_ALLOCATE(sizeof(int) * n_cells);
    cells->faces = BM_ALLOCATE(sizeof(int) * max_cell_faces * n_cells);
    cells->x = BM_ALLOCATE(sizeof(double) * DIM * n_cells);
    cells->volume = BM_ALLOCATE(sizeof(double) * n_cells);
    cells->dx = BM_ALLOCATE(sizeof(double) * DIM * n_cells);

    set_value_int_n(0, n_cells, cells->n_vertices);
    set_value_int_n(0, n_cells, cells->n_faces);

    cells->stride = BM_ALLOCATE(sizeof(hsize_t) * n_cells);
    if (is_parallel())
    {
        Partition_t *partition = mesh->partition;
        int *partition_cells = partition->partition_cells;
        int *partition_receives = partition->partition_receives;
        int n_partition_cells = partition->n_partition_cells;
        int n_partition_receives = partition->n_partition_receives;

        for (int i = 0; i < n_partition_cells; ++i)
            cells->stride[i] = partition_cells[i];

        for (int i = 0; i < n_partition_receives; ++i)
            cells->stride[n_partition_cells + i] = partition_receives[i];
    }
    else
    {
        for (int i = 0; i < n_cells; ++i)
            cells->stride[i] = i;
    }

    return cells;
}

/*******************************************************************************
 * @brief Allocate the faces
 * @param mesh
 * @param n_faces
 * @param max_face_vertices
 * @return Faces_t*
 ******************************************************************************/
Faces_t *allocate_faces(Mesh_t *mesh, int n_faces, int max_face_vertices)
{
    mesh->faces = BM_ALLOCATE(sizeof(Faces_t));
    Faces_t *faces = mesh->faces;

    faces->n_global_faces = n_faces;
    faces->n_faces = n_faces;
    faces->max_face_vertices = max_face_vertices;

    if (is_parallel())
    {
        faces->n_faces = mesh->partition->n_partition_faces;
        n_faces = faces->n_faces;
    }

    faces->type = BM_ALLOCATE(sizeof(int) * n_faces);
    faces->n_vertices = BM_ALLOCATE(sizeof(int) * n_faces);
    faces->vertices = BM_ALLOCATE(sizeof(int) * max_face_vertices * n_faces);
    faces->cells = BM_ALLOCATE(sizeof(int) * FACE_CELLS * n_faces);
    faces->boundary = BM_ALLOCATE(sizeof(int) * n_faces);
    faces->area = BM_ALLOCATE(sizeof(double) * n_faces);
    faces->lambda = BM_ALLOCATE(sizeof(double) * n_faces);
    faces->x = BM_ALLOCATE(sizeof(double) * DIM * n_faces);
    faces->n = BM_ALLOCATE(sizeof(double) * DIM * n_faces);
    faces->t1 = BM_ALLOCATE(sizeof(double) * DIM * n_faces);
    faces->t2 = BM_ALLOCATE(sizeof(double) * DIM * n_faces);

    set_value_int_n(0, n_faces, faces->n_vertices);

    faces->n_internal_faces = 0;
    faces->n_boundary_faces = 0;
    faces->dist_cell_1 = NULL;
    faces->dist_cell_2 = NULL;
    faces->internal_faces = NULL;
    faces->boundary_faces = NULL;

    faces->stride = BM_ALLOCATE(sizeof(hsize_t) * n_faces);
    if (is_parallel())
    {
        Partition_t *partition = mesh->partition;
        int *partition_faces = partition->partition_faces;
        int n_partition_faces = partition->n_partition_faces;

        for (int i = 0; i < n_partition_faces; ++i)
            faces->stride[i] = partition_faces[i];
    }
    else
    {
        for (int i = 0; i < n_faces; ++i)
            faces->stride[i] = i;
    }

    return faces;
}

/*******************************************************************************
 * @brief Allocate the mesh
 * @return Mesh_t*
 ******************************************************************************/
Mesh_t *allocate_mesh()
{
    Mesh_t *tmp = BM_ALLOCATE(sizeof(Mesh_t));

    tmp->dimension = 0;
    tmp->is_partitioned = 0;

    tmp->partition = NULL;
    tmp->vertices = NULL;
    tmp->cells = NULL;
    tmp->boundaries = NULL;
    tmp->faces = NULL;
    tmp->regions = NULL;

    tmp->local_volume = 0.0;
    tmp->global_volume = 0.0;

    return tmp;
}

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_partitions
 * @param n_partition_cells
 * @param n_partition_boundaries
 * @param n_partition_faces
 * @param n_partition_sends
 * @param n_partition_receives
 * @return Partition_t*
 ******************************************************************************/
Partition_t *allocate_partition(Mesh_t *mesh,
                                int n_partitions,
                                int n_partition_cells,
                                int n_partition_boundaries,
                                int n_partition_faces,
                                int n_partition_sends,
                                int n_partition_receives)
{
    mesh->partition = BM_ALLOCATE(sizeof(Partition_t));
    Partition_t *partition = mesh->partition;

    partition->n_partitions = n_partitions;
    partition->n_partition_cells = n_partition_cells;
    partition->n_partition_boundaries = n_partition_boundaries;
    partition->n_partition_faces = n_partition_faces;
    partition->n_partition_sends = n_partition_sends;
    partition->n_partition_receives = n_partition_receives;

    partition->partition_cells = BM_ALLOCATE(sizeof(int) * n_partition_cells);
    partition->partition_boundaries =
        BM_ALLOCATE(sizeof(int) * n_partition_boundaries);
    partition->partition_faces = BM_ALLOCATE(sizeof(int) * n_partition_faces);
    partition->partition_sends = BM_ALLOCATE(sizeof(int) * n_partition_sends);
    partition->partition_sends_pid =
        BM_ALLOCATE(sizeof(int) * n_partition_sends);
    partition->partition_receives =
        BM_ALLOCATE(sizeof(int) * n_partition_receives);
    partition->partition_receives_pid =
        BM_ALLOCATE(sizeof(int) * n_partition_receives);

    partition->n_partition_sends_to = BM_ALLOCATE(sizeof(int) * n_partitions);
    partition->partition_sends_to =
        BM_ALLOCATE(sizeof(int) * n_partition_sends * n_partitions);
    partition->n_partition_receives_from =
        BM_ALLOCATE(sizeof(int) * n_partitions);
    partition->partition_receives_from =
        BM_ALLOCATE(sizeof(int) * n_partition_receives * n_partitions);

    set_value_int_n(0, n_partitions, partition->n_partition_sends_to);
    set_value_int_n(0, n_partitions, partition->n_partition_receives_from);

    return partition;
}

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_regions
 * @param max_name_length
 * @return Regions_t*
 ******************************************************************************/
Regions_t *allocate_regions(Mesh_t *mesh, int n_regions, int max_name_length)
{
    mesh->regions = BM_ALLOCATE(sizeof(Regions_t));
    Regions_t *regions = mesh->regions;

    regions->n_regions = n_regions;
    regions->max_name_length = max_name_length;

    regions->name =
        allocate_hdf5_string_buffer(n_regions, max_name_length, NULL);
    regions->is_boundary = BM_ALLOCATE(sizeof(int) * n_regions);

    regions->type = NULL;
    regions->function_id = NULL;
    regions->phi_total = NULL;

    return regions;
}

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_vertices
 * @return Vertices_t*
 ******************************************************************************/
Vertices_t *allocate_vertices(Mesh_t *mesh, int n_vertices)
{
    mesh->vertices = BM_ALLOCATE(sizeof(Vertices_t));
    Vertices_t *vertices = mesh->vertices;

    vertices->n_vertices = n_vertices;

    vertices->x = BM_ALLOCATE(sizeof(double) * n_vertices * DIM);

    return vertices;
}

/*******************************************************************************
 * @brief Deallocate the boundaries
 * @param boundaries
 ******************************************************************************/
void deallocate_boundaries(Boundaries_t *boundaries)
{
    if (boundaries == NULL)
        return;

    BM_DEALLOCATE(boundaries->id);
    BM_DEALLOCATE(boundaries->type);
    BM_DEALLOCATE(boundaries->n_vertices);
    BM_DEALLOCATE(boundaries->vertices);
    BM_DEALLOCATE(boundaries->face);
    BM_DEALLOCATE(boundaries->distance);
    BM_DEALLOCATE(boundaries->n);
    BM_DEALLOCATE(boundaries->t1);
    BM_DEALLOCATE(boundaries->t2);

    BM_DEALLOCATE(boundaries->stride);
}

/*******************************************************************************
 * @brief Deallocate the cells
 * @param cells
 ******************************************************************************/
void deallocate_cells(Cells_t *cells)
{
    if (cells == NULL)
        return;

    BM_DEALLOCATE(cells->id);
    BM_DEALLOCATE(cells->type);
    BM_DEALLOCATE(cells->n_vertices);
    BM_DEALLOCATE(cells->vertices);
    BM_DEALLOCATE(cells->n_faces);
    BM_DEALLOCATE(cells->faces);
    BM_DEALLOCATE(cells->x);
    BM_DEALLOCATE(cells->volume);
    BM_DEALLOCATE(cells->dx);

    BM_DEALLOCATE(cells->stride);
}

/*******************************************************************************
 * @brief Deallocate the faces
 * @param faces
 ******************************************************************************/
void deallocate_faces(Faces_t *faces)
{
    if (faces == NULL)
        return;

    BM_DEALLOCATE(faces->type);
    BM_DEALLOCATE(faces->n_vertices);
    BM_DEALLOCATE(faces->vertices);
    BM_DEALLOCATE(faces->cells);
    BM_DEALLOCATE(faces->boundary);
    BM_DEALLOCATE(faces->area);
    BM_DEALLOCATE(faces->lambda);
    BM_DEALLOCATE(faces->x);
    BM_DEALLOCATE(faces->n);
    BM_DEALLOCATE(faces->t1);
    BM_DEALLOCATE(faces->t2);

    BM_DEALLOCATE(faces->stride);

    BM_DEALLOCATE(faces->dist_cell_1);
    BM_DEALLOCATE(faces->dist_cell_2);
    BM_DEALLOCATE(faces->internal_faces);
    BM_DEALLOCATE(faces->boundary_faces);
}

/*******************************************************************************
 * @brief Deallocate the mesh
 * @param mesh
 ******************************************************************************/
void deallocate_mesh(Mesh_t *mesh)
{
    if (mesh == NULL)
        return;

    deallocate_partition(mesh->partition);
    BM_DEALLOCATE(mesh->partition);

    deallocate_vertices(mesh->vertices);
    BM_DEALLOCATE(mesh->vertices);

    deallocate_cells(mesh->cells);
    BM_DEALLOCATE(mesh->cells);

    deallocate_boundaries(mesh->boundaries);
    BM_DEALLOCATE(mesh->boundaries);

    deallocate_faces(mesh->faces);
    BM_DEALLOCATE(mesh->faces);

    deallocate_regions(mesh->regions);
    BM_DEALLOCATE(mesh->regions);
}

/*******************************************************************************
 * @brief Deallocate the partition
 * @param partition
 ******************************************************************************/
void deallocate_partition(Partition_t *partition)
{
    if (partition == NULL)
        return;

    BM_DEALLOCATE(partition->partition_cells);
    BM_DEALLOCATE(partition->partition_boundaries);
    BM_DEALLOCATE(partition->partition_faces);
    BM_DEALLOCATE(partition->partition_sends);
    BM_DEALLOCATE(partition->partition_sends_pid);
    BM_DEALLOCATE(partition->partition_receives);
    BM_DEALLOCATE(partition->partition_receives_pid);

    BM_DEALLOCATE(partition->n_partition_sends_to);
    BM_DEALLOCATE(partition->partition_sends_to);
    BM_DEALLOCATE(partition->n_partition_receives_from);
    BM_DEALLOCATE(partition->partition_receives_from);
}

/*******************************************************************************
 * @brief Deallocate the regions
 * @param regions
 ******************************************************************************/
void deallocate_regions(Regions_t *regions)
{
    if (regions == NULL)
        return;

    deallocate_hdf5_string_buffer(regions->name);
    BM_DEALLOCATE(regions->name);
    BM_DEALLOCATE(regions->is_boundary);
}

/*******************************************************************************
 * @brief Deallocate the vertices
 * @param vertices
 ******************************************************************************/
void deallocate_vertices(Vertices_t *vertices)
{
    if (vertices == NULL)
        return;

    BM_DEALLOCATE(vertices->x);
}