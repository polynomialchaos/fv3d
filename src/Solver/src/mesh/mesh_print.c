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
 * @brief Print the boundaries
 * @param boundaries
 ******************************************************************************/
void print_boundaries(Boundaries_t *boundaries)
{
    if (boundaries == NULL)
        return;

    BM_PRINT("BOUNDARIES\n");
    BM_PRINT("n_global_boundaries   = %d\n", boundaries->n_global_boundaries);
    BM_PRINT("n_boundaries          = %d\n", boundaries->n_boundaries);
    BM_PRINT("max_boundary_vertices = %d\n", boundaries->max_boundary_vertices);
}

/*******************************************************************************
 * @brief Print the cellls
 * @param celss
 ******************************************************************************/
void print_cells(Cells_t *cells)
{
    if (cells == NULL)
        return;

    BM_PRINT("CELLS\n");
    BM_PRINT("n_global_cells    = %d\n", cells->n_global_cells);
    BM_PRINT("n_local_cells     = %d\n", cells->n_local_cells);
    BM_PRINT("n_domain_cells    = %d\n", cells->n_domain_cells);
    BM_PRINT("max_cell_vertices = %d\n", cells->max_cell_vertices);
    BM_PRINT("max_cell_faces    = %d\n", cells->max_cell_faces);
}

/*******************************************************************************
 * @brief Print the faces
 * @param faces
 ******************************************************************************/
void print_faces(Faces_t *faces)
{
    if (faces == NULL)
        return;

    BM_PRINT("FACES\n");
    BM_PRINT("n_global_faces    = %d\n", faces->n_global_faces);
    BM_PRINT("n_faces           = %d\n", faces->n_faces);
    BM_PRINT("max_face_vertices = %d\n", faces->max_face_vertices);
}

/*******************************************************************************
 * @brief Print the mesh
 * @param mesh
 ******************************************************************************/
void print_mesh(Mesh_t *mesh)
{
    if (mesh == NULL)
        return;

    BM_PRINT("MESH\n");
    BM_PRINT("dimension      = %d\n", mesh->dimension);
    BM_PRINT("is_partitioned = %d\n", mesh->is_partitioned);

    print_partition(mesh->partition);
    print_vertices(mesh->vertices);
    print_cells(mesh->cells);
    print_boundaries(mesh->boundaries);
    print_faces(mesh->faces);
    print_regions(mesh->regions);

    BM_PRINT("local_volume  = %e\n", mesh->local_volume);
    BM_PRINT("global_volume = %e\n", mesh->global_volume);
}

/*******************************************************************************
 * @brief Print the partition
 * @param partition
 ******************************************************************************/
void print_partition(Partition_t *partition)
{
    if (partition == NULL)
        return;

    BM_PRINT("PARTITION\n");
    BM_PRINT("n_partitions           = %d\n", partition->n_partitions);
    BM_PRINT("n_partition_cells      = %d\n", partition->n_partition_cells);
    BM_PRINT("n_partition_boundaries = %d\n",
             partition->n_partition_boundaries);
    BM_PRINT("n_partition_faces      = %d\n", partition->n_partition_faces);
    BM_PRINT("n_partition_sends      = %d\n", partition->n_partition_sends);
    BM_PRINT("n_partition_receives   = %d\n", partition->n_partition_receives);
}

/*******************************************************************************
 * @brief Print the regions
 * @param regions
 ******************************************************************************/
void print_regions(Regions_t *regions)
{
    if (regions == NULL)
        return;

    BM_PRINT("REGIONS\n");
    BM_PRINT("n_regions   = %d\n", regions->n_regions);
    BM_PRINT("flow_region = %d\n", regions->flow_region);

    for (int i = 0; i < regions->n_regions; ++i)
        BM_PRINT("%d: %s\n", i, regions->name[i]);
}

/*******************************************************************************
 * @brief Print the vertices
 * @param vertices
 ******************************************************************************/
void print_vertices(Vertices_t *vertices)
{
    if (vertices == NULL)
        return;

    BM_PRINT("VERTICES\n");
    BM_PRINT("n_vertices = %d\n", vertices->n_vertices);
}
