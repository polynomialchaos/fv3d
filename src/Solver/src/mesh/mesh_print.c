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

    PRINTF("BOUNDARIES\n");
    PRINTF("n_global_boundaries   = %d\n", boundaries->n_global_boundaries);
    PRINTF("n_boundaries          = %d\n", boundaries->n_boundaries);
    PRINTF("max_boundary_vertices = %d\n", boundaries->max_boundary_vertices);
}

/*******************************************************************************
 * @brief Print the cellls
 * @param celss
 ******************************************************************************/
void print_cells(Cells_t *cells)
{
    if (cells == NULL)
        return;

    PRINTF("CELLS\n");
    PRINTF("n_global_cells    = %d\n", cells->n_global_cells);
    PRINTF("n_local_cells     = %d\n", cells->n_local_cells);
    PRINTF("n_domain_cells    = %d\n", cells->n_domain_cells);
    PRINTF("max_cell_vertices = %d\n", cells->max_cell_vertices);
    PRINTF("max_cell_faces    = %d\n", cells->max_cell_faces);
}

/*******************************************************************************
 * @brief Print the faces
 * @param faces
 ******************************************************************************/
void print_faces(Faces_t *faces)
{
    if (faces == NULL)
        return;

    PRINTF("FACES\n");
    PRINTF("n_global_faces    = %d\n", faces->n_global_faces);
    PRINTF("n_faces           = %d\n", faces->n_faces);
    PRINTF("max_face_vertices = %d\n", faces->max_face_vertices);
}

/*******************************************************************************
 * @brief Print the mesh
 * @param mesh
 ******************************************************************************/
void print_mesh(Mesh_t *mesh)
{
    if (mesh == NULL)
        return;

    PRINTF("MESH\n");
    PRINTF("dimension      = %d\n", mesh->dimension);
    PRINTF("is_partitioned = %d\n", mesh->is_partitioned);

    print_partition(mesh->partition);
    print_vertices(mesh->vertices);
    print_cells(mesh->cells);
    print_boundaries(mesh->boundaries);
    print_faces(mesh->faces);
    print_regions(mesh->regions);

    PRINTF("local_volume  = %e\n", mesh->local_volume);
    PRINTF("global_volume = %e\n", mesh->global_volume);
}

/*******************************************************************************
 * @brief Print the partition
 * @param partition
 ******************************************************************************/
void print_partition(Partition_t *partition)
{
    if (partition == NULL)
        return;

    PRINTF("PARTITION\n");
    PRINTF("n_partitions           = %d\n", partition->n_partitions);
    PRINTF("n_partition_cells      = %d\n", partition->n_partition_cells);
    PRINTF("n_partition_boundaries = %d\n", partition->n_partition_boundaries);
    PRINTF("n_partition_faces      = %d\n", partition->n_partition_faces);
    PRINTF("n_partition_sends      = %d\n", partition->n_partition_sends);
    PRINTF("n_partition_receives   = %d\n", partition->n_partition_receives);
}

/*******************************************************************************
 * @brief Print the regions
 * @param regions
 ******************************************************************************/
void print_regions(Regions_t *regions)
{
    if (regions == NULL)
        return;

    PRINTF("REGIONS\n");
    PRINTF("n_regions   = %d\n", regions->n_regions);
    PRINTF("flow_region = %d\n", regions->flow_region);

    for (int i = 0; i < regions->n_regions; ++i)
        PRINTF("%d: %s\n", i, regions->name[i]);
}

/*******************************************************************************
 * @brief Print the vertices
 * @param vertices
 ******************************************************************************/
void print_vertices(Vertices_t *vertices)
{
    if (vertices == NULL)
        return;

    PRINTF("VERTICES\n");
    PRINTF("n_vertices = %d\n", vertices->n_vertices);
}
