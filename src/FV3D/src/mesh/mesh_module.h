/*******************************************************************************
 * @file mesh_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef MESH_MODULE_H
#define MESH_MODULE_H

#include "fv3d/fv3d_module.h"

#define FACE_CELLS 2
#define DIM 3

typedef struct Partition
{
    int n_partitions;
    int n_partition_cells;
    int n_partition_boundaries;
    int n_partition_faces;
    int n_partition_sends;
    int n_partition_receives;

    int *partition_cells;
    int *partition_boundaries;
    int *partition_faces;
    int *partition_sends;
    int *partition_sends_pid;
    int *partition_receives;
    int *partition_receives_pid;

    int *n_partition_sends_to;
    int *partition_sends_to;
    int *n_partition_receives_from;
    int *partition_receives_from;
} Partition_t;

typedef struct Vertices
{
    int n_vertices;

    double *x;
} Vertices_t;

typedef struct Cells
{
    int n_global_cells;
    int n_local_cells;
    int n_domain_cells;
    int max_cell_vertices;
    int max_cell_faces;

    int *id;
    int *type;
    int *n_vertices;
    int *vertices;
    int *n_faces;
    int *faces;
    double *x;
    double *volume;
    double *dx;

    hsize_t *stride;
} Cells_t;

typedef struct Boundaries
{
    int n_global_boundaries;
    int n_boundaries;
    int max_boundary_vertices;

    int *id;
    int *type;
    int *n_vertices;
    int *vertices;
    int *face;
    double *distance;
    double *n;
    double *t1;
    double *t2;

    hsize_t *stride;
} Boundaries_t;

typedef struct Faces
{
    int n_global_faces;
    int n_faces;
    int max_face_vertices;

    int *type;
    int *n_vertices;
    int *vertices;
    int *cells;
    int *boundary;
    double *area;
    double *lambda;
    double *x;
    double *n;
    double *t1;
    double *t2;

    hsize_t *stride;

    double *dist_cell_1;
    double *dist_cell_2;

    int n_internal_faces;
    int *internal_faces;
    int n_boundary_faces;
    int *boundary_faces;
} Faces_t;

typedef struct Regions
{
    int n_regions;
    int max_name_length;

    string_t *name;
    int *is_boundary;

    int flow_region;

    int *type;
    int *function_id;
    double *phi_total;
} Regions_t;

typedef struct Mesh
{
    int dimension;
    int is_partitioned;

    Partition_t *partition;
    Vertices_t *vertices;
    Cells_t *cells;
    Boundaries_t *boundaries;
    Faces_t *faces;
    Regions_t *regions;

    double local_volume;
    double global_volume;
} Mesh_t;

extern Mesh_t *global_mesh;

/*******************************************************************************
 * @brief Allocate the boundaries
 * @param mesh
 * @param n_boundaries
 * @param max_boundary_vertices
 * @return Boundaries_t*
 ******************************************************************************/
Boundaries_t *allocate_boundaries(Mesh_t *mesh, int n_boundaries, int max_boundary_vertices);

/*******************************************************************************
 * @brief Allocate the cells
 * @param mesh
 * @param n_cells
 * @param max_cell_vertices
 * @param max_cell_faces
 * @return Cells_t*
 ******************************************************************************/
Cells_t *allocate_cells(Mesh_t *mesh, int n_cells, int max_cell_vertices, int max_cell_faces);

/*******************************************************************************
 * @brief Allocate the faces
 * @param mesh
 * @param n_faces
 * @param max_face_vertices
 * @return Faces_t*
 ******************************************************************************/
Faces_t *allocate_faces(Mesh_t *mesh, int n_faces, int max_face_vertices);

/*******************************************************************************
 * @brief Allocate the mesh
 * @return Mesh_t*
 ******************************************************************************/
Mesh_t *allocate_mesh();

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
Partition_t *allocate_partition(Mesh_t *mesh, int n_partitions, int n_partition_cells, int n_partition_boundaries,
                                int n_partition_faces, int n_partition_sends, int n_partition_receives);

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_regions
 * @param max_name_length
 * @return Regions_t*
 ******************************************************************************/
Regions_t *allocate_regions(Mesh_t *mesh, int n_regions, int max_name_length);

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_vertices
 * @return Vertices_t*
 ******************************************************************************/
Vertices_t *allocate_vertices(Mesh_t *mesh, int n_vertices);

/*******************************************************************************
 * @brief Deallocate the boundaries
 * @param boundaries
 ******************************************************************************/
void deallocate_boundaries(Boundaries_t **boundaries);

/*******************************************************************************
 * @brief Deallocate the cells
 * @param cells
 ******************************************************************************/
void deallocate_cells(Cells_t **cells);

/*******************************************************************************
 * @brief Deallocate the faces
 * @param faces
 ******************************************************************************/
void deallocate_faces(Faces_t **faces);

/*******************************************************************************
 * @brief Deallocate the mesh
 * @param mesh
 ******************************************************************************/
void deallocate_mesh(Mesh_t **mesh);

/*******************************************************************************
 * @brief Deallocate the partition
 * @param partition
 ******************************************************************************/
void deallocate_partition(Partition_t **partition);

/*******************************************************************************
 * @brief Deallocate the regions
 * @param regions
 ******************************************************************************/
void deallocate_regions(Regions_t **regions);

/*******************************************************************************
 * @brief Deallocate the vertices
 * @param vertices
 ******************************************************************************/
void deallocate_vertices(Vertices_t **vertices);

/*******************************************************************************
 * @brief Define mesh
 ******************************************************************************/
void mesh_define();

/*******************************************************************************
 * @brief Finalize mesh
 ******************************************************************************/
void mesh_finalize();

/*******************************************************************************
 * @brief Initialize mesh
 ******************************************************************************/
void mesh_initialize();

/*******************************************************************************
 * @brief Print the boundaries
 * @param boundaries
 ******************************************************************************/
void print_boundaries(Boundaries_t *boundaries);

/*******************************************************************************
 * @brief Print the cellls
 * @param celss
 ******************************************************************************/
void print_cells(Cells_t *cells);

/*******************************************************************************
 * @brief Print the faces
 * @param faces
 ******************************************************************************/
void print_faces(Faces_t *faces);

/*******************************************************************************
 * @brief Print the mesh
 * @param mesh
 ******************************************************************************/
void print_mesh_info(Mesh_t *mesh);

/*******************************************************************************
 * @brief Print the partition
 * @param partition
 ******************************************************************************/
void print_partition(Partition_t *partition);

/*******************************************************************************
 * @brief Print the regions
 * @param regions
 ******************************************************************************/
void print_regions(Regions_t *regions);

/*******************************************************************************
 * @brief Print the vertices
 * @param vertices
 ******************************************************************************/
void print_vertices(Vertices_t *vertices);

/*******************************************************************************
 * @brief Read the mesh file
 * @param mesh
 ******************************************************************************/
void read_mesh_file(Mesh_t *mesh);

/*******************************************************************************
 * @brief Remap the mesh (global/local)
 * @param mesh
 ******************************************************************************/
void remap_local_mesh(Mesh_t *mesh);

#endif /* MESH_MODULE_H */