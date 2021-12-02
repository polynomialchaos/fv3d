/*******************************************************************************
 * @file solver_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef SOLVER_MODULE_H
#define SOLVER_MODULE_H

#include "basec/basec_module.h"

#define SMALL 1e-16 /** Numerical small value */

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

extern Mesh_t *solver_mesh; /** Global mesh */

/*******************************************************************************
 * @brief Free mesh
 ******************************************************************************/
void free_mesh();

/*******************************************************************************
 * @brief Free solver
 ******************************************************************************/
void free_solver();

/*******************************************************************************
 * @brief Return the simulation title
 ******************************************************************************/
cstring_t get_simulation_title();

/*******************************************************************************
 * @brief Initialize mesh
 ******************************************************************************/
void init_mesh(cstring_t mesh_file);

/*******************************************************************************
 * @brief Initialize solver
 ******************************************************************************/
void init_solver(cstring_t title);

/*******************************************************************************
 * @brief Print mesh information
 * @param mesh
 ******************************************************************************/
void print_mesh_info();

#endif /* SOLVER_MODULE_H */