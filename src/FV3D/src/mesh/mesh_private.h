//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef MESH_PRIVATE_H
#define MESH_PRIVATE_H

#include "fv3d_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------
#define FACE_CELLS 2
#define DIM 3

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
typedef struct Partition_t
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
    int n_cells;
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
} Cells_t;

typedef struct Boundaries
{
    int n_global_boundaries;
    int n_local_boundaries;
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
} Boundaries_t;

typedef struct Faces
{
    int n_global_faces;
    int n_local_faces;
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

    double total_volume;
} Mesh_t;

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_define();

Mesh_t *allocate_mesh();
void print_mesh_info( Mesh_t *mesh );
void deallocate_mesh( Mesh_t **mesh );

void allocate_partition( Partition_t *partition, int n_partitions,
    int n_cells, int n_boundaries, int n_faces, int n_sends, int n_receives );
void allocate_vertices( Vertices_t *vertices, int n_vertices );
void allocate_cells( Cells_t *cells, int n_cells, int max_vertices, int max_faces );
void allocate_boundaries( Boundaries_t *boundaries, int n_boundaries, int max_vertices );
void allocate_faces( Faces_t *faces, int n_faces, int max_vertices );
void allocate_regions( Regions_t *regions, int n_regions, int max_name_length );

#endif /* MESH_PRIVATE_H */