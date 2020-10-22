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

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
typedef struct Partition_t
{
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
    double *x;
} Vertices_t;

typedef struct Cells
{
    int *id;
    int *type;
    int *n_vertices;
    int *vertices;
    int *n_faces;
    int *faces;
    double *x;
    double *volume;
    double *ds;
} Cells_t;

typedef struct Boundaries
{
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
} Faces_t;

typedef struct Regions
{
    string_t *name;
    double *phi;
    int *type;
    int *function_id;
} Regions_t;

typedef struct Mesh
{
    int dimension;
    int is_partitioned;

    int n_partitions;
    int n_partition_cells;
    int n_partition_boundaries;
    int n_partition_faces;
    int n_partition_sends;
    int n_partition_receives;
    Partition_t partition;

    int n_vertices;
    Vertices_t vertices;

    int n_global_cells;
    int n_local_cells;
    int n_cells;
    int max_cell_vertices;
    int max_cell_faces;
    Cells_t cells;

    int n_global_boundaries;
    int n_local_boundaries;
    int n_boundaries;
    int max_boundary_vertices;
    Boundaries_t boundaries;

    int n_global_faces;
    int n_local_faces;
    int n_faces;
    int max_face_vertices;
    Faces_t faces;

    int n_regions;
    Regions_t regions;

    // double *dist_cell_1;
    // double *dist_cell_2;

    // double total_volume;
    // int flow_region;

    // int n_internal_faces;
    // int *internal_faces;
    // int n_boundary_faces;
    // int *boundary_faces;
} Mesh_t;

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_define();

void free_mesh( Mesh_t *mesh );
void allocate_partition( Partition_t *partition, int n_partitions,
    int n_cells, int n_boundaries, int n_faces, int n_sends, int n_receives );
void allocate_vertices( Vertices_t *vertices, int n_vertices );
void allocate_cells( Cells_t *cells, int n_cells, int max_vertices, int max_faces );
void allocate_boundaries( Boundaries_t *boundaries, int n_boundaries, int max_vertices );
void allocate_faces( Faces_t *faces, int n_faces, int max_vertices );
void allocate_regions( Regions_t *regions, int n_regions, int length );

#endif /* MESH_PRIVATE_H */