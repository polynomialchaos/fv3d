/*******************************************************************************
 * @file solver_variable.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef SOVLER_VARIABLE_H
#define SOVLER_VARIABLE_H

#include "basec/basec_module.h"

#define SMALL 1e-16  /** Numerical small value */
#define FACE_CELLS 2 /** Number of cells per face (left/right neighbour) */
#define DIM 3        /** Dimension of position, normal and tangential vecotrs */

/*******************************************************************************
 * @brief Boundaries structure
 ******************************************************************************/
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

/*******************************************************************************
 * @brief Cells structure
 ******************************************************************************/
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

/*******************************************************************************
 * @brief Faces structure
 ******************************************************************************/
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

/*******************************************************************************
 * @brief Partition structure
 ******************************************************************************/
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

/*******************************************************************************
 * @brief Regions structure
 ******************************************************************************/
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

/*******************************************************************************
 * @brief Vertices structure
 ******************************************************************************/
typedef struct Vertices
{
    int n_vertices;

    double *x;
} Vertices_t;

/*******************************************************************************
 * @brief Mesh structure
 ******************************************************************************/
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

/*******************************************************************************
 * @brief Variable structure
 ******************************************************************************/
typedef struct Variable
{
    string_t name;
} Variable_t;

/*******************************************************************************
 * @brief Variables structure
 ******************************************************************************/
typedef struct Variables
{
    int n_sol_variables;
    int n_dep_variables;
    int n_tot_variables;

    Variable_t *sol_variables;
    Variable_t *dep_variables;
    Variable_t **tot_variables;
} Variables_t;

/*******************************************************************************
 * @brief Limiter type enumeration
 ******************************************************************************/
typedef enum LimiterType
{
    BarthJespersenn, /** Barth-Jeprsenn's limiter */
    NoLimter         /** No limiter (=1) */
} limiter_type_t;

/*******************************************************************************
 * @brief Reconstruction type enumeration
 ******************************************************************************/
typedef enum ReconstructionType
{
    FirstOrder, /** First-order reconstruction */
    Linear      /** Second-order (linear) reconstruction */
} reconstruction_type_t;

/*******************************************************************************
 * @brief Limiter type enumeration
 ******************************************************************************/
typedef enum TimediscType
{
    Explicit, /** Explicit timestep */
    Implicit  /** Implicit timestep */
} timedisc_type_t;

/*******************************************************************************
 * @brief Explicit scheme enumeration
 ******************************************************************************/
typedef enum ExplicitScheme
{
    EulerExplicit, /** Euler explicit scheme */
    RungeKutta33,  /** Runge-Kutta 3rd order scheme */
    RungeKutta45   /** Runge-Kutta 4th order scheme */
} explicit_scheme_t;

/*******************************************************************************
 * @brief Implicit scheme enumeration
 ******************************************************************************/
typedef enum ImplicitScheme
{
    EulerImplicit, /** Euler implicit scheme */
    BDF2,          /** Backward Differentiation  Bashforth secon oder scheme */
} implicit_scheme_t;

/*******************************************************************************
 * @brief Implicit solver enumeration
 ******************************************************************************/
typedef enum ImplicitSolver
{
    GMRes,    /** GMRes solver */
    BiCGStab, /** BiCGStab solver */
} implicit_solver_t;

#endif /* SOVLER_VARIABLE_H */