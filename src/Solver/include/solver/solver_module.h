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

typedef struct Variable
{
    string_t name;
} Variable_t;

typedef struct Variables
{
    int n_sol_variables;
    int n_dep_variables;
    int n_tot_variables;

    Variable_t *sol_variables;
    Variable_t *dep_variables;
    Variable_t **tot_variables;
} Variables_t;

extern Variables_t *solver_variables;

extern int solver_i_output_data;
extern bool_t solver_do_output_data;

extern int solver_use_restart;
extern int solver_iter_restart;
extern double solver_t_restart;

extern double *solver_residual;

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

typedef void (*void_calc_flux_ft)();
extern void_calc_flux_ft calc_flux_function_pointer;

typedef void (*void_calc_exact_ft)(int id, double t, double *x, double *phi);
extern void_calc_exact_ft calc_exact_function_pointer;

typedef double (*double_limiter_ft)(int i_cell, int i_var, double slope);
extern double_limiter_ft limiter_function_pointer;

typedef void (*void_reconstruction_ft)();
extern void_reconstruction_ft reconstruction_function_pointer;

typedef void (*void_update_ft)(double t);
extern void_update_ft update_function_pointer;

typedef void (*void_update_gradients_ft)();
extern void_update_gradients_ft update_gradients_function_pointer;

extern double *solver_phi_total;
extern double *solver_grad_phi_total_x;
extern double *solver_grad_phi_total_y;
extern double *solver_grad_phi_total_z;

extern double *solver_phi_total_left;
extern double *solver_phi_total_right;

extern double *solver_phi_dt;
extern double *solver_flux;

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

typedef double (*double_calc_timestep_ft)();
extern double_calc_timestep_ft calc_time_step_function_pointer;

extern int solver_explicit_active;
extern int solver_implicit_active;

/*******************************************************************************
 * @brief Add a solution variable
 * @param name
 * @return int
 ******************************************************************************/
int add_sol_variable(string_t name);

/*******************************************************************************
 * @brief Add a dependent variable
 * @param name
 * @return int
 ******************************************************************************/
int add_dep_variable(string_t name);

/*******************************************************************************
 * @brief Calculate the global residual
 * @param dt
 ******************************************************************************/
void calc_global_residual(double dt);

/*******************************************************************************
 * @brief Create a file header
 * @param file_name
 ******************************************************************************/
void create_file_header(cstring_t file_name);

/*******************************************************************************
 * @brief The finite volume time derivative
 * @param t
 ******************************************************************************/
void finite_volume_time_derivative(double t);

/*******************************************************************************
 * @brief Finalize analyze
 ******************************************************************************/
void free_analyze();

/*******************************************************************************
 * @brief Free equation
 ******************************************************************************/
void free_equation();

/*******************************************************************************
 * @brief Free explicit
 ******************************************************************************/
void free_explicit();

/*******************************************************************************
 * @brief Free finite_volume
 ******************************************************************************/
void free_finite_volume();

/*******************************************************************************
 * @brief Free implicit timedisc
 ******************************************************************************/
void free_implicit();

/*******************************************************************************
 * @brief Free limiter
 ******************************************************************************/
void free_limiter();

/*******************************************************************************
 * @brief Free mesh
 ******************************************************************************/
void free_mesh();

/*******************************************************************************
 * @brief Free reconstruction
 ******************************************************************************/
void free_reconstruction();

/*******************************************************************************
 * @brief Free output
 ******************************************************************************/
void free_output();

/*******************************************************************************
 * @brief Free restart
 ******************************************************************************/
void free_restart();

/*******************************************************************************
 * @brief Free solver
 ******************************************************************************/
void free_solver();

/*******************************************************************************
 * @brief Free timedisc
 ******************************************************************************/
void free_timedisc();

/*******************************************************************************
 * @brief Return the simulation title
 ******************************************************************************/
cstring_t get_simulation_title();

/*******************************************************************************
 * @brief Initialize analyze
 ******************************************************************************/
void init_analyze();

/*******************************************************************************
 * @brief Initialize equation
 ******************************************************************************/
void init_equation();

/*******************************************************************************
 * @brief Initialize explicit
 * @param explicit_scheme
 ******************************************************************************/
void init_explicit(explicit_scheme_t explicit_scheme);

/*******************************************************************************
 * @brief Initialize finite_volume
 ******************************************************************************/
void init_finite_volume();

/*******************************************************************************
 * @brief Initialize implicit timedisc
 * @param implicit_scheme
 * @param implicit_solver
 * @param max_iter_inner
 * @param tolerance_lsoe
 * @param max_iter_lsoe
 * @param max_krylov_dims
 * @param max_krylov_restarts
 ******************************************************************************/
void init_implicit(implicit_scheme_t implicit_scheme,
                   implicit_solver_t implicit_solver,
                   int max_iter_inner,
                   double tolerance_lsoe, int max_iter_lsoe,
                   int max_krylov_dims, int max_krylov_restarts);

/*******************************************************************************
 * @brief Initialize limiter
 * @param limiter_type
 ******************************************************************************/
void init_limiter(limiter_type_t limiter_type);

/*******************************************************************************
 * @brief Initialize mesh
 * @param mesh_file
 ******************************************************************************/
void init_mesh(cstring_t mesh_file);

/*******************************************************************************
 * @brief Initialize reconstruction
 * @param reconstruction_type
 ******************************************************************************/
void init_reconstruction(reconstruction_type_t reconstruction_type);

/*******************************************************************************
 * @brief Initialize output
 * @param i_output_data
 ******************************************************************************/
void init_output(int i_output_data);

/*******************************************************************************
 * @brief Initialize restart
 * @param use_restart
 * @param restart_file
 ******************************************************************************/
void init_restart(bool_t use_restart, cstring_t restart_file);

/*******************************************************************************
 * @brief Initialize solver
 * @param title
 ******************************************************************************/
void init_solver(cstring_t title);

/*******************************************************************************
 * @brief Initialize timedisc
 * @param timedisc_type
 * @param max_iter
 * @param is_transient
 * @param abort_residual
 * @param t_start
 * @param t_end
 ******************************************************************************/
void init_timedisc(timedisc_type_t timedisc_type, int max_iter,
                   bool_t is_transient, double abort_residual,
                   double t_start, double t_end);

/*******************************************************************************
 * @brief Print mesh information
 * @param mesh
 ******************************************************************************/
void print_mesh_info();

/*******************************************************************************
 * @brief Print the variables
 ******************************************************************************/
void print_variables();

/*******************************************************************************
 * @brief Read restart data
 * @param restart_file
 ******************************************************************************/
void read_restart_data(cstring_t restart_file);

/*******************************************************************************
 * @brief Set the viscous timestep flag
 * @param is_viscous_dt
 ******************************************************************************/
void set_viscous_dt(bool_t is_viscous_dt);

/*******************************************************************************
 * @brief Time discretizazion routine
 ******************************************************************************/
void timedisc();

/*******************************************************************************
 * @brief Write data to file
 * @param iter
 * @param t
 ******************************************************************************/
void write_output(int iter, double t);

#endif /* SOLVER_MODULE_H */