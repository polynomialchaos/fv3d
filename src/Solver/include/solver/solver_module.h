/*******************************************************************************
 * @file solver_module.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#ifndef SOLVER_MODULE_H
#define SOLVER_MODULE_H

#include "basec/basec_module.h"
#include "solver/solver_type.h"

extern Mesh_t *solver_mesh;           /** Global mesh */
extern Variables_t *solver_variables; /** Global variables */
extern data_t *solver_data;           /** Global solver data */

/*******************************************************************************
 * @brief Add a solution variable
 * @param name
 * @return int
 ******************************************************************************/
int add_sol_variable(cstring_t name);

/*******************************************************************************
 * @brief Add a dependent variable
 * @param name
 * @return int
 ******************************************************************************/
int add_dep_variable(cstring_t name);

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
 * @param time
 ******************************************************************************/
void finite_volume_time_derivative(double time);

/*******************************************************************************
 * @brief Free analyze
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
 * @brief Return is_explicit flag
 * bool_t
 ******************************************************************************/
bool_t is_explicit();

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
 * @brief Set the flux calculation routine
 * @param fun_ptr
 ******************************************************************************/
void set_calc_flux(void_calc_flux_ft fun_ptr);

/*******************************************************************************
 * @brief Set the timestep calculation routine
 * @param fun_ptr
 ******************************************************************************/
void set_calc_timestep(double_calc_timestep_ft fun_ptr);

/*******************************************************************************
 * @brief Set the exact function routine
 * @param fun_ptr
 ******************************************************************************/
void set_exact_function(void_calc_exact_ft fun_ptr);

/*******************************************************************************
 * @brief Set the update routine
 * @param fun_ptr
 ******************************************************************************/
void set_update(void_update_ft fun_ptr);

/*******************************************************************************
 * @brief Set the update gradients routine
 * @param fun_ptr
 ******************************************************************************/
void set_update_gradients(void_update_gradients_ft fun_ptr);

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