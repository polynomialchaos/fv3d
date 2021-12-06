/*******************************************************************************
 * @file restart.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "restart_private.h"

int solver_use_restart = 0;
int solver_iter_restart = 0;
double solver_t_restart = 0.0;

double *phi_total_restart = NULL;
double *phi_dt_restart = NULL;
double **phi_old_restart = NULL;

int n_stages_restart = 0;

/*******************************************************************************
 * @brief Free restart
 ******************************************************************************/
void free_restart()
{
    DEALLOCATE(phi_total_restart);
    DEALLOCATE(phi_dt_restart);
    DEALLOCATE(phi_old_restart);
}

/*******************************************************************************
 * @brief Initialize restart
 * @param use_restart
 * @param restart_file
 ******************************************************************************/
void init_restart(bool_t use_restart, cstring_t restart_file)
{
    solver_use_restart = use_restart;

    if (solver_use_restart == BTRU)
        read_restart_data(restart_file);
}

/*******************************************************************************
 * @brief Read restart data
 * @param restart_file
 ******************************************************************************/
void read_restart_data(cstring_t restart_file)
{
    Cells_t *cells = solver_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_tot_variables = solver_variables->n_tot_variables;
    int n_sol_variables = solver_variables->n_sol_variables;

    hid_t file_id = open_hdf5_file(restart_file);

    GET_HDF5_ATTRIBUTE(file_id, "iter", HDF5Int, &solver_iter_restart);
    GET_HDF5_ATTRIBUTE(file_id, "t", HDF5Double, &solver_t_restart);

    phi_total_restart = ALLOCATE(sizeof(double) * n_tot_variables * n_domain_cells);
    phi_dt_restart = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);

    if (is_parallel())
    {
        {
            hsize_t dims[2] = {n_domain_cells, n_tot_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_tot_variables};

            GET_HDF5_DATASET_SELECT_N_M(file_id, "solver_phi_total", HDF5Double,
                                        count, dims, offset[0], cells->stride, n_domain_cells, phi_total_restart);
        }

        {
            hsize_t dims[2] = {n_domain_cells, n_sol_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_sol_variables};

            GET_HDF5_DATASET_SELECT_N_M(file_id, "solver_phi_dt", HDF5Double,
                                        count, dims, offset[0], cells->stride, n_domain_cells, phi_dt_restart);
        }

        // if (n_bdf_stages > 0)
        // {
        //     GET_HDF5_ATTRIBUTE(file_id, "n_stages", HDF5Int, &n_stages_restart);
        //     CHECK_EXPRESSION(n_stages_restart == n_bdf_stages);

        //     phi_old_restart = ALLOCATE(sizeof(double *) * n_stages_restart);
        //     for (int i = 0; i < n_stages_restart; ++i)
        //         phi_old_restart[i] = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);

        //     for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
        //     {
        //         char iter_string[256];
        //         sprintf(iter_string, "%d", i_stage);
        //         string_t tmp = allocate_strcat("phi_old:", iter_string);

        //         hsize_t dims[2] = {n_domain_cells, n_sol_variables};
        //         hsize_t offset[2] = {0, 0};
        //         hsize_t count[2] = {n_domain_cells, n_sol_variables};

        //         GET_HDF5_DATASET_SELECT_N_M(file_id, tmp, HDF5Double,
        //                                     count, dims, offset[0], cells->stride, n_domain_cells, phi_old_restart[i_stage]);

        //         DEALLOCATE(tmp);
        //     }
        // }
    }
    else
    {
        {
            hsize_t dims[2] = {n_domain_cells, n_tot_variables};
            GET_HDF5_DATASET_N_M(file_id, "solver_phi_total", HDF5Double, dims, phi_total_restart);
        }

        {
            hsize_t dims[2] = {n_domain_cells, n_sol_variables};
            GET_HDF5_DATASET_N_M(file_id, "solver_phi_dt", HDF5Double, dims, phi_dt_restart);
        }

        // if (n_bdf_stages > 0)
        // {
        //     GET_HDF5_ATTRIBUTE(file_id, "n_stages", HDF5Int, &n_stages_restart);
        //     CHECK_EXPRESSION(n_stages_restart == n_bdf_stages);

        //     phi_old_restart = ALLOCATE(sizeof(double *) * n_stages_restart);
        //     for (int i = 0; i < n_stages_restart; ++i)
        //         phi_old_restart[i] = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);

        //     for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
        //     {
        //         char iter_string[256];
        //         sprintf(iter_string, "%d", i_stage);
        //         string_t tmp = allocate_strcat("phi_old:", iter_string);

        //         hsize_t dims[2] = {n_domain_cells, n_sol_variables};

        //         GET_HDF5_DATASET_N_M(file_id, tmp, HDF5Double, dims, phi_old_restart[i_stage]);

        //         DEALLOCATE(tmp);
        //     }
        // }
    }

    close_hdf5_file(file_id);

    copy_n(phi_total_restart, n_tot_variables * n_domain_cells, solver_phi_total);
    copy_n(phi_dt_restart, n_sol_variables * n_domain_cells, solver_phi_dt);

    // for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
    // {
    //     copy_n(phi_old_restart[i_stage], n_sol_variables * n_domain_cells, phi_old[i_stage]);
    // }

    DEALLOCATE(phi_total_restart);
    DEALLOCATE(phi_dt_restart);
    DEALLOCATE(phi_old_restart);
}