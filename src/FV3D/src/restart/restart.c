/*******************************************************************************
 * @file restart.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "fv3d_private.h"
#include "restart_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"
#include "output/output_module.h"
#include "fv/fv_module.h"
#include "timedisc/implicit/implicit_module.h"

int use_restart = 0;
int restart_iter = 0;

int iter_restart = 0;
double t_restart = 0;

double *phi_total_restart = NULL;
double *phi_dt_restart = NULL;
double **phi_old_restart = NULL;

int n_stages_restart = 0;

/*******************************************************************************
 * @brief Read restart data
 ******************************************************************************/
void read_restart_data()
{
    Cells_t *cells = global_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_tot_variables = all_variables->n_tot_variables;
    int n_sol_variables = all_variables->n_sol_variables;

    char file_prefix[256];
    sprintf(file_prefix, ".%09d.h5", restart_iter);
    string_t output_file = allocate_strcat(title, file_prefix);
    create_file_header(output_file);

    hid_t file_id = open_hdf5_file(output_file);

    GET_HDF5_ATTRIBUTE(file_id, "iter", HDF5Int, &iter_restart);
    GET_HDF5_ATTRIBUTE(file_id, "t", HDF5Double, &t_restart);

    phi_total_restart = ALLOCATE(sizeof(double) * n_tot_variables * n_domain_cells);
    phi_dt_restart = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);

    if (is_parallel())
    {
        {
            hsize_t dims[2] = {n_domain_cells, n_tot_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_tot_variables};

            GET_HDF5_DATASET_SELECT_N_M(file_id, "phi_total", HDF5Double,
                                        count, dims, offset[0], cells->stride, n_domain_cells, phi_total_restart);
        }

        {
            hsize_t dims[2] = {n_domain_cells, n_sol_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_sol_variables};

            GET_HDF5_DATASET_SELECT_N_M(file_id, "phi_dt", HDF5Double,
                                        count, dims, offset[0], cells->stride, n_domain_cells, phi_dt_restart);
        }

        if (n_bdf_stages > 0)
        {
            GET_HDF5_ATTRIBUTE(file_id, "n_stages", HDF5Int, &n_stages_restart);
            CHECK_EXPRESSION(n_stages_restart == n_bdf_stages);

            phi_old_restart = ALLOCATE(sizeof(double *) * n_stages_restart);
            for (int i = 0; i < n_stages_restart; ++i)
                phi_old_restart[i] = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);

            for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
            {
                char iter_string[256];
                sprintf(iter_string, "%d", i_stage);
                string_t tmp = allocate_strcat("phi_old:", iter_string);

                hsize_t dims[2] = {n_domain_cells, n_sol_variables};
                hsize_t offset[2] = {0, 0};
                hsize_t count[2] = {n_domain_cells, n_sol_variables};

                GET_HDF5_DATASET_SELECT_N_M(file_id, tmp, HDF5Double,
                                            count, dims, offset[0], cells->stride, n_domain_cells, phi_old_restart[i_stage]);

                DEALLOCATE(tmp);
            }
        }
    }
    else
    {
        {
            hsize_t dims[2] = {n_domain_cells, n_tot_variables};
            GET_HDF5_DATASET_N_M(file_id, "phi_total", HDF5Double, dims, phi_total_restart);
        }

        {
            hsize_t dims[2] = {n_domain_cells, n_sol_variables};
            GET_HDF5_DATASET_N_M(file_id, "phi_dt", HDF5Double, dims, phi_dt_restart);
        }

        if (n_bdf_stages > 0)
        {
            GET_HDF5_ATTRIBUTE(file_id, "n_stages", HDF5Int, &n_stages_restart);
            CHECK_EXPRESSION(n_stages_restart == n_bdf_stages);

            phi_old_restart = ALLOCATE(sizeof(double *) * n_stages_restart);
            for (int i = 0; i < n_stages_restart; ++i)
                phi_old_restart[i] = ALLOCATE(sizeof(double) * n_sol_variables * n_domain_cells);

            for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
            {
                char iter_string[256];
                sprintf(iter_string, "%d", i_stage);
                string_t tmp = allocate_strcat("phi_old:", iter_string);

                hsize_t dims[2] = {n_domain_cells, n_sol_variables};

                GET_HDF5_DATASET_N_M(file_id, tmp, HDF5Double, dims, phi_old_restart[i_stage]);

                DEALLOCATE(tmp);
            }
        }
    }

    close_hdf5_file(file_id);

    copy_n(phi_total_restart, n_tot_variables * n_domain_cells, phi_total);
    copy_n(phi_dt_restart, n_sol_variables * n_domain_cells, phi_dt);

    for (int i_stage = 0; i_stage < n_stages_restart; ++i_stage)
    {
        copy_n(phi_old_restart[i_stage], n_sol_variables * n_domain_cells, phi_old[i_stage]);
    }

    DEALLOCATE(phi_total_restart);
    DEALLOCATE(phi_dt_restart);
    DEALLOCATE(phi_old_restart);
    DEALLOCATE(output_file);
}

/*******************************************************************************
 * @brief Define restart
 ******************************************************************************/
void restart_define()
{
    REGISTER_INITIALIZE_ROUTINE(restart_initialize);
    REGISTER_FINALIZE_ROUTINE(restart_finalize);

    SET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart, "The flag to start from restart", NULL, 0);
    SET_PARAMETER("Restart/restart_iter", DigitParameter, &restart_iter, "The restart iteration", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize restart
 ******************************************************************************/
void restart_finalize()
{
    DEALLOCATE(phi_total_restart);
    DEALLOCATE(phi_dt_restart);
    DEALLOCATE(phi_old_restart);
}

/*******************************************************************************
 * @brief Initialize restart
 ******************************************************************************/
void restart_initialize()
{
    GET_PARAMETER("Restart/use_restart", LogicalParameter, &use_restart);
    GET_PARAMETER("Restart/restart_iter", DigitParameter, &restart_iter);

    if (use_restart == 1)
        read_restart_data();
}
