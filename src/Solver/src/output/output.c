/*******************************************************************************
 * @file output.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <string.h>
#include "output_private.h"
// #include "fv3d_private.h"
// #include "mesh/mesh_module.h"
//
// #include "fv/fv_module.h"
// #include "timedisc/implicit/implicit_module.h"

int solver_i_output_data = -1;
bool_t solver_do_output_data = BFLS;

/*******************************************************************************
 * @brief Create a file header
 ******************************************************************************/
void create_file_header(cstring_t file_name)
{
    hid_t file_id = create_hdf5_file(file_name);

    SET_HDF5_ATTRIBUTE(file_id, "n_tot_variables",
                       HDF5Int, &solver_variables->n_tot_variables);
    SET_HDF5_ATTRIBUTE(file_id, "n_sol_variables",
                       HDF5Int, &solver_variables->n_sol_variables);

    {
        size_t max_len = 0;
        for (int i = 0; i < solver_variables->n_tot_variables; ++i)
            max_len = MAX(max_len, strlen(
                                       solver_variables->tot_variables[i]->name));

        string_t *tmp = allocate_hdf5_string_buffer(
            solver_variables->n_tot_variables, max_len + 1, NULL);
        for (int i = 0; i < solver_variables->n_tot_variables; ++i)
            strcpy(tmp[i], solver_variables->tot_variables[i]->name);

        hsize_t dims[2] = {solver_variables->n_tot_variables, max_len + 1};
        SET_HDF5_DATASET_N(file_id, "tot_variables", HDF5String, tmp, dims[0]);

        deallocate_hdf5_string_buffer(tmp);
        DEALLOCATE(tmp);
    }

    {
        size_t max_len = 0;
        for (int i = 0; i < solver_variables->n_sol_variables; ++i)
            max_len = MAX(max_len, strlen(
                                       (&solver_variables->sol_variables[i])->name));

        string_t *tmp = allocate_hdf5_string_buffer(
            solver_variables->n_sol_variables, max_len + 1, NULL);
        for (int i = 0; i < solver_variables->n_sol_variables; ++i)
            strcpy(tmp[i], (&solver_variables->sol_variables[i])->name);

        hsize_t dims[2] = {solver_variables->n_sol_variables, max_len + 1};
        SET_HDF5_DATASET_N(file_id, "sol_variables", HDF5String, tmp, dims[0]);

        deallocate_hdf5_string_buffer(tmp);
        DEALLOCATE(tmp);
    }

    close_hdf5_file(file_id);
}

/*******************************************************************************
 * @brief Free output
 ******************************************************************************/
void free_output()
{
}

/*******************************************************************************
 * @brief Initialize output
 * @param i_output_data
 ******************************************************************************/
void init_output(int i_output_data)
{
    solver_i_output_data = i_output_data;
    solver_do_output_data = (solver_i_output_data != 0);
}

/*******************************************************************************
 * @brief Write data to file
 * @param iter
 * @param t
 ******************************************************************************/
void write_output(int iter, double t)
{
    // Cells_t *cells = solver_mesh->cells;
    // int n_domain_cells = cells->n_domain_cells;
    // int n_global_cells = cells->n_global_cells;
    // int n_tot_variables = solver_variables->n_tot_variables;
    // int n_sol_variables = solver_variables->n_sol_variables;

    char file_prefix[256];
    sprintf(file_prefix, ".%09d.h5", iter);
    cstring_t title = get_simulation_title();
    string_t output_file = allocate_strcat(title, file_prefix);
    create_file_header(output_file);

    hid_t file_id = open_hdf5_file(output_file);

    SET_HDF5_ATTRIBUTE(file_id, "iter", HDF5Int, &iter);
    SET_HDF5_ATTRIBUTE(file_id, "t", HDF5Double, &t);

    // if (is_parallel())
    // {
    //     {
    //         hsize_t dims[2] = {n_global_cells, n_sol_variables};
    //         hsize_t offset[2] = {0, 0};
    //         hsize_t count[2] = {n_domain_cells, n_sol_variables};
    //         SET_HDF5_DATASET_SELECT_N_M(file_id, "phi_total", HDF5Double, phi_total,
    //                                     count, dims, offset[0], cells->stride, n_domain_cells);
    //     }

    //     {
    //         hsize_t dims[2] = {n_global_cells, n_sol_variables};
    //         hsize_t offset[2] = {0, 0};
    //         hsize_t count[2] = {n_domain_cells, n_sol_variables};
    //         SET_HDF5_DATASET_SELECT_N_M(file_id, "phi_dt", HDF5Double, phi_dt,
    //                                     count, dims, offset[0], cells->stride, n_domain_cells);
    //     }

    //     if (n_bdf_stages > 0)
    //     {
    //         SET_HDF5_ATTRIBUTE(file_id, "n_stages", HDF5Int, &n_bdf_stages);

    //         for (int i_stage = 0; i_stage < n_bdf_stages; ++i_stage)
    //         {
    //             char stage_string[256];
    //             sprintf(stage_string, "%d", i_stage);
    //             string_t tmp = allocate_strcat("phi_old:", stage_string);

    //             hsize_t dims[2] = {n_global_cells, n_sol_variables};
    //             hsize_t offset[2] = {0, 0};
    //             hsize_t count[2] = {n_domain_cells, n_sol_variables};
    //             SET_HDF5_DATASET_SELECT_N_M(file_id, tmp, HDF5Double, phi_old[i_stage],
    //                                         count, dims, offset[0], cells->stride, n_domain_cells);

    //             DEALLOCATE(tmp);
    //         }
    //     }
    // }
    // else
    // {
    //     {
    //         hsize_t dims[2] = {n_domain_cells, n_tot_variables};
    //         SET_HDF5_DATASET_N_M(file_id, "phi_total", HDF5Double, phi_total, dims);
    //     }

    //     {
    //         hsize_t dims[2] = {n_domain_cells, n_sol_variables};
    //         SET_HDF5_DATASET_N_M(file_id, "phi_dt", HDF5Double, phi_dt, dims);
    //     }

    //     if (n_bdf_stages > 0)
    //     {
    //         SET_HDF5_ATTRIBUTE(file_id, "n_stages", HDF5Int, &n_bdf_stages);

    //         for (int i_stage = 0; i_stage < n_bdf_stages; ++i_stage)
    //         {
    //             char stage_string[256];
    //             sprintf(stage_string, "%d", i_stage);
    //             string_t tmp = allocate_strcat("phi_old:", stage_string);

    //             hsize_t dims[2] = {n_domain_cells, n_sol_variables};

    //             SET_HDF5_DATASET_N_M(file_id, tmp, HDF5Double, phi_old[i_stage], dims);

    //             DEALLOCATE(tmp);
    //         }
    //     }
    // }

    close_hdf5_file(file_id);
    DEALLOCATE(output_file);
}