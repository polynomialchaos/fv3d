//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "fv3d_private.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"
#include "output_module.h"
#include "fv/fv_module.h"
#include "timedisc/implicit/implicit_module.h"



int i_output_data = -1;
int do_output_data = 0;
string_t output_file = NULL;

void output_initialize();
void output_finalize();

void output_define()
{
    REGISTER_INITIALIZE_ROUTINE(output_initialize);
    REGISTER_FINALIZE_ROUTINE(output_finalize);

    SET_PARAMETER("Output/i_output_data", DigitParameter, &i_output_data,
                  "The output file frequency  (-1 ... first/solutions/last, 0 ... disable)", NULL, 0);
}

void output_initialize()
{
    GET_PARAMETER("Output/i_output_data", DigitParameter, &i_output_data);
    do_output_data = (i_output_data != 0);

    output_file = allocate_strcat(title, ".h5");
}

void output_finalize()
{
    DEALLOCATE(output_file);
}

void create_file_header()
{
    hid_t file_id = create_hdf5_file(output_file);

    SET_HDF5_ATTRIBUTE(file_id, "n_tot_variables", HDF5Int, &all_variables->n_tot_variables);
    SET_HDF5_ATTRIBUTE(file_id, "n_sol_variables", HDF5Int, &all_variables->n_sol_variables);

    {
        size_t max_len = 0;
        for (int i = 0; i < all_variables->n_tot_variables; ++i)
            max_len = MAX(max_len, strlen(all_variables->tot_variables[i]->name));

        string_t *tmp = allocate_hdf5_string_buffer(all_variables->n_tot_variables, max_len + 1, NULL);
        for (int i = 0; i < all_variables->n_tot_variables; ++i)
            strcpy(tmp[i], all_variables->tot_variables[i]->name);

        hsize_t dims[2] = {all_variables->n_tot_variables, max_len + 1};
        SET_HDF5_DATASET_N(file_id, "tot_variables", HDF5String, tmp, dims[0]);

        deallocate_hdf5_string_buffer(tmp);
        DEALLOCATE(tmp);
    }

    {
        size_t max_len = 0;
        for (int i = 0; i < all_variables->n_sol_variables; ++i)
            max_len = MAX(max_len, strlen((&all_variables->sol_variables[i])->name));

        string_t *tmp = allocate_hdf5_string_buffer(all_variables->n_sol_variables, max_len + 1, NULL);
        for (int i = 0; i < all_variables->n_sol_variables; ++i)
            strcpy(tmp[i], (&all_variables->sol_variables[i])->name);

        hsize_t dims[2] = {all_variables->n_sol_variables, max_len + 1};
        SET_HDF5_DATASET_N(file_id, "sol_variables", HDF5String, tmp, dims[0]);

        deallocate_hdf5_string_buffer(tmp);
        DEALLOCATE(tmp);
    }

    hid_t group_id = create_hdf5_group(file_id, "SOLUTIONS");
    close_hdf5_group(group_id);

    close_hdf5_file(file_id);
}

void write_output(int iter, double t)
{
    Cells_t *cells = global_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_global_cells = cells->n_global_cells;
    int n_tot_variables = all_variables->n_tot_variables;
    int n_sol_variables = all_variables->n_sol_variables;

    if (is_valid_hdf5_file(output_file) == 0)
        create_file_header();

    char iter_string[10];
    sprintf(iter_string, "%09d", iter);

    hid_t file_id = open_hdf5_file(output_file);

    hid_t group_id = open_hdf5_group(file_id, "SOLUTIONS");

    hid_t solution_id = create_hdf5_group(group_id, iter_string);

    SET_HDF5_ATTRIBUTE(solution_id, "iter", HDF5Int, &iter);
    SET_HDF5_ATTRIBUTE(solution_id, "t", HDF5Double, &t);

    if (is_parallel())
    {
        {
            hsize_t dims[2] = {n_global_cells, n_sol_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_sol_variables};
            SET_HDF5_DATASET_SELECT_N_M(solution_id, "phi_total", HDF5Double, phi_total,
                                        count, dims, offset[0], cells->stride, n_domain_cells);
        }

        {
            hsize_t dims[2] = {n_global_cells, n_sol_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_sol_variables};
            SET_HDF5_DATASET_SELECT_N_M(solution_id, "phi_dt", HDF5Double, phi_dt,
                                        count, dims, offset[0], cells->stride, n_domain_cells);
        }

        if (n_bdf_stages > 0)
        {
            SET_HDF5_ATTRIBUTE(solution_id, "n_stages", HDF5Int, &n_bdf_stages);

            for (int i_stage = 0; i_stage < n_bdf_stages; ++i_stage)
            {
                char stage_string[256];
                sprintf(stage_string, "%d", i_stage);
                string_t tmp = allocate_strcat("phi_old:", stage_string);

                hsize_t dims[2] = {n_global_cells, n_sol_variables};
                hsize_t offset[2] = {0, 0};
                hsize_t count[2] = {n_domain_cells, n_sol_variables};
                SET_HDF5_DATASET_SELECT_N_M(solution_id, tmp, HDF5Double, phi_old[i_stage],
                                            count, dims, offset[0], cells->stride, n_domain_cells);

                DEALLOCATE(tmp);
            }
        }
    }
    else
    {
        {
            hsize_t dims[2] = {n_domain_cells, n_tot_variables};
            SET_HDF5_DATASET_N_M(solution_id, "phi_total", HDF5Double, phi_total, dims);
        }

        {
            hsize_t dims[2] = {n_domain_cells, n_sol_variables};
            SET_HDF5_DATASET_N_M(solution_id, "phi_dt", HDF5Double, phi_dt, dims);
        }

        if (n_bdf_stages > 0)
        {
            SET_HDF5_ATTRIBUTE(solution_id, "n_stages", HDF5Int, &n_bdf_stages);

            for (int i_stage = 0; i_stage < n_bdf_stages; ++i_stage)
            {
                char stage_string[256];
                sprintf(stage_string, "%d", i_stage);
                string_t tmp = allocate_strcat("phi_old:", stage_string);

                hsize_t dims[2] = {n_domain_cells, n_sol_variables};

                SET_HDF5_DATASET_N_M(solution_id, tmp, HDF5Double, phi_old[i_stage], dims);

                DEALLOCATE(tmp);
            }
        }
    }

    close_hdf5_group(solution_id);

    close_hdf5_group(group_id);

    if (exists_hdf5_link(file_id, "SOLUTION"))
        delete_hdf5_link(file_id, "SOLUTION");
    string_t tmp = allocate_strcat("SOLUTIONS/", iter_string);
    create_hdf5_soft_link(file_id, "SOLUTION", tmp);
    DEALLOCATE(tmp);

    close_hdf5_file(file_id);
}