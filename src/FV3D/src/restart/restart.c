//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "restart_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"
#include "output/output_module.h"
#include "fv/fv_module.h"
#include "timedisc/implicit/implicit_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int use_restart = 0;

int iter_restart = 0;
double t_restart = 0;

double *phi_total_restart = NULL;
double *phi_dt_restart = NULL;
double **phi_old_restart = NULL;

int n_stages_restart = 0;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void restart_initialize();
void restart_finalize();

void read_restart_data();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void restart_define()
{
    register_initialize_routine(restart_initialize);
    register_finalize_routine(restart_finalize);

    set_parameter("Restart/use_restart", ParameterBool, &use_restart, "The flag to start from restart", NULL, 0);
}

void restart_initialize()
{
    get_parameter("Restart/use_restart", ParameterBool, &use_restart);

    if (use_restart == 1)
        read_restart_data();
    else
        create_file_header();
}

void restart_finalize()
{
    deallocate(phi_total_restart);
    deallocate(phi_dt_restart);
    deallocate(phi_old_restart);
}

void read_restart_data()
{
    Cells_t *cells = global_mesh->cells;
    int n_domain_cells = cells->n_domain_cells;
    int n_tot_variables = all_variables->n_tot_variables;
    int n_sol_variables = all_variables->n_sol_variables;

    hid_t file_id = open_hdf5_file(output_file);

    hid_t last_id = open_hdf5_group(file_id, "SOLUTION");

    get_hdf5_attribute(last_id, "iter", HDF5Int, &iter_restart);
    get_hdf5_attribute(last_id, "t", HDF5Double, &t_restart);

    phi_total_restart = allocate(sizeof(double) * n_tot_variables * n_domain_cells);
    phi_dt_restart = allocate(sizeof(double) * n_sol_variables * n_domain_cells);

    if (get_is_parallel())
    {
        {
            hsize_t dims[2] = {n_domain_cells, n_tot_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_tot_variables};

            get_hdf5_dataset_select_n_m(last_id, "phi_total", HDF5Double, phi_total_restart,
                                        2, dims, NULL, offset, count, cells->stride, n_domain_cells);
        }

        {
            hsize_t dims[2] = {n_domain_cells, n_sol_variables};
            hsize_t offset[2] = {0, 0};
            hsize_t count[2] = {n_domain_cells, n_sol_variables};

            get_hdf5_dataset_select_n_m(last_id, "phi_dt", HDF5Double, phi_dt_restart,
                                        2, dims, NULL, offset, count, cells->stride, n_domain_cells);
        }

        if (n_bdf_stages > 0)
        {
            get_hdf5_attribute(last_id, "n_stages", HDF5Int, &n_stages_restart);
            check_error((n_stages_restart == n_bdf_stages));

            phi_old_restart = allocate(sizeof(double *) * n_stages_restart);
            for (int i = 0; i < n_stages_restart; i++)
                phi_old_restart[i] = allocate(sizeof(double) * n_sol_variables * n_domain_cells);

            for (int i_stage = 0; i_stage < n_stages_restart; i_stage++)
            {
                char iter_string[10];
                sprintf(iter_string, "%d", i_stage);
                string_t tmp = allocate_strcat("phi_old:", iter_string);

                hsize_t dims[2] = {n_domain_cells, n_sol_variables};
                hsize_t offset[2] = {0, 0};
                hsize_t count[2] = {n_domain_cells, n_sol_variables};

                get_hdf5_dataset_select_n_m(last_id, tmp, HDF5Double, phi_old_restart[i_stage],
                                            2, dims, NULL, offset, count, cells->stride, n_domain_cells);

                deallocate(tmp);
            }
        }
    }
    else
    {
        {
            hsize_t dims[2] = {n_domain_cells, n_tot_variables};
            get_hdf5_dataset_n_m(last_id, "phi_total", HDF5Double, phi_total_restart, 2, dims);
        }

        {
            hsize_t dims[2] = {n_domain_cells, n_sol_variables};
            get_hdf5_dataset_n_m(last_id, "phi_dt", HDF5Double, phi_dt_restart, 2, dims);
        }

        if (n_bdf_stages > 0)
        {
            get_hdf5_attribute(last_id, "n_stages", HDF5Int, &n_stages_restart);
            check_error((n_stages_restart == n_bdf_stages));

            phi_old_restart = allocate(sizeof(double *) * n_stages_restart);
            for (int i = 0; i < n_stages_restart; i++)
                phi_old_restart[i] = allocate(sizeof(double) * n_sol_variables * n_domain_cells);

            for (int i_stage = 0; i_stage < n_stages_restart; i_stage++)
            {
                char iter_string[10];
                sprintf(iter_string, "%d", i_stage);
                string_t tmp = allocate_strcat("phi_old:", iter_string);

                hsize_t dims[2] = {n_domain_cells, n_sol_variables};

                get_hdf5_dataset_n_m(last_id, tmp, HDF5Double, phi_old_restart[i_stage], 2, dims);

                deallocate(tmp);
            }
        }
    }

    close_hdf5_group(last_id);

    close_hdf5_file(file_id);

    copy_n(phi_total_restart, phi_total, n_tot_variables * n_domain_cells);
    copy_n(phi_dt_restart, phi_dt, n_sol_variables * n_domain_cells);

    for (int i_stage = 0; i_stage < n_stages_restart; i_stage++)
    {
        copy_n(phi_old_restart[i_stage], phi_old[i_stage], n_sol_variables * n_domain_cells);
    }

    deallocate(phi_total_restart);
    deallocate(phi_dt_restart);
    deallocate(phi_old_restart);
}