//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "fv_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
double_limiter_fp_t limiter_function_pointer = NULL;

string_t limiter_name = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void limiter_initialize();
void limiter_finalize();

double limiter_none(int i_cell, int i_var, double slope);
double limiter_barth_jespersenn(int i_cell, int i_var, double slope);

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void limiter_define()
{
    register_initialize_routine(limiter_initialize);
    register_finalize_routine(limiter_finalize);

    string_t tmp_opt[] = {"Barth-Jespersenn", "None"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    set_parameter("FV/Limiter/limiter", ParameterString, &tmp,
                  "The limiter method", &tmp_opt, tmp_opt_n);
}

void limiter_initialize()
{
    get_parameter("FV/Limiter/limiter", ParameterString, &limiter_name);

    if (is_equal(limiter_name, "None"))
    {
        limiter_function_pointer = limiter_none;
    }
    else if (is_equal(limiter_name, "Barth-Jespersenn"))
    {
        limiter_function_pointer = limiter_barth_jespersenn;
    }
    else
    {
        check_error(0);
    }
}

void limiter_finalize()
{
    limiter_function_pointer = NULL;

    deallocate(limiter_name);
}

double limiter_none(int i_cell, int i_var, double slope)
{
#ifdef DEBUG
    u_unused(i_cell);
    u_unused(i_var);
    u_unused(slope);
#endif /* DEBUG */
    return 1.0;
}

double limiter_barth_jespersenn(int i_cell, int i_var, double slope)
{
    Cells_t *cells = global_mesh->cells;
    Faces_t *faces = global_mesh->faces;
    int max_cell_faces = cells->max_cell_faces;
    int n_tot_variables = all_variables->n_tot_variables;
    int *cf = &cells->faces[i_cell * max_cell_faces];

    double phi_min = phi_total[i_cell * n_tot_variables + i_var];
    double phi_max = phi_total[i_cell * n_tot_variables + i_var];

    for (int i = 0; i < cells->n_faces[i_cell]; i++)
    {
        int *fc = &faces->cells[cf[i] * FACE_CELLS];

        if (fc[0] == i_cell)
        {
            phi_min = u_min(phi_min, phi_total[fc[1] * n_tot_variables + i_var]);
            phi_max = u_max(phi_max, phi_total[fc[1] * n_tot_variables + i_var]);
        }
        else
        {
            phi_min = u_min(phi_min, phi_total[fc[0] * n_tot_variables + i_var]);
            phi_max = u_max(phi_max, phi_total[fc[0] * n_tot_variables + i_var]);
        }
    }

    double tmp = 1.0;

    for (int i = 0; i < cells->n_faces[i_cell]; i++)
    {
        double y = 1.0;
        if (slope > 0)
        {
            y = (phi_max - phi_total[i_cell * n_tot_variables + i_var]) / slope;
        }
        else if (slope < 0)
        {
            y = (phi_min - phi_total[i_cell * n_tot_variables + i_var]) / slope;
        }

        tmp = u_min(tmp, y);
    }

    return tmp;
}