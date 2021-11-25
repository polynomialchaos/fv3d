/*******************************************************************************
 * @file limiter.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "fv_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"



double_limiter_fp_t limiter_function_pointer = NULL;

string_t limiter_name = NULL;

void limiter_initialize();
void limiter_finalize();

double limiter_none(int i_cell, int i_var, double slope);
double limiter_barth_jespersenn(int i_cell, int i_var, double slope);

void limiter_define()
{
    REGISTER_INITIALIZE_ROUTINE(limiter_initialize);
    REGISTER_FINALIZE_ROUTINE(limiter_finalize);

    string_t tmp_opt[] = {"Barth-Jespersenn", "None"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    SET_PARAMETER("FV/Limiter/limiter", StringParameter, &tmp,
                  "The limiter method", &tmp_opt, tmp_opt_n);
}

void limiter_initialize()
{
    GET_PARAMETER("FV/Limiter/limiter", StringParameter, &limiter_name);

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
        CHECK_EXPRESSION(0);
    }
}

void limiter_finalize()
{
    limiter_function_pointer = NULL;

    DEALLOCATE(limiter_name);
}

double limiter_none(int i_cell, int i_var, double slope)
{
#ifdef DEBUG
    UNUSED(i_cell);
    UNUSED(i_var);
    UNUSED(slope);
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

    for (int i = 0; i < cells->n_faces[i_cell]; ++i)
    {
        int *fc = &faces->cells[cf[i] * FACE_CELLS];

        if (fc[0] == i_cell)
        {
            phi_min = MIN(phi_min, phi_total[fc[1] * n_tot_variables + i_var]);
            phi_max = MAX(phi_max, phi_total[fc[1] * n_tot_variables + i_var]);
        }
        else
        {
            phi_min = MIN(phi_min, phi_total[fc[0] * n_tot_variables + i_var]);
            phi_max = MAX(phi_max, phi_total[fc[0] * n_tot_variables + i_var]);
        }
    }

    double tmp = 1.0;

    for (int i = 0; i < cells->n_faces[i_cell]; ++i)
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

        tmp = MIN(tmp, y);
    }

    return tmp;
}