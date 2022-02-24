/*******************************************************************************
 * @file boundary.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "navier_stokes_module.h"
#include "mesh/mesh_module.h"
#include "fv/fv_module.h"

string_t boundary_type_strings[BoundaryTypeMax] = {
    "flow",
    "inflow",
    "outflow",
    "wall-adiabatic",
    "wall-isothermal",
    "wall-slip",
    "symmetry",
    "state",
    "function"};

/*******************************************************************************
 * @brief Define boundary
 ******************************************************************************/
void boundary_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(boundary_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(boundary_finalize);
}

/*******************************************************************************
 * @brief Finalize boundary
 ******************************************************************************/
void boundary_finalize()
{
    if (solver_mesh != NULL)
    {
        Regions_t *regions = solver_mesh->regions;

        if (regions != NULL)
        {
            BM_DEALLOCATE(regions->type);
            BM_DEALLOCATE(regions->function_id);
            BM_DEALLOCATE(regions->phi_total);
        }
    }
}

/*******************************************************************************
 * @brief Initialize boundary
 ******************************************************************************/
void boundary_initialize()
{
    Regions_t *regions = solver_mesh->regions;
    int n_regions = regions->n_regions;
    int n_tot_variables = solver_variables->n_tot_variables;

    regions->type = BM_ALLOCATE(sizeof(int) * n_regions);
    regions->function_id = BM_ALLOCATE(sizeof(int) * n_regions);
    regions->phi_total =
        BM_ALLOCATE(sizeof(double) * n_tot_variables * n_regions);

    for (int i = 0; i < n_regions; ++i)
    {
        string_t path = allocate_strcat(
            "Equation/Navier-Stokes/Boundary/", regions->name[i]);
        BM_CHECK_EXPRESSION(BM_PARAMETER_EXISTS(path) == BC_TRUE);

        string_t type = allocate_strcat(path, "/type");
        string_t tmp;
        BM_GET_PARAMETER(type, StringParameter, &tmp);

        if (is_equal(tmp, boundary_type_strings[BoundaryFlow]))
        {
            BM_CHECK_EXPRESSION(regions->flow_region == i);
            regions->type[i] = BoundaryFlow;
            parse_primitive_state(
                path, &regions->phi_total[i * n_tot_variables]);
            prim_to_con(&regions->phi_total[i * n_tot_variables]);
        }
        else if (is_equal(tmp, boundary_type_strings[BoundaryInflow]))
        {
            regions->type[i] = BoundaryInflow;
            parse_primitive_state(
                path, &regions->phi_total[i * n_tot_variables]);
            prim_to_con(&regions->phi_total[i * n_tot_variables]);
        }
        else if (is_equal(tmp, boundary_type_strings[BoundaryOutflow]))
        {
            regions->type[i] = BoundaryOutflow;
        }
        else if (is_equal(tmp, boundary_type_strings[BoundaryAdiabaticWall]))
        {
            regions->type[i] = BoundaryAdiabaticWall;
        }
        else if (is_equal(tmp, boundary_type_strings[BoundaryIsothermalWall]))
        {
            regions->type[i] = BoundaryIsothermalWall;
            string_t tmp_T = allocate_strcat(path, "/T");
            BM_GET_PARAMETER(tmp_T, NumberParameter,
                             &regions->phi_total[i * n_tot_variables + ip_T]);
            BM_DEALLOCATE(tmp_T);
        }
        else if (is_equal(tmp, boundary_type_strings[BoundarySlipWall]))
        {
            regions->type[i] = BoundarySlipWall;
        }
        else if (is_equal(tmp, boundary_type_strings[BoundarySymmetry]))
        {
            regions->type[i] = BoundarySymmetry;
        }
        else if (is_equal(tmp, boundary_type_strings[BoundaryState]))
        {
            regions->type[i] = BoundaryState;
            parse_primitive_state(
                path, &regions->phi_total[i * n_tot_variables]);
            prim_to_con(&regions->phi_total[i * n_tot_variables]);
        }
        else if (is_equal(tmp, boundary_type_strings[BoundaryFunction]))
        {
            regions->type[i] = BoundaryFunction;
            string_t tmp_f = allocate_strcat(path, "/function_id");
            BM_GET_PARAMETER(tmp_f, DigitParameter, &regions->function_id[i]);
            BM_CHECK_EXPRESSION(regions->function_id[i] >= BoundaryTypeMax);
            BM_DEALLOCATE(tmp_f);
        }
        else
        {
            BM_CHECK_EXPRESSION(0);
        }

        BM_DEALLOCATE(tmp);
        BM_DEALLOCATE(type);
        BM_DEALLOCATE(path);
    }
}

/*******************************************************************************
 * @brief Parse primitive state parameter
 ******************************************************************************/
void parse_primitive_state(cstring_t prefix, double *phi)
{
    /* string_t tmp_rho    = allocate_strcat( prefix, "/rho" ); */
    string_t tmp_u = allocate_strcat(prefix, "/u");
    string_t tmp_v = allocate_strcat(prefix, "/v");
    string_t tmp_w = allocate_strcat(prefix, "/w");
    string_t tmp_p = allocate_strcat(prefix, "/p");
    string_t tmp_T = allocate_strcat(prefix, "/T");

    /* int has_rho = BM_PARAMETER_EXISTS( tmp_rho ); */
    int has_u = BM_PARAMETER_EXISTS(tmp_u);
    int has_v = BM_PARAMETER_EXISTS(tmp_v);
    int has_w = BM_PARAMETER_EXISTS(tmp_w);
    int has_p = BM_PARAMETER_EXISTS(tmp_p);
    int has_T = BM_PARAMETER_EXISTS(tmp_T);

    set_value_n(0.0, solver_variables->n_tot_variables, phi);

    /* check the provided data and fill the arrays */
    if (has_u && has_v && has_w && has_p && has_T)
    {
        BM_GET_PARAMETER(tmp_u, NumberParameter, &phi[ip_u]);
        BM_GET_PARAMETER(tmp_v, NumberParameter, &phi[ip_v]);
        BM_GET_PARAMETER(tmp_w, NumberParameter, &phi[ip_w]);
        BM_GET_PARAMETER(tmp_p, NumberParameter, &phi[ip_p]);
        BM_GET_PARAMETER(tmp_T, NumberParameter, &phi[ip_T]);
        phi[ic_rho] = calc_ig_rho(phi[ip_p], phi[ip_T], R_mix);
    }
    else
    {
        BM_CHECK_EXPRESSION(0);
    }

    /* BM_DEALLOCATE( tmp_rho ); */
    BM_DEALLOCATE(tmp_u);
    BM_DEALLOCATE(tmp_v);
    BM_DEALLOCATE(tmp_w);
    BM_DEALLOCATE(tmp_p);
    BM_DEALLOCATE(tmp_T);
}

/*******************************************************************************
 * @brief Update boundaries
 ******************************************************************************/
void update_boundaries(double t)
{
    Cells_t *cells = solver_mesh->cells;
    Boundaries_t *boundaries = solver_mesh->boundaries;
    Faces_t *faces = solver_mesh->faces;
    Regions_t *regions = solver_mesh->regions;
    int n_local_cells = cells->n_local_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_tot_variables = solver_variables->n_tot_variables;

    double tmp_u;
    double tmp_v;
    double tmp_w;

    for (int i = 0; i < n_boundaries; ++i)
    {
        int bf = boundaries->face[i];
        int bc = faces->cells[bf * FACE_CELLS];
        int id = boundaries->id[i];
        double *n = &faces->n[bf * DIM];
        double *t1 = &faces->t1[bf * DIM];
        double *t2 = &faces->t2[bf * DIM];

        double *phi_total_i =
            &solver_data->phi_total[(n_local_cells + i) * n_tot_variables];
        copy_n(&solver_data->phi_total[bc * n_tot_variables],
               n_tot_variables, phi_total_i);

        switch (regions->type[id])
        {
        case BoundaryInflow:
        case BoundaryState:
            copy_n(&regions->phi_total[id * n_tot_variables],
                   n_tot_variables, phi_total_i);
            break;
        case BoundaryOutflow:
            break;
        case BoundaryAdiabaticWall:
            /* rotate into local coordinate system */
            phi_total_i[ip_u] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], n, DIM);
            phi_total_i[ip_v] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], t1, DIM);
            phi_total_i[ip_w] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], t2, DIM);

            /* apply boundary specific behaviour */
            phi_total_i[ip_p] = calc_riemann_p(phi_total_i);
            set_value_n(0.0, DIM, &phi_total_i[ip_u]);
            phi_total_i[ic_rho] =
                calc_ig_rho(phi_total_i[ip_p], phi_total_i[ip_T], R_mix);

            /* rotate back to global coordinate system */
            tmp_u = phi_total_i[ip_u] * n[0] +
                    phi_total_i[ip_v] * t1[0] + phi_total_i[ip_w] * t2[0];
            tmp_v = phi_total_i[ip_u] * n[1] +
                    phi_total_i[ip_v] * t1[1] + phi_total_i[ip_w] * t2[1];
            tmp_w = phi_total_i[ip_u] * n[2] +
                    phi_total_i[ip_v] * t1[2] + phi_total_i[ip_w] * t2[2];

            phi_total_i[ip_u] = tmp_u;
            phi_total_i[ip_v] = tmp_v;
            phi_total_i[ip_w] = tmp_w;
            break;
        case BoundaryIsothermalWall:
            /* rotate into local coordinate system */
            phi_total_i[ip_u] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], n, DIM);
            phi_total_i[ip_v] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], t1, DIM);
            phi_total_i[ip_w] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], t2, DIM);

            /* apply boundary specific behaviour */
            phi_total_i[ip_p] = calc_riemann_p(phi_total_i);
            set_value_n(0.0, DIM, &phi_total_i[ip_u]);
            phi_total_i[ip_T] = regions->phi_total[id * n_tot_variables + ip_T];
            phi_total_i[ic_rho] =
                calc_ig_rho(phi_total_i[ip_p], phi_total_i[ip_T], R_mix);

            /* rotate back to global coordinate system */
            tmp_u = phi_total_i[ip_u] * n[0] +
                    phi_total_i[ip_v] * t1[0] + phi_total_i[ip_w] * t2[0];
            tmp_v = phi_total_i[ip_u] * n[1] +
                    phi_total_i[ip_v] * t1[1] + phi_total_i[ip_w] * t2[1];
            tmp_w = phi_total_i[ip_u] * n[2] +
                    phi_total_i[ip_v] * t1[2] + phi_total_i[ip_w] * t2[2];

            phi_total_i[ip_u] = tmp_u;
            phi_total_i[ip_v] = tmp_v;
            phi_total_i[ip_w] = tmp_w;
            break;
        case BoundarySlipWall:
        case BoundarySymmetry:
            /* rotate into local coordinate system */
            phi_total_i[ip_u] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], n, DIM);
            phi_total_i[ip_v] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], t1, DIM);
            phi_total_i[ip_w] = dot_n(
                &solver_data->phi_total[bc * n_tot_variables + ip_u], t2, DIM);

            /* apply boundary specific behaviour */
            phi_total_i[ip_p] = calc_riemann_p(phi_total_i);
            phi_total_i[ip_u] = 0.0;
            phi_total_i[ip_T] =
                calc_ig_T(phi_total_i[ip_p], phi_total_i[ic_rho], R_mix);

            /* rotate back to global coordinate system */
            tmp_u = phi_total_i[ip_u] * n[0] +
                    phi_total_i[ip_v] * t1[0] + phi_total_i[ip_w] * t2[0];
            tmp_v = phi_total_i[ip_u] * n[1] +
                    phi_total_i[ip_v] * t1[1] + phi_total_i[ip_w] * t2[1];
            tmp_w = phi_total_i[ip_u] * n[2] +
                    phi_total_i[ip_v] * t1[2] + phi_total_i[ip_w] * t2[2];

            phi_total_i[ip_u] = tmp_u;
            phi_total_i[ip_v] = tmp_v;
            phi_total_i[ip_w] = tmp_w;
            break;
        case BoundaryFunction:
            calc_ns_exact_function(
                regions->function_id[id], t, &faces->x[bf * DIM], phi_total_i);
            con_to_prim(phi_total_i);
            break;
        default:
            BM_CHECK_EXPRESSION(0);
            break;
        }

        prim_to_con(phi_total_i);
    }
}

/*******************************************************************************
 * @brief Update gradient boundaries
 ******************************************************************************/
void update_gradients_boundaries(double time)
{
#ifdef DEBUG
    BM_UNUSED(time);
#endif /* DEBUG */
    Cells_t *cells = solver_mesh->cells;
    Boundaries_t *boundaries = solver_mesh->boundaries;
    Faces_t *faces = solver_mesh->faces;
    Regions_t *regions = solver_mesh->regions;
    int n_local_cells = cells->n_local_cells;
    int n_boundaries = boundaries->n_boundaries;
    int n_tot_variables = solver_variables->n_tot_variables;

    double tmp_x[n_tot_variables];
    double tmp_y[n_tot_variables];
    double tmp_z[n_tot_variables];

    for (int i = 0; i < n_boundaries; ++i)
    {
        int bf = boundaries->face[i];
        int bc = faces->cells[bf * FACE_CELLS];
        int id = boundaries->id[i];
        double *n = &faces->n[bf * DIM];
        double *t1 = &faces->t1[bf * DIM];
        double *t2 = &faces->t2[bf * DIM];

        double *grad_phi_total_x_i =
            &solver_data->grad_phi_total_x
                 [(n_local_cells + i) * n_tot_variables];
        double *grad_phi_total_y_i =
            &solver_data->grad_phi_total_y
                 [(n_local_cells + i) * n_tot_variables];
        double *grad_phi_total_z_i =
            &solver_data->grad_phi_total_z
                 [(n_local_cells + i) * n_tot_variables];

        copy_n(&solver_data->grad_phi_total_x[bc * n_tot_variables],
               n_tot_variables, grad_phi_total_x_i);
        copy_n(&solver_data->grad_phi_total_y[bc * n_tot_variables],
               n_tot_variables, grad_phi_total_y_i);
        copy_n(&solver_data->grad_phi_total_z[bc * n_tot_variables],
               n_tot_variables, grad_phi_total_z_i);

        switch (regions->type[id])
        {
        case BoundaryInflow:
        case BoundaryState:
        case BoundaryOutflow:
        case BoundaryAdiabaticWall:
        case BoundaryIsothermalWall:
        case BoundaryFunction:
            break;
        case BoundarySlipWall:
        case BoundarySymmetry:
            /* rotate neighbour cell gradient into local coordinates */
            for (int j = 0; j < n_tot_variables; ++j)
            {
                tmp_x[j] = grad_phi_total_x_i[j] * n[0] +
                           grad_phi_total_y_i[j] * n[1] +
                           grad_phi_total_z_i[j] * n[2];
                tmp_y[j] = grad_phi_total_x_i[j] * t1[0] +
                           grad_phi_total_y_i[j] * t1[1] +
                           grad_phi_total_z_i[j] * t1[2];
                tmp_z[j] = grad_phi_total_x_i[j] * t2[0] +
                           grad_phi_total_y_i[j] * t2[1] +
                           grad_phi_total_z_i[j] * t2[2];
            }

            /* apply boundary specific behaviour */
            set_value_n(0.0, n_tot_variables, tmp_x);

            /* rotate neighbour cell gradient back from local coordinates */
            for (int j = 0; j < n_tot_variables; ++j)
            {
                grad_phi_total_x_i[j] =
                    tmp_x[j] * n[0] + tmp_y[j] * t1[0] + tmp_z[j] * t2[0];
                grad_phi_total_y_i[j] =
                    tmp_x[j] * n[1] + tmp_y[j] * t1[1] + tmp_z[j] * t2[1];
                grad_phi_total_z_i[j] =
                    tmp_x[j] * n[2] + tmp_y[j] * t1[2] + tmp_z[j] * t2[2];
            }
            break;
        default:
            BM_CHECK_EXPRESSION(0);
            break;
        }
    }
}
