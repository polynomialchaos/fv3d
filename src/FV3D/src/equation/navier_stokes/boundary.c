//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "navier_stokes_module.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"
#include "fv/fv_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------
enum BoundaryType {
    BoundaryFlow, BoundaryInflow, BoundaryOutflow,
    BoundaryAdiabaticWall, BoundaryIsothermalWall, BoundarySlipWall,
    BoundarySymmetry, BoundaryState, BoundaryFunction,
    BoundaryTypeMax
};

string_t boundary_type_strings[BoundaryTypeMax] =
{
    "FLOW",
    "INFLOW",
    "OUTFLOW",
    "ADIABATIC_WALL",
    "ISOTHERMAL_WALL",
    "SLIP_WALL",
    "SYMMETRY",
    "STATE",
    "FUNCTION"
};

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void boundary_initialize();
void boundary_finalize();

void parse_primitive_state( const_string_t prefix, double *phi );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void boundary_define()
{
    register_initialize_routine( boundary_initialize );
    register_finalize_routine( boundary_finalize );
}

void boundary_initialize()
{
    Regions_t *regions  = global_mesh->regions;
    int n_regions       = regions->n_regions;
    int n_tot_variables = all_variables->n_tot_variables;

    regions->type           = allocate( sizeof( int ) * n_regions );
    regions->function_id    = allocate( sizeof( int ) * n_regions );
    regions->phi_total      = allocate( sizeof( double ) * n_tot_variables * n_regions );

    for ( int i = 0; i < n_regions; i++ )
    {
        string_t path = allocate_strcat( "Equation/Navier-Stokes/Boundary/", regions->name[i] );
        check_error( (parameter_exists( path ) == 1) );

        string_t type = allocate_strcat( path, "/type" );
        string_t tmp;
        get_parameter( type, ParameterString, &tmp );

        if (is_equal( tmp, boundary_type_strings[BoundaryFlow] ))
        {
            check_error( (regions->flow_region == i) );
            regions->type[i] = BoundaryFlow;
            parse_primitive_state( path, &regions->phi_total[i*n_tot_variables] );
            prim_to_con( &regions->phi_total[i*n_tot_variables] );
        }
        else if (is_equal( tmp, boundary_type_strings[BoundaryInflow] ))
        {
            regions->type[i] = BoundaryInflow;
            parse_primitive_state( path, &regions->phi_total[i*n_tot_variables] );
            prim_to_con( &regions->phi_total[i*n_tot_variables] );
        }
        else if (is_equal( tmp, boundary_type_strings[BoundaryOutflow] ))
        {
            regions->type[i] = BoundaryOutflow;
        }
        else if (is_equal( tmp, boundary_type_strings[BoundaryAdiabaticWall] ))
        {
            regions->type[i] = BoundaryAdiabaticWall;
        }
        else if (is_equal( tmp, boundary_type_strings[BoundaryIsothermalWall] ))
        {
            regions->type[i]    = BoundaryIsothermalWall;
            string_t tmp_T      = allocate_strcat( path, "/T" );
            get_parameter( tmp_T, ParameterNumber, &regions->phi_total[i*n_tot_variables+ip_T] );
            deallocate( tmp_T );
        }
        else if (is_equal( tmp, boundary_type_strings[BoundarySlipWall] ))
        {
            regions->type[i] = BoundarySlipWall;
        }
        else if (is_equal( tmp, boundary_type_strings[BoundarySymmetry] ))
        {
            regions->type[i] = BoundarySymmetry;
        }
        else if (is_equal( tmp, boundary_type_strings[BoundaryState] ))
        {
            regions->type[i] = BoundaryState;
            parse_primitive_state( path, &regions->phi_total[i*n_tot_variables] );
            prim_to_con( &regions->phi_total[i*n_tot_variables] );
        }
        else if (is_equal( tmp, boundary_type_strings[BoundaryFunction] ))
        {
            regions->type[i]    = BoundaryFunction;
            string_t tmp_f      = allocate_strcat( path, "/function_id" );
            get_parameter( tmp_f, ParameterDigit, &regions->function_id[i] );
            deallocate( tmp_f );
        }
        else
        {
            check_error( 0 );
        }

        deallocate( tmp );
        deallocate( type );
        deallocate( path );
    }
}

void boundary_finalize()
{
    Regions_t *regions = global_mesh->regions;

    deallocate( regions->type );
    deallocate( regions->function_id );
    deallocate( regions->phi_total );
}

void update_boundaries( double t )
{
    Cells_t *cells              = global_mesh->cells;
    Boundaries_t *boundaries    = global_mesh->boundaries;
    Faces_t *faces              = global_mesh->faces;
    Regions_t *regions          = global_mesh->regions;

    int n_local_cells           = cells->n_local_cells;
    int n_boundaries            = boundaries->n_boundaries;
    int n_tot_variables         = all_variables->n_tot_variables;

    double tmp_u;
    double tmp_v;
    double tmp_w;

    for ( int i = 0; i < n_boundaries; i++ )
    {
        int bf  = boundaries->face[i];
        int bc  = faces->cells[bf*FACE_CELLS];
        int id  = boundaries->id[i];

        double *phi_total_i = &phi_total[(n_local_cells + i) * n_tot_variables];
        copy_n( &phi_total[bc * n_tot_variables], phi_total_i, n_tot_variables );

        switch (regions->type[id])
        {
            case BoundaryInflow:
            case BoundaryState:
                copy_n( &regions->phi_total[id*n_tot_variables], phi_total_i, n_tot_variables );
                break;
            case BoundaryOutflow:
                break;
            case BoundaryAdiabaticWall:
                // rotate into local coordinate system
                phi_total_i[ip_u]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->n[bf*DIM], DIM );
                phi_total_i[ip_v]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->t1[bf*DIM], DIM );
                phi_total_i[ip_w]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->t2[bf*DIM], DIM );

                //apply boundary specific behaviour
                phi_total_i[ip_p]   = calc_riemann_p( phi_total_i );
                set_value_n( 0.0, &phi_total_i[ip_u], DIM );
                phi_total_i[ic_rho] = calc_ig_rho( phi_total_i[ip_p], phi_total_i[ip_T], R_mix );

                // rotate back to global coordinate system
                tmp_u   = phi_total_i[ip_u] * faces->n[bf*DIM] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM] + phi_total_i[ip_w] * faces->t2[bf*DIM];
                tmp_v   = phi_total_i[ip_u] * faces->n[bf*DIM+1] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM+1] + phi_total_i[ip_w] * faces->t2[bf*DIM+1];
                tmp_w   = phi_total_i[ip_u] * faces->n[bf*DIM+2] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM+2] + phi_total_i[ip_w] * faces->t2[bf*DIM+2];

                phi_total_i[ip_u]   = tmp_u;
                phi_total_i[ip_v]   = tmp_v;
                phi_total_i[ip_w]   = tmp_w;
                break;
            case BoundaryIsothermalWall:
                // rotate into local coordinate system
                phi_total_i[ip_u]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->n[bf*DIM], DIM );
                phi_total_i[ip_v]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->t1[bf*DIM], DIM );
                phi_total_i[ip_w]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->t2[bf*DIM], DIM );

                //apply boundary specific behaviour
                phi_total_i[ip_p]   = calc_riemann_p( phi_total_i );
                set_value_n( 0.0, &phi_total_i[ip_u], DIM );
                phi_total_i[ip_T]   = regions->phi_total[id*n_tot_variables+ip_T];
                phi_total_i[ic_rho] = calc_ig_rho( phi_total_i[ip_p], phi_total_i[ip_T], R_mix );

                // rotate back to global coordinate system
                tmp_u   = phi_total_i[ip_u] * faces->n[bf*DIM] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM] + phi_total_i[ip_w] * faces->t2[bf*DIM];
                tmp_v   = phi_total_i[ip_u] * faces->n[bf*DIM+1] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM+1] + phi_total_i[ip_w] * faces->t2[bf*DIM+1];
                tmp_w   = phi_total_i[ip_u] * faces->n[bf*DIM+2] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM+2] + phi_total_i[ip_w] * faces->t2[bf*DIM+2];

                phi_total_i[ip_u]   = tmp_u;
                phi_total_i[ip_v]   = tmp_v;
                phi_total_i[ip_w]   = tmp_w;
                break;
            case BoundarySlipWall:
            case BoundarySymmetry:
                // rotate into local coordinate system
                phi_total_i[ip_u]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->n[bf*DIM], DIM );
                phi_total_i[ip_v]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->t1[bf*DIM], DIM );
                phi_total_i[ip_w]   = dot_n( &phi_total[bc*n_tot_variables+ip_u], &faces->t2[bf*DIM], DIM );

                //apply boundary specific behaviour
                phi_total_i[ip_p]   = calc_riemann_p( phi_total_i );
                phi_total_i[ip_u]   = 0.0;
                phi_total_i[ip_T]   = calc_ig_T( phi_total_i[ip_p], phi_total_i[ic_rho], R_mix );

                // rotate back to global coordinate system
                tmp_u   = phi_total_i[ip_u] * faces->n[bf*DIM] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM] + phi_total_i[ip_w] * faces->t2[bf*DIM];
                tmp_v   = phi_total_i[ip_u] * faces->n[bf*DIM+1] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM+1] + phi_total_i[ip_w] * faces->t2[bf*DIM+1];
                tmp_w   = phi_total_i[ip_u] * faces->n[bf*DIM+2] +
                    phi_total_i[ip_v] * faces->t1[bf*DIM+2] + phi_total_i[ip_w] * faces->t2[bf*DIM+2];

                phi_total_i[ip_u]   = tmp_u;
                phi_total_i[ip_v]   = tmp_v;
                phi_total_i[ip_w]   = tmp_w;
                break;
            case BoundaryFunction:
                calc_exact_function_pointer( regions->function_id[id], t, &faces->x[bf*DIM], phi_total_i );
                con_to_prim( phi_total_i );
                break;
            default:
                check_error( 0 );
                break;
        }

        prim_to_con( phi_total_i );
    }
}

void update_gradients_boundaries( double t )
{
    Cells_t *cells              = global_mesh->cells;
    Boundaries_t *boundaries    = global_mesh->boundaries;
    Faces_t *faces              = global_mesh->faces;
    Regions_t *regions          = global_mesh->regions;

    int n_local_cells           = cells->n_local_cells;
    int n_boundaries            = boundaries->n_boundaries;
    int n_tot_variables         = all_variables->n_tot_variables;

    double tmp_x[n_tot_variables];
    double tmp_y[n_tot_variables];
    double tmp_z[n_tot_variables];

    for ( int i = 0; i < n_boundaries; i++ )
    {
        int bf  = boundaries->face[i];
        int bc  = faces->cells[bf*FACE_CELLS];
        int id  = boundaries->id[i];

        double *grad_phi_total_x_i = &grad_phi_total_x[(n_local_cells + i) * n_tot_variables];
        double *grad_phi_total_y_i = &grad_phi_total_y[(n_local_cells + i) * n_tot_variables];
        double *grad_phi_total_z_i = &grad_phi_total_z[(n_local_cells + i) * n_tot_variables];
        copy_n( &grad_phi_total_x[bc*n_tot_variables], grad_phi_total_x_i, n_tot_variables );
        copy_n( &grad_phi_total_y[bc*n_tot_variables], grad_phi_total_y_i, n_tot_variables );
        copy_n( &grad_phi_total_z[bc*n_tot_variables], grad_phi_total_z_i, n_tot_variables );

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
                // rotate neighbour cell gradient into local coordinates
                for ( int j = 0; j < n_tot_variables; j++ )
                {
                    tmp_x[j] = grad_phi_total_x_i[j] * faces->n[bf*DIM] +
                        grad_phi_total_y_i[j] * faces->n[bf*DIM+1] + grad_phi_total_z_i[j] * faces->n[bf*DIM+2];
                    tmp_y[j] = grad_phi_total_x_i[j] * faces->t1[bf*DIM] +
                        grad_phi_total_y_i[j] * faces->t1[bf*DIM+1] + grad_phi_total_z_i[j] * faces->t1[bf*DIM+2];
                    tmp_z[j] = grad_phi_total_x_i[j] * faces->t2[bf*DIM] +
                        grad_phi_total_y_i[j] * faces->t2[bf*DIM+1] + grad_phi_total_z_i[j] * faces->t2[bf*DIM+2];
                }

                //apply boundary specific behaviour
                for ( int j = 0; j < n_tot_variables; j++ )
                {
                    tmp_x[j] = 0.0;
                }

                // rotate neighbour cell gradient back from local coordinates
                for ( int j = 0; j < n_tot_variables; j++ )
                {
                    grad_phi_total_x_i[j] = tmp_x[j] * faces->n[bf*DIM] +
                        tmp_y[j] * faces->t1[bf*DIM] + tmp_z[j] * faces->t2[bf*DIM];
                    grad_phi_total_y_i[j] = tmp_x[j] * faces->n[bf*DIM+1] +
                        tmp_y[j] * faces->t1[bf*DIM+1] + tmp_z[j] * faces->t2[bf*DIM+1];
                    grad_phi_total_z_i[j] = tmp_x[j] * faces->n[bf*DIM+2] +
                        tmp_y[j] * faces->t1[bf*DIM+2] + tmp_z[j] * faces->t2[bf*DIM+2];
                }
                break;
            default:
                check_error( 0 );
                break;
        }
    }
}

void parse_primitive_state( const_string_t prefix, double *phi )
{
    string_t tmp_rho    = allocate_strcat( prefix, "/rho" );
    string_t tmp_u      = allocate_strcat( prefix, "/u" );
    string_t tmp_v      = allocate_strcat( prefix, "/v" );
    string_t tmp_w      = allocate_strcat( prefix, "/w" );
    string_t tmp_p      = allocate_strcat( prefix, "/p" );
    string_t tmp_T      = allocate_strcat( prefix, "/T" );

    int has_rho = parameter_exists( tmp_rho );
    int has_u   = parameter_exists( tmp_u );
    int has_v   = parameter_exists( tmp_v );
    int has_w   = parameter_exists( tmp_w );
    int has_p   = parameter_exists( tmp_p );
    int has_T   = parameter_exists( tmp_T );

    set_value_n( 0.0, phi, all_variables->n_tot_variables );

    // check the provided data and fill the arrays
    if (has_u && has_v && has_w && has_p && has_T)
    {
        get_parameter( tmp_u, ParameterNumber, &phi[ip_u] );
        get_parameter( tmp_v, ParameterNumber, &phi[ip_v] );
        get_parameter( tmp_w, ParameterNumber, &phi[ip_w] );
        get_parameter( tmp_p, ParameterNumber, &phi[ip_p] );
        get_parameter( tmp_T, ParameterNumber, &phi[ip_T] );
        phi[ic_rho] = calc_ig_rho( phi[ip_p], phi[ip_T], R_mix );
    }
    else
    {
        check_error( 0 );
    }

    deallocate( tmp_rho );
    deallocate( tmp_u );
    deallocate( tmp_v );
    deallocate( tmp_w );
    deallocate( tmp_p );
    deallocate( tmp_T );
}