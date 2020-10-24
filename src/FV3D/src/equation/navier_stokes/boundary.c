//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "navier_stokes_module.h"
#include "equation/equation_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------
enum BoundaryType {
    BoundaryFlow, BoundaryInflow, BoundaryOutflow, BoundaryAdiabaticWall,
    BoundaryIsothermalWall, BoundarySlipWall, BoundarySymmetry, BoundaryState, BoundaryFunction,
    BoundaryTypeMax
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
        // do i = 1, n_regions
        //     prefix = 'Equation/Navier-Stokes/Boundary/' // strip( regions(i)%name )
        //     if( .not. parameter_exists( prefix ) ) &
        //         call add_error( __LINE__, __FILE__, &
        //             'Missing boundary specification (' // set_string( prefix ) // ')!' )

        //     call get_parameter( prefix // '/type', bnd_type )

        //     allocate( regions(i)%phi(n_tot_variables) ); regions(i)%phi = 0.0

        //     select case( to_lower( strip( bnd_type ) ) )
        //         case ( bnd_type_names(BND_FLOW) )
        //             regions(i)%type     = BND_FLOW
        //             flow_region         = i
        //             regions(i)%phi      = prim_to_con( parse_primitive_state( prefix ) )
        //         case ( bnd_type_names(BND_INFLOW) )
        //             regions(i)%type     = BND_INFLOW
        //             regions(i)%phi      = prim_to_con( parse_primitive_state( prefix ) )
        //         case ( bnd_type_names(BND_OUTFLOW) )
        //             regions(i)%type     = BND_OUTFLOW
        //         case ( bnd_type_names(BND_AD_WALL) )
        //             regions(i)%type     = BND_AD_WALL
        //         case ( bnd_type_names(BND_IT_WALL) )
        //             regions(i)%type     = BND_IT_WALL
        //             call get_parameter( prefix // '/T', regions(i)%phi(IP_T) )
        //         case ( bnd_type_names(BND_SLIP_WALL) )
        //             regions(i)%type     = BND_SLIP_WALL
        //         case ( bnd_type_names(BND_SYMM) )
        //             regions(i)%type     = BND_SYMM
        //         case ( bnd_type_names(BND_STATE) )
        //             regions(i)%phi      = prim_to_con( parse_primitive_state( prefix ) )
        //         case ( bnd_type_names(BND_FUNC) )
        //             regions(i)%type = BND_FUNC
        //             call get_parameter( prefix // '/function_id', regions(i)%function_id )
        //         case default
        //             call add_error( __LINE__, __FILE__, &
        //                 'Unknown boundary type provided (' // set_string( bnd_type ) // ')!' )
        //     end select
        // end do
}

void boundary_finalize()
{
}

void update_boundaries( double t )
{
       //     do i = 1, n_boundaries
        //     j   = n_cells + n_partition_receives + i
        //     bf  = boundaries(i)%face
        //     fc  = faces(bf)%cells(1)
        //     ir  = boundaries(i)%id

        //     phi_total(:,j)  = phi_total(:,fc)

        //     select case ( regions(ir)%type )
        //         case ( BND_INFLOW )
        //             phi_total(:,j)      = regions(ir)%phi
        //         case ( BND_AD_WALL )
        //             ! rotate into local coordinate system
        //             phi_total(IP_U,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%n )
        //             phi_total(IP_V,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%t1 )
        //             phi_total(IP_W,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%t2 )

        //             ! apply boundary specific behaviour
        //             phi_total(IP_P,j)   = pressure_riemann( phi_total(:,j) )
        //             phi_total(IP_UVW,j) = 0.0
        //             phi_total(IC_RHO,j) = ig_density( phi_total(IP_P,j), phi_total(IP_T,j), R_mix )

        //             ! rotate back to global coordinate system
        //             phi_total(IP_UVW,j) = phi_total(IP_U,j) * faces(bf)%n + &
        //                 phi_total(IP_V,j) * faces(bf)%t1 + phi_total(IP_W,j) * faces(bf)%t2
        //         case ( BND_IT_WALL )
        //             ! rotate into local coordinate system
        //             phi_total(IP_U,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%n )
        //             phi_total(IP_V,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%t1 )
        //             phi_total(IP_W,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%t2 )

        //             ! apply boundary specific behaviour
        //             phi_total(IP_P,j)   = pressure_riemann( phi_total(:,j) )
        //             phi_total(IP_UVW,j) = 0.0
        //             phi_total(IP_T,j)   = regions(ir)%phi(IP_T)
        //             phi_total(IC_RHO,j) = ig_density( phi_total(IP_P,j), phi_total(IP_T,j), R_mix )

        //             ! rotate back to global coordinate system
        //             phi_total(IP_UVW,j) = phi_total(IP_U,j) * faces(bf)%n + &
        //                 phi_total(IP_V,j) * faces(bf)%t1 + phi_total(IP_W,j) * faces(bf)%t2
        //         case ( BND_SLIP_WALL, BND_SYMM )
        //             ! rotate into local coordinate system
        //             phi_total(IP_U,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%n )
        //             phi_total(IP_V,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%t1 )
        //             phi_total(IP_W,j)   = dot_product( phi_total(IP_UVW,j), faces(bf)%t2 )

        //             ! apply boundary specific behaviour
        //             phi_total(IP_P,j)   = pressure_riemann( phi_total(:,j) )
        //             phi_total(IP_U,j)   = 0.0
        //             phi_total(IP_T,j)   = ig_temperature( phi_total(IP_P,j), phi_total(IC_RHO,j), R_mix )

        //             ! rotate back to global coordinate system
        //             phi_total(IP_UVW,j) = phi_total(IP_U,j) * faces(bf)%n + &
        //                 phi_total(IP_V,j) * faces(bf)%t1 + phi_total(IP_W,j) * faces(bf)%t2
        //         case ( BND_STATE )
        //             phi_total(IP,j)     = regions(ir)%phi(IP)
        //         case ( BND_FUNC )
        //             call exact_func_routine( regions(ir)%function_id, t, faces(bf)%x, phi_total(:,j) )
        //             phi_total(:,j) = con_to_prim( phi_total(:,j) )
        //         case ( BND_OUTFLOW )
        //         case default
        //             call add_error( __LINE__, __FILE__, &
        //                 'Unsupported region type provided (' // set_string( regions(ir)%type ) // ')!' )
        //     end select

        //     phi_total(:,j) = prim_to_con( phi_total(:,j) )
        // end do
}

void update_gradients_boundaries( double t )
{
    //    do i = 1, n_boundaries
    //         j   = n_cells + n_partition_receives + i
    //         bf  = boundaries(i)%face
    //         bc  = faces(bf)%cells(1)

    //         ! rotate neighbour cell gradient into local coordinates
    //         grad_phi_total_x_r  = grad_phi_total_x(:,bc) * faces(bf)%n(1) + &
    //             grad_phi_total_y(:,bc) * faces(bf)%n(2) + grad_phi_total_z(:,bc) * faces(bf)%n(3)

    //         grad_phi_total_y_r  = grad_phi_total_x(:,bc) * faces(bf)%t1(1) + &
    //             grad_phi_total_y(:,bc) * faces(bf)%t1(2) + grad_phi_total_z(:,bc) * faces(bf)%t1(3)

    //         grad_phi_total_z_r  = grad_phi_total_x(:,bc) * faces(bf)%t2(1) + &
    //             grad_phi_total_y(:,bc) * faces(bf)%t2(2) + grad_phi_total_z(:,bc) * faces(bf)%t2(3)

    //         select case ( regions(boundaries(i)%id)%type )
    //             case ( BND_SLIP_WALL, BND_SYMM )
    //                 grad_phi_total_x_r(:) = 0.0
    //         end select

    //         ! rotate neighbour cell gradient back from local coordinates
    //         grad_phi_total_x(:,j) = grad_phi_total_x_r * faces(bf)%n(1) + &
    //             grad_phi_total_y_r * faces(bf)%t1(1) + grad_phi_total_z_r * faces(bf)%t2(1)

    //         grad_phi_total_y(:,j) = grad_phi_total_x_r * faces(bf)%n(2) + &
    //             grad_phi_total_y_r * faces(bf)%t1(2) + grad_phi_total_z_r * faces(bf)%t2(2)

    //         grad_phi_total_z(:,j) = grad_phi_total_x_r * faces(bf)%n(3) + &
    //             grad_phi_total_y_r * faces(bf)%t1(3) + grad_phi_total_z_r * faces(bf)%t2(3)
    //     end do
}

void parse_primitive_state( const_string_t prefix, double *phi_total_i )
{
    string_t tmp_rho    = allocate_strcat( prefix, "rho" );
    string_t tmp_u      = allocate_strcat( prefix, "u" );
    string_t tmp_v      = allocate_strcat( prefix, "v" );
    string_t tmp_w      = allocate_strcat( prefix, "w" );
    string_t tmp_p      = allocate_strcat( prefix, "p" );
    string_t tmp_T      = allocate_strcat( prefix, "T" );

    int has_rho = parameter_exists( tmp_rho );
    int has_u   = parameter_exists( tmp_u );
    int has_v   = parameter_exists( tmp_v );
    int has_w   = parameter_exists( tmp_w );
    int has_p   = parameter_exists( tmp_p );
    int has_T   = parameter_exists( tmp_T );

    set_value_n( 0.0, phi_total_i, all_variables->n_tot_variables );

    // check the provided data and fill the arrays
    if (has_rho && has_u && has_v && has_w && has_p)
    {
        // get_parameter( tmp_rho, res[ic_rho] );
        // get_parameter( tmp_u, res[ip_u] );
        // get_parameter( tmp_v, res[ip_v] );
        // get_parameter( tmp_w, res[ip_w] );
        // get_parameter( tmp_p, res[ip_p] );
        // ip_T =
        // get_parameter( tmp_T, res[] );

        // call get_parameter( prefix // '/rho', res(IC_RHO) )
        // call get_parameter( prefix // '/v_x', res(IP_U) )
        // call get_parameter( prefix // '/v_y', res(IP_V) )
        // call get_parameter( prefix // '/v_z', res(IP_W) )
        // call get_parameter( prefix // '/p', res(IP_P) )
        // res(IP_T) = ig_temperature( res(IP_P), res(IC_RHO), R_mix )
    }
    else if (has_rho && has_u && has_v && has_w && has_T)
    {
        // call get_parameter( prefix // '/rho', res(IC_RHO) )
        // call get_parameter( prefix // '/v_x', res(IP_U) )
        // call get_parameter( prefix // '/v_y', res(IP_V) )
        // call get_parameter( prefix // '/v_z', res(IP_W) )
        // call get_parameter( prefix // '/T', res(IP_T) )
        // res(IP_P) = ig_pressure( res(IC_RHO), res(IP_T), R_mix )
    }
    else if (has_u && has_v && has_w && has_p && has_T)
    {
        // call get_parameter( prefix // '/v_x', res(IP_U) )
        // call get_parameter( prefix // '/v_y', res(IP_V) )
        // call get_parameter( prefix // '/v_z', res(IP_W) )
        // call get_parameter( prefix // '/p', res(IP_P) )
        // call get_parameter( prefix // '/T', res(IP_T) )
        // res(IC_RHO) = ig_density( res(IP_P), res(IP_T), R_mix )
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