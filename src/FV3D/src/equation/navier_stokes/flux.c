//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "navier_stokes_module.h"
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
string_t flux_scheme_name = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void flux_initialize();
void flux_finalize();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void flux_define()
{
    register_initialize_routine( flux_initialize );
    register_finalize_routine( flux_finalize );

    string_t tmp_opt[] = {"AUSM", "Rusanov"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    set_parameter( "Equation/Navier-Stokes/Flux/flux_scheme", ParameterString, &tmp,
        "The Riemann solver", &tmp_opt, tmp_opt_n );
}

void flux_initialize()
{
    get_parameter( "Equation/Navier-Stokes/Flux/flux_scheme", ParameterString, &flux_scheme_name );

        // select case( to_lower( strip( flux_scheme ) ) )
        //     case( 'rusanov' )
        //         flux_scheme_routine => riemann_rusanonv
        //     case( 'ausm' )
        //         flux_scheme_routine => riemann_ausm
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown Riemann solver selected (' // set_string( flux_scheme ) // ')!' )
        // end select
}

void flux_finalize()
{
    deallocate( flux_scheme_name );
}

void calc_flux()
{
        // do ii = 1, n_internal_faces
        //     i   = internal_faces(ii)
        //     fc  = faces(i)%cells

        //     ! rotate master and slave to face local coordinates
        //     phi_master_r[ic_rho]    = phi_total_left(IC_RHO,i)
        //     phi_master_r[ic_rho_u]  = dot_product( phi_total_left(IC_RHO_UVW,i), faces(i)%n )
        //     phi_master_r[ic_rho_v]  = dot_product( phi_total_left(IC_RHO_UVW,i), faces(i)%t1 )
        //     phi_master_r[ic_rho_w]  = dot_product( phi_total_left(IC_RHO_UVW,i), faces(i)%t2 )
        //     phi_master_r[ic_rho_e]  = phi_total_left(IC_RHO_E,i)
        //     phi_master_r            = con_to_prim( phi_master_r )

        //     phi_slave_r[ic_rho]     = phi_total_right(IC_RHO,i)
        //     phi_slave_r[ic_rho_u]   = dot_product( phi_total_right(IC_RHO_UVW,i), faces(i)%n )
        //     phi_slave_r[ic_rho_v]   = dot_product( phi_total_right(IC_RHO_UVW,i), faces(i)%t1 )
        //     phi_slave_r[ic_rho_w]   = dot_product( phi_total_right(IC_RHO_UVW,i), faces(i)%t2 )
        //     phi_slave_r[ic_rho_e]   = phi_total_right(IC_RHO_E,i)
        //     phi_slave_r             = con_to_prim( phi_slave_r )

        //     ! calculate convective flux
        //     call flux_scheme_routine( phi_master_r, phi_slave_r, flux_c )

        //     ! rotate convective flux back to global coordiantes and add to flux
        //     flux(1,i)   = flux_c(1)
        //     flux(2:4,i) = flux_c(2) * faces(i)%n + flux_c(3) * faces(i)%t1 + flux_c(4) * faces(i)%t2
        //     flux(5,i)   = flux_c(5)

        //     ! calculate diffusive flux
        //     call viscous_flux( phi_total_left(:,i), grad_phi_total_x(:,fc(1)), grad_phi_total_y(:,fc(1)), &
        //         grad_phi_total_z(:,fc(1)), phi_total_right(:,i), grad_phi_total_x(:,fc(2)), &
        //             grad_phi_total_y(:,fc(2)), grad_phi_total_z(:,fc(2)), flux_d_x, flux_d_y, flux_d_z )

        //     ! add diffusive flux to flux
        //     flux(:,i)   = flux(:,i) + flux_d_x * faces(i)%n(1) + flux_d_y * faces(i)%n(2) + flux_d_z * faces(i)%n(3)
        // end do

        // do ii = 1, n_flux_faces
        //     i   = flux_faces(ii)
        //     fc  = faces(i)%cells
        //     ir  = boundaries(faces(i)%boundary)%id

        //     ! rotate master and slave to face local coordinates
        //     phi_master_r[ic_rho]    = phi_total_left(IC_RHO,i)
        //     phi_master_r[ic_rho_u]  = dot_product( phi_total_left(IC_RHO_UVW,i), faces(i)%n )
        //     phi_master_r[ic_rho_v]  = dot_product( phi_total_left(IC_RHO_UVW,i), faces(i)%t1 )
        //     phi_master_r[ic_rho_w]  = dot_product( phi_total_left(IC_RHO_UVW,i), faces(i)%t2 )
        //     phi_master_r[ic_rho_e]  = phi_total_left(IC_RHO_E,i)
        //     phi_master_r            = con_to_prim( phi_master_r )

        //     phi_slave_r[ic_rho]     = phi_total_right(IC_RHO,i)
        //     phi_slave_r[ic_rho_u]   = dot_product( phi_total_right(IC_RHO_UVW,i), faces(i)%n )
        //     phi_slave_r[ic_rho_v]   = dot_product( phi_total_right(IC_RHO_UVW,i), faces(i)%t1 )
        //     phi_slave_r[ic_rho_w]   = dot_product( phi_total_right(IC_RHO_UVW,i), faces(i)%t2 )
        //     phi_slave_r[ic_rho_e]   = phi_total_right(IC_RHO_E,i)
        //     phi_slave_r             = con_to_prim( phi_slave_r )

        //     ! calculate convective flux
        //     select case ( regions(ir)%type )
        //         case default
        //             call flux_scheme_routine( phi_master_r, phi_slave_r, flux_c )
        //     end select

        //     ! rotate convective flux back to global coordiantes and add to flux
        //     flux(1,i)   = flux_c(1)
        //     flux(2:4,i) = flux_c(2) * faces(i)%n + flux_c(3) * faces(i)%t1 + flux_c(4) * faces(i)%t2
        //     flux(5,i)   = flux_c(5)

        //     ! calculate diffusive flux
        //     select case ( regions(ir)%type )
        //         case default
        //             call viscous_flux( phi_total_left(:,i), grad_phi_total_x(:,fc(1)), grad_phi_total_y(:,fc(1)), &
        //                 grad_phi_total_z(:,fc(1)), phi_total_right(:,i), grad_phi_total_x(:,fc(2)), &
        //                     grad_phi_total_y(:,fc(2)), grad_phi_total_z(:,fc(2)), flux_d_x, flux_d_y, flux_d_z )
        //     end select

        //     ! add diffusive flux to flux
        //     flux(:,i) = flux(:,i) + flux_d_x * faces(i)%n(1) + flux_d_y * faces(i)%n(2) + flux_d_z * faces(i)%n(3)
        // end do

}

void riemann_rusanonv( double *phi_l, double *phi_r, double *f )
{
    int n_sol_variables = all_variables->n_sol_variables;

    double c_l      = sqrt( kappa * phi_l[ip_p] / phi_l[ic_rho] );
    double c_r      = sqrt( kappa * phi_r[ip_p] / phi_r[ic_rho] );
    double eigval   = u_max( u_abs( phi_l[ip_u] ) + c_l, u_abs( phi_r[ip_u] ) + c_r );

    double f_l[n_sol_variables];
    double f_r[n_sol_variables];
    eval_euler_flux_1d( phi_l, f_l );
    eval_euler_flux_1d( phi_r, f_r );

    for ( int j = 0; j < n_sol_variables; j++ )
        f[j] = 0.5 * (f_l[j] + f_r[j]) - 0.5 * eigval * (phi_r[j] - phi_l[j]);
}

void riemann_ausm( double *phi_l, double *phi_r, double *f )
{
        // c_l = sqrt( kappa * phi_l[ip_p] / phi_l[ic_rho])
        // c_r = sqrt( kappa * phi_r[ip_p] / phi_r[ic_rho])

        // M_l = phi_l[ip_u] / c_l
        // M_r = phi_r[ip_u] / c_r

        // H_l = (phi_l[ic_rho_e] + phi_l[ip_p]) / phi_l[ic_rho]
        // H_r = (phi_r[ic_rho_e] + phi_r[ip_p]) / phi_r[ic_rho]

        // ! Positive M and p in the LEFT cell.
        // if( M_l .le. -1.0 ) then
        //     M_p = 0.0
        //     P_p = 0.0
        // else if( M_l .lt. 1.0 ) then
        //     M_p = 0.25 * (M_l + 1.0) * (M_l + 1.0)
        //     P_p = 0.25 * phi_l[ip_p] * (1.0 + M_l) * (1.0 + M_l) * (2.0 - M_l) ! or use P_p = half*(1.0+M_l)*phi_l[ip_v]
        // else
        //     M_p = M_l
        //     P_p = phi_l[ip_p]
        // end if

        // ! Negative M and p in the RIGHT cell.
        // if( M_r .le. -1.0 ) then
        //     M_m = M_r
        //     P_m = phi_r[ip_p]
        // else if( M_r .lt. 1.0 ) then
        //     M_m = -0.25 * (M_r - 1.0) * (M_r - 1.0)
        //     P_m =  0.25 * phi_r[ip_p] * (1.0 - M_r) * (1.0 - M_r) * (2.0 + M_r) ! or use P_m = half*(1.0-M_r)*phi_r[ip_v]
        // else
        //     M_m = 0.0
        //     P_m = 0.0
        // end if

        // ! Positive Part of Flux evaluated in the left cell.
        // f_p(1)  = u_max( 0.0, M_p + M_m ) * c_l * phi_l[ic_rho]
        // f_p(2)  = u_max( 0.0, M_p + M_m ) * c_l * phi_l[ic_rho] * phi_l[ip_u] + P_p
        // f_p(3)  = u_max( 0.0, M_p + M_m ) * c_l * phi_l[ic_rho] * phi_l[ip_v]
        // f_p(4)  = u_max( 0.0, M_p + M_m ) * c_l * phi_l[ic_rho] * phi_l[ip_w]
        // f_p(5)  = u_max( 0.0, M_p + M_m ) * c_l * phi_l[ic_rho] * H_l

        // ! Negative Part of Flux evaluated in the right cell.
        // f_m(1)  = u_min( 0.0, M_p + M_m ) * c_r * phi_r[ic_rho]
        // f_m(2)  = u_min( 0.0, M_p + M_m ) * c_r * phi_r[ic_rho] * phi_r[ip_u] + P_m
        // f_m(3)  = u_min( 0.0, M_p + M_m ) * c_r * phi_r[ic_rho] * phi_r[ip_v]
        // f_m(4)  = u_min( 0.0, M_p + M_m ) * c_r * phi_r[ic_rho] * phi_r[ip_w]
        // f_m(5)  = u_min( 0.0, M_p + M_m ) * c_r * phi_r[ic_rho] * H_r

        // f = f_p + f_m
}

void viscous_flux( double *phi_l, double *grad_phi_x_l, double *grad_phi_y_l, double *grad_phi_z_l,
        double *phi_r, double *grad_phi_x_r, double *grad_phi_y_r, double *grad_phi_z_r, double *f, double *g, double *h )
{
    // phi_l, grad_phi_x_l, grad_phi_y_l, grad_phi_z_l, &
    //     phi_r, grad_phi_x_r, grad_phi_y_r, grad_phi_z_r, f, g, h
        // call eval_viscous_flux_1d( phi_l, grad_phi_x_l, grad_phi_y_l, grad_phi_z_l, f_l, g_l, h_l )
        // call eval_viscous_flux_1d( phi_r, grad_phi_x_r, grad_phi_y_r, grad_phi_z_r, f_r, g_r, h_r )

        // f   = 0.5 * (f_l + f_r)
        // g   = 0.5 * (g_l + g_r)
        // h   = 0.5 * (h_l + h_r)
}

void eval_euler_flux_1d( double *phi, double *f )
{

        // f(1)    = phi[ic_rho_u]                             ! rho * u
        // f(2)    = phi[ic_rho_u] * phi[ip_u] + phi[ip_p]     ! rho * u * u + p
        // f(3)    = phi[ic_rho_u] * phi[ip_v]                 ! rho * u * v
        // f(4)    = phi[ic_rho_u] * phi[ip_w]                 ! rho * u * w
        // f(5)    = (phi[ic_rho_e] + phi[ip_p]) * phi[ip_u]   ! (rho * e + p) * u

}

void eval_viscous_flux_1d( double *phi, double *grad_phi_x, double *grad_phi_y, double *grad_phi_z,
        double *f, double *g, double *h )
{

        // tau_xx = mu_mix * ( s_43 * grad_phi_x(2) - &
        //     s_23 * grad_phi_y(3) - s_23 * grad_phi_z(4))        !  4/3 * mu * u_x - 2/3 * mu * v_y - 2/3 * mu * w_z
        // tau_yy = mu_mix * (-s_23 * grad_phi_x(2) + &
        //     s_43 * grad_phi_y(3) - s_23 * grad_phi_z(4))        ! -2/3 * mu * u_x + 4/3 * mu * v_y - 2/3 * mu * w_z
        // tau_zz = mu_mix * (-s_23 * grad_phi_x(2) - &
        //     s_23 * grad_phi_y(3) + s_43 * grad_phi_z(4))        ! -2/3 * mu * u_x - 2/3 * mu * v_y + 4/3 * mu * w_z

        // tau_xy = mu_mix * (grad_phi_y(2) + grad_phi_x(3))       ! mu * (u_y + v_x)
        // tau_xz = mu_mix * (grad_phi_z(2) + grad_phi_x(4))       ! mu * (u_z + w_x)
        // tau_yz = mu_mix * (grad_phi_z(3) + grad_phi_y(4))       ! mu * (y_z + w_y)

        // f(1)    =  0.0
        // f(2)    = -tau_xx                                       ! -4/3 * mu * u_x + 2/3 * mu * (v_y + w_z)
        // f(3)    = -tau_xy                                       ! -mu * (u_y + v_x)
        // f(4)    = -tau_xz                                       ! -mu * (u_z + w_x)
        // f(5)    = -tau_xx * phi[ip_u] - tau_xy * phi[ip_v] - &
        //     tau_xz * phi[ip_w] - lambda * grad_phi_x[ip_T]      ! -(tau_xx * phi + tau_xy * v + tau_xz * w - q_x) q_x=-lambda * T_x

        // g(1)    =  0.0
        // g(2)    = -tau_xy                                       ! -mu * (u_y + v_x)
        // g(3)    = -tau_yy                                       ! -4/3 * mu * v_y + 2/3 * mu * (u_x + w_z)
        // g(4)    = -tau_yz                                       ! -mu * (y_z + w_y)
        // g(5)    = -tau_xy * phi[ip_u] - tau_yy * phi[ip_v] - &
        //     tau_yz * phi[ip_w] - lambda * grad_phi_y[ip_T]      ! -(tau_yx * phi + tau_yy * v + tau_yz * w - q_y) q_y=-lambda * T_y

        // h(1)    =  0.0
        // h(2)    = -tau_xz                                       ! -mu * (u_z + w_x)
        // h(3)    = -tau_yz                                       ! -mu * (y_z + w_y)
        // h(4)    = -tau_zz                                       ! -4/3 * mu * w_z + 2/3 * mu * (u_x + v_y)
        // h(5)    = -tau_xz * phi[ip_u] - tau_yz * phi[ip_v] - &
        //     tau_zz * phi[ip_w] - lambda * grad_phi_z[ip_T]      ! -(tau_zx * phi + tau_zy * v + tau_zz * w - q_z) q_z=-lambda * T_z

}