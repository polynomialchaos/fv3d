
//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "implicit_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int implicit_active = 0;
int n_iter_inner    = 0;
int n_iter_lsoe     = 0;

string_t implicit_scheme_name   = NULL;
string_t method_name            = NULL;
int max_iter_inner              = 100;
string_t solver_name            = NULL;
string_t jacobian_type_name     = NULL;
double tolerance_lsoe           = 1e-12;
int max_iter_lsoe               = 100;
int max_krylov_dims             = 15;
int max_krylov_restarts         = 2;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void implicit_initialize();
void implicit_finalize();
void matrix_vector_numerical( double *x, double *b, int n, int m );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void implicit_define()
{
    register_initialize_routine( implicit_initialize );
    register_finalize_routine( implicit_finalize );

    string_t tmp_opt[] = {"BDF-2", "Euler"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    string_t tmp2_opt[] = {"Newton"};
    int tmp2_opt_n = sizeof( tmp2_opt ) / sizeof( string_t );
    string_t tmp2 = tmp2_opt[0];

    string_t tmp3_opt[] = {"Numerical"};
    int tmp3_opt_n = sizeof( tmp3_opt ) / sizeof( string_t );
    string_t tmp3 = tmp3_opt[0];

    string_t tmp4_opt[] = {"BiCGStab", "GMRes"};
    int tmp4_opt_n = sizeof( tmp4_opt ) / sizeof( string_t );
    string_t tmp4 = tmp4_opt[0];

    set_parameter( "TimeDisc/Implicit/scheme", ParameterString, &tmp, "The implicit timestep scheme", &tmp_opt, tmp_opt_n );
    set_parameter( "TimeDisc/Implicit/method", ParameterString, &tmp2,
        "The method to solve the non-linear system of equations", &tmp2_opt, tmp2_opt_n );
    set_parameter( "TimeDisc/Implicit/max_iter_inner", ParameterDigit, &max_iter_inner,
        "The maximum number of inner iterations", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/solver", ParameterString, &tmp4,
        "The linear system of equations solver", &tmp4_opt, tmp4_opt_n );
    set_parameter( "TimeDisc/Implicit/jacobian_type", ParameterString, &tmp3,
        "The type of jacobian generation", &tmp3_opt, tmp3_opt_n );
    set_parameter( "TimeDisc/Implicit/tolerance_lsoe", ParameterNumber, &tolerance_lsoe, "he linear solver tolerance", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/max_iter_lsoe", ParameterDigit, &max_iter_lsoe,
        "The linear solver maximum number of iterations", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/max_krylov_dims", ParameterDigit, &max_krylov_dims,
        "The maximum Krylov space dimension in GMRes solver", NULL, 0 );
    set_parameter( "TimeDisc/Implicit/max_krylov_restarts", ParameterDigit, &max_krylov_restarts,
        "The maximum restarts performed in GMRes solver", NULL, 0 );
}

void implicit_initialize()
{
    if (implicit_active == 0) return;

    get_parameter( "TimeDisc/Implicit/scheme", ParameterString, &implicit_scheme_name );
    get_parameter( "TimeDisc/Implicit/method", ParameterString, &method_name );
    get_parameter( "TimeDisc/Implicit/max_iter", ParameterDigit, &max_iter_inner );
    get_parameter( "TimeDisc/Implicit/solver", ParameterString, &solver_name );
    get_parameter( "TimeDisc/Implicit/jacobian_type", ParameterString, &jacobian_type_name );
    get_parameter( "TimeDisc/Implicit/tolerance_lsoe", ParameterNumber, &tolerance_lsoe );
    get_parameter( "TimeDisc/Implicit/max_iter_lsoe", ParameterDigit, &max_iter_lsoe );
    get_parameter( "TimeDisc/Implicit/max_krylov_dims", ParameterDigit, &max_krylov_dims );
    get_parameter( "TimeDisc/Implicit/max_krylov_restarts", ParameterDigit, &max_krylov_restarts );

        // select case( to_lower( strip( scheme ) ) )
        //     case( 'euler' )
        //         n_stages = 1
        //         allocate( phi_old(n_variables,n_cells,n_stages) ); phi_old=0.0
        //         allocate( bdf_a(n_stages) ); bdf_a = bdf_a_euler
        //         bdf_b = bdf_b_euler
        //     case( 'bdf-2' )
        //         n_stages = 2
        //         allocate( phi_old(n_variables,n_cells,n_stages) ); phi_old=0.0
        //         allocate( bdf_a(n_stages) ); bdf_a = bdf_a_bdf2
        //         bdf_b = bdf_b_bdf2
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown implicit scheme selected (' // set_string( scheme ) // ')!' )
        // end select

        // select case( to_lower( strip( method ) ) )
        //     case( 'newton' )
        //         time_step_routine => time_step_newton
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown implicit method selected (' // set_string( method ) // ')!' )
        // end select

        // select case( to_lower( strip( jacobian_type ) ) )
        //     case( 'numerical' )
        //         matrix_vector_d => matrix_vector_numerical
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown jacobian type selected (' // set_string( jacobian_type ) // ')!' )
        // end select

        // select case( to_lower( strip( solver ) ) )
        //     case( 'gmres' )
        //         solve_lsoe => GMRes_d
        //     case( 'bicgstab' )
        //         solve_lsoe => BiCGStab_d
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown linear solver selected (' // set_string( solver ) // ')!' )
        // end select

        // allocate( Y_n(n_variables,n_cells) ); Y_n = 0.0
        // allocate( Y_dt_n(n_variables,n_cells) ); Y_dt_n = 0.0

}

void implicit_finalize()
{
    deallocate( implicit_scheme_name );
    deallocate( method_name);
    deallocate( solver_name );
    deallocate( jacobian_type_name );

    // _DEALLOCATE( phi_old )
    // _DEALLOCATE( bdf_a )

    // _DEALLOCATE( Y_n )
    // _DEALLOCATE( Y_dt_n )
}

void time_step_newton( double t, double dt, int iter )
{
        // ! timestep
        // loc_dt  = dt
        // loc_t   = t + loc_dt

        // ! discretization parameters
        // loc_n_stages = n_stages; if( iter < n_stages ) loc_n_stages = 1
        // loc_bdf_a = bdf_a; if( iter < n_stages ) loc_bdf_a = bdf_a_euler
        // loc_bdf_b = bdf_b; if( iter < n_stages ) loc_bdf_b = bdf_b_euler

        // ! store the old state before calculating the FVTimeDerivative
        // do i_stage = n_stages, 2, -1
        //     phi_old(:,:,i_stage) = phi_old(:,:,i_stage-1)
        // end do
        // phi_old(:,:,1) = phi(:,:n_cells)

        // ! fill inital values for newton iteration
        // Y_n  = phi(:,:n_cells)
        // Y_dt_n = phi_dt(:,:n_cells)

        // ! calculate the inital error for newton abort criterion
        // fY_n = -(Y_n - phi_old(:,:,1)) / loc_dt + loc_bdf_b * Y_dt_n
        // do i_stage = 1, loc_n_stages
        //     fY_n = fY_n - loc_bdf_a(i_stage) / loc_dt * phi_old(:,:,i_stage)
        // end do

        // err_fY_0 = sqrt( array_dot_product( fY_n, fY_n ) )

        // do n_iter_inner = 1, max_iter
        //     ! Jac * dY = fY_n => dY ... Jacobian is determined via finite difference. fY_n = phi - dt * RHS
        //     n_iter_lsoe = max_iter_lsoe; residual_lsoe = tolerance_lsoe
        //     call solve_lsoe( n_variables, n_cells, fY_n, dY, n_iter_lsoe, residual_lsoe )

        //     Y_n = Y_n + dY  ! Y^(n+1) = Y^(n) + (Y^(n+1)-Y^(n))
        //     phi(:,:n_cells) = Y_n

        //     call fv_time_derivative( loc_t )

        //     Y_dt_n = phi_dt(:,:n_cells)
        //     fY_n = -(Y_n - phi_old(:,:,1)) / loc_dt + loc_bdf_b * Y_dt_n
        //     do i_stage = 1, loc_n_stages
        //         fY_n = fY_n - loc_bdf_a(i_stage) / loc_dt * phi_old(:,:,i_stage)
        //     end do

        //     err_fY = sqrt( array_dot_product( fY_n, fY_n ) )
        //     if( err_fY .le. err_fY_0 ) return
        //     if( .not. transient ) return
        // end do

        // ! check if iteration is lower than maximum (not converged)
        // call add_error( __LINE__, __FILE__, &
        //     'Newton not converged (error:' // set_string( err_fY ) // &
        //     ',limit:' // set_string( err_fY_0 ) // ')!' )

}

void matrix_vector_numerical( double *x, double *b, int n, int m )
{
        // do i = 1, n
        //     epsFD   = sqrt( dot_product( Y_n(i,:), Y_n(i,:) ) ) * rEps0
        //     phi(:,:n_cells) = Y_n

        //     ! positive slope
        //     phi(i,:n_cells) = Y_n(i,:) + 0.5 * epsFD
        //     call fv_time_derivative( loc_t )
        //     jacobian(:,i,:) = phi_dt(:,:n_cells)

        //     ! negative slope
        //     phi(i,:n_cells) = Y_n(i,:) - 0.5 * epsFD
        //     call fv_time_derivative( loc_t )
        //     jacobian(:,i,:) = -(jacobian(:,i,:) - phi_dt(:,:n_cells))

        //     ! combine slopes
        //     jacobian(:,i,:) = jacobian(:,i,:) * loc_bdf_b / (epsFD + _SMALL_)
        //     jacobian(i,i,:) = jacobian(i,i,:) + 1.0 / loc_dt
        // end do

        // do i = 1, m
        //     b(:,i) = matmul( jacobian(:,:,i), x(:,i) )
        // end do
}