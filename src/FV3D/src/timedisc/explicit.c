//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "explicit_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
int explicit_active = 0;

string_t explicit_scheme_name = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void explicit_initialize();
void explicit_finalize();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void explicit_define()
{
    register_initialize_routine( explicit_initialize );
    register_finalize_routine( explicit_finalize );

    string_t tmp_opt[] = {"RK3-3", "RK4-5", "Euler"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    set_parameter( "TimeDisc/Explicit/scheme", ParameterString, &tmp, "The explicit timestep scheme", &tmp_opt, tmp_opt_n );
}

void explicit_initialize()
{
    if (explicit_active == 0) return;

    get_parameter( "TimeDisc/Explicit/scheme", ParameterString, &explicit_scheme_name );

        // time_step_routine => time_step_lserkw2

        // select case( to_lower( strip( scheme ) ) )
        //     case( 'euler' )
        //         n_rk_stages = 1
        //         allocate(rk_a(1:n_rk_stages))
        //         allocate(rk_b(1:n_rk_stages))
        //         allocate(rk_g(1:n_rk_stages))

        //         rk_a = [0.]
        //         rk_b = [0.]
        //         rk_g = [1.]
        //     case( 'rk3-3' )
        //         n_rk_stages = 3
        //         allocate(rk_a(1:n_rk_stages))
        //         allocate(rk_b(1:n_rk_stages))
        //         allocate(rk_g(1:n_rk_stages))

        //         rk_a = [0., -5./9., -153./128.]
        //         rk_b = [0., 1./3., 3./4.]
        //         rk_g = [1./3., 15./16., 8./15.]
        //     case( 'rk4-5' )
        //         n_rk_stages = 5
        //         allocate(rk_a(1:n_rk_stages))
        //         allocate(rk_b(1:n_rk_stages))
        //         allocate(rk_g(1:n_rk_stages))

        //         rk_a = [0., -567301805773.0/1357537059087.0, -2404267990393.0/2016746695238.0, &
        //             -3550918686646.0/2091501179385.0, -1275806237668.0/842570457699.0]
        //         rk_b = [0., 1432997174477.0/9575080441755.0, 2526269341429.0/6820363962896.0, &
        //             2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0]
        //         rk_g = [1432997174477.0/9575080441755.0, 5161836677717.0/13612068292357.0, &
        //             1720146321549.0/2090206949498.0, 3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0]
        //     case default
        //         call add_error( __LINE__, __FILE__, &
        //             'Unknown explicit scheme selected (' // set_string( scheme ) // ')!' )
        // end select

}

void explicit_finalize()
{
    deallocate( explicit_scheme_name );
        //    _DEALLOCATE( rk_a )
        // _DEALLOCATE( rk_b )
        // _DEALLOCATE( rk_g )
}

void time_step_lserkw2( double t, double dt, int iter )
{
        // phi_dt_tmp = 0.00

        // ! first stage
        // t_stage = t ! + dt * rk_b(i_stage) = 0
        // call fv_time_derivative( t_stage )

        // phi_dt_tmp = phi_dt(:,:n_cells) ! + phi_dt_tmp * rk_a(i_stage) = 0
        // phi(:,:n_cells) = phi(:,:n_cells) + phi_dt_tmp * dt * rk_g(1)

        // ! 2nd to n_rk_stages
        // do i_stage = 2, n_rk_stages
        //     t_stage = t + dt * rk_b(i_stage)
        //     call fv_time_derivative( t_stage )

        //     phi_dt_tmp = phi_dt(:,:n_cells) + phi_dt_tmp * rk_a(i_stage)
        //     phi(:,:n_cells) = phi(:,:n_cells) + phi_dt_tmp * dt * rk_g(i_stage)
        // end do

}