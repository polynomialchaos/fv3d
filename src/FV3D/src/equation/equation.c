//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "equation_private.h"
#include "navier_stokes/navier_stokes_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
string_t equation_name = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_initialize();
void equation_finalize();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void equation_define()
{
    register_initialize_routine( equation_initialize );
    register_finalize_routine( equation_finalize );

    string_t tmp_opt[] = {"Navier-Stokes"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    set_parameter( "Equation/equation", ParameterString, &tmp,
        "The limiter method.", &tmp_opt, tmp_opt_n );

    navier_stokes_define();
}

void equation_initialize()
{
    get_parameter( "Equation/equation", ParameterString, &equation_name );

    if (is_equal( equation_name, "Navier-Stokes" ))
    {
        navier_stokes_active = 1;
    }
    else
    {
        check_error( 0 );
    }
}

void equation_finalize()
{
    deallocate( equation_name );

        //     _DEALLOCATE( variables )

        // nullify( exact_func_routine )
        // nullify( update_routine )
        // nullify( update_gradients_routine )
        // nullify( calc_time_step_routine )
        // nullify( calc_flux_routine )
}

int add_variable( string_t name, int is_active )
{
        // if( allocated( variables ) ) then
        //     if( .not. variables(n_tot_variables)%is_active .and. is_active ) &
        //         call add_error( __LINE__, __FILE__, &
        //             'Cannot add active variable after inactive variable (' // set_string( name ) // ')!' )

        //     allocate( tmp( n_tot_variables ) )
        //     do i = 1, n_tot_variables
        //         tmp(i) = variables(i)
        //     end do

        //     _DEALLOCATE( variables )
        // end if

        // if( is_active ) n_variables = n_variables + 1
        // n_tot_variables = n_tot_variables + 1
        // allocate( variables(n_tot_variables) )

        // do i = 1, n_tot_variables - 1
        //     variables(i) = tmp(i)
        // end do

        // variables(n_tot_variables)%name       = name
        // variables(n_tot_variables)%is_active  = is_active

        // res = n_tot_variables
    return 9;
}