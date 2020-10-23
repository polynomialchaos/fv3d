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

double limiter_none( int i_cell, int i_var, double slope );
double limiter_barth_jespersenn( int i_cell, int i_var, double slope );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void limiter_define()
{
    register_initialize_routine( limiter_initialize );
    register_finalize_routine( limiter_finalize );

    string_t tmp_opt[] = {"Barth-Jespersenn", "None"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    set_parameter( "FV/Limiter/limiter", ParameterString, &tmp,
        "The limiter method", &tmp_opt, tmp_opt_n );
}

void limiter_initialize()
{
    get_parameter( "FV/Limiter/limiter", ParameterString, &limiter_name );

    if (is_equal( limiter_name, "None" ))
    {
        limiter_function_pointer = limiter_none;
    }
    else if (is_equal( limiter_name, "Barth-Jespersenn" ))
    {
        limiter_function_pointer = limiter_barth_jespersenn;
    }
    else
    {
        check_error( 0 );
    }
}

void limiter_finalize()
{
    deallocate( limiter_name );
    limiter_function_pointer = NULL;
}

double limiter_none( int i_cell, int i_var, double slope )
{
#ifdef DEBUG
    u_unused( i_cell );
    u_unused( i_var );
    u_unused( slope );
#endif /* DEBUG */
    return 1.0;
}

double limiter_barth_jespersenn( int i_cell, int i_var, double slope )
{
    int n_tot_variables = all_variables->n_tot_variables;
    double u_min        = phi_total[i_cell*n_tot_variables+i_var];
    double u_max        = phi_total[i_cell*n_tot_variables+i_var];

    return 1.0;

    //         u_min = phi_total(:,i_cell)
//         u_max = phi_total(:,i_cell)

//         do i = 1, cells(i_cell)%n_faces
//             j = cells(i_cell)%faces(i)

//             if( faces(j)%cells(1) .eq. i_cell ) then
//                 u_min = min( u_min, phi_total(:,faces(j)%cells(2)) )
//                 u_max = max( u_max, phi_total(:,faces(j)%cells(2)) )
//             else
//                 u_min = min( u_min, phi_total(:,faces(j)%cells(1)) )
//                 u_max = max( u_max, phi_total(:,faces(j)%cells(1)) )
//             end if
//         end do

//         lim = 1.0

//         do i = 1, cells(i_cell)%n_faces
//             j = cells(i_cell)%faces(i)

//             do k = 1, n_tot_variables
//                 if( slope(k) .gt. 0 ) then
//                     y = (u_max(k) - phi_total(k,i_cell)) / slope(k)
//                 else if( slope(k) .lt. 0 ) then
//                     y = (u_min(k) - phi_total(k,i_cell)) / slope(k)
//                 else
//                     y = 1.0
//                 end if

//                 lim(k) = min( lim(k), y )
//             end do
//         end do
}