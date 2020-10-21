//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "limiter_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
void_limiter_fp_t limiter_routine = NULL;

string_t limiter_name = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void limiter_initialize();
void limiter_finalize();
void limiter_none( size_t i, double *slope, double *lim, int n );
void limiter_barth_jespersenn( size_t i, double *slope, double *lim, int n );

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
        "The limiter method.", &tmp_opt, tmp_opt_n );
}

void limiter_initialize()
{
    get_parameter( "FV/Limiter/limiter", ParameterString, &limiter_name );

    if (is_equal( limiter_name, "None" ))
    {
        limiter_routine = limiter_none;
    }
    else if (is_equal( limiter_name, "Barth-Jespersenn" ))
    {
        limiter_routine = limiter_barth_jespersenn;
    }
    else
    {
        check_error( 0 );
    }
}

void limiter_finalize()
{
    deallocate( limiter_name );
    limiter_routine = NULL;
}

void limiter_none( size_t i, double *slope, double *lim, int n )
{
    set_value_n( 1.0, lim, n );
}

void limiter_barth_jespersenn( size_t i, double *slope, double *lim, int n )
{
    set_value_n( 1.0, lim, n );
}
//  !###############################################################################################################################
//     !> Calculate Barth-Jespersenn's limiter.
//     !-------------------------------------------------------------------------------------------------------------------------------
//     subroutine limiter_barth_jespersenn( i_cell, slope, lim )
//         use mod_equation_vars,  only: n_tot_variables
//         use mod_mesh_vars,      only: cells, faces
//         use mod_fv_vars,        only: phi_total
//         implicit none
//         !---------------------------------------------------------------------------------------------------------------------------
//         integer,    intent(in)  :: i_cell
//         real,       intent(in)  :: slope(n_tot_variables)
//         real,       intent(out) :: lim(n_tot_variables)
//         !---------------------------------------------------------------------------------------------------------------------------
//         integer         :: i, j, k
//         real            :: y, u_min(n_tot_variables), u_max(n_tot_variables)
//         !---------------------------------------------------------------------------------------------------------------------------

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

//     end subroutine limiter_barth_jespersenn

// end module mod_fv_limiter