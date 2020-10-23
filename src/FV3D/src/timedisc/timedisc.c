
//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "timedisc_private.h"
#include "explicit/explicit_private.h"
#include "implicit/implicit_private.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
string_t time_step_name = NULL;
int max_iter = 1000;
int transient = 1;
double abort_residual = 1e-10;
double t_start = 0.0;
double t_end = 0.2;

void_timedisc_fp_t time_step = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void timedisc_initialize();
void timedisc_finalize();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void timedisc_define()
{
    register_initialize_routine( timedisc_initialize );
    register_finalize_routine( timedisc_finalize );

    string_t tmp_opt[] = {"Explicit", "Implicit"};
    int tmp_opt_n = sizeof( tmp_opt ) / sizeof( string_t );
    string_t tmp = tmp_opt[0];

    set_parameter( "TimeDisc/time_step", ParameterString, &tmp, "The timestep mehtod", &tmp_opt, tmp_opt_n );
    set_parameter( "TimeDisc/max_iter", ParameterDigit, &max_iter, "The maximum number of iterations", NULL, 0 );
    set_parameter( "TimeDisc/transient", ParameterBool, &transient, "The flag wheter to be transient or steady-state", NULL, 0 );
    set_parameter( "TimeDisc/abort_residual", ParameterNumber, &abort_residual, "The abort residual", NULL, 0 );
    set_parameter( "TimeDisc/t_start", ParameterNumber, &t_start, "The start time", NULL, 0 );
    set_parameter( "TimeDisc/t_end", ParameterNumber, &t_end, "The end time", NULL, 0 );

    explicit_define();
    implicit_define();
}

void timedisc_initialize()
{
    get_parameter( "TimeDisc/time_step", ParameterString, &time_step_name );
    get_parameter( "TimeDisc/max_iter", ParameterDigit, &max_iter );
    get_parameter( "TimeDisc/transient", ParameterBool, &transient );
    get_parameter( "TimeDisc/abort_residual", ParameterNumber, &abort_residual );
    get_parameter( "TimeDisc/t_start", ParameterNumber, &t_start );
    get_parameter( "TimeDisc/t_end", ParameterNumber, &t_end );

    if (transient == 0) t_end = DOUBLE_MAX;

    if (is_equal( time_step_name, "Explicit" ))
    {
        explicit_active = 1;
    }
    else if (is_equal( time_step_name, "Implicit" ))
    {
        implicit_active = 1;
    }
    else
    {
        check_error( 0 );
    }
}

void timedisc_finalize()
{
    deallocate( time_step_name );
    time_step = NULL;
        // nullify( time_step_routine )
}

void timedisc()
{
        // ! check for errors
        // call check_abort()

        // do_output = .false.
        // do_finalize = .false.

        // iter = 0
        // if( use_restart ) iter = iter_restart

        // t = t_start
        // if( use_restart ) t = time_restart

        // call print_residual_header()
        // if( do_output_data .and. .not. use_restart ) call write_output( iter, t )

        // do while( .true. )
        //     ! check for errors
        //     call check_abort()

        //     ! calculate time step dt
        //     call calc_time_step_routine( dt_loc )
        //     call allreduce_mpi( dt_loc, dt, MPI_MIN )

        //     if( t+dt .ge. t_end ) then
        //         dt = t_end-t
        //         do_output = .true.
        //         do_finalize = .true.
        //     end if

        //     ! the timestep to be called
        //     call time_step_routine( t, dt, iter )
        //     call calc_global_residual( dt )

        //     ! check for NAN
        //     if( any( isnan( residual ) ) ) &
        //         call add_error( __LINE__, __FILE__, &
        //             'Detected NaN value(s) in residual!' )

        //     ! check for infinity
        //     if( any( residual .gt. infinity ) ) &
        //         call add_error( __LINE__, __FILE__, &
        //             'Detected infinite value(s) in residual!' )

        //     t = t + dt
        //     iter = iter + 1

        //     ! steady-state simulation
        //     if( ( .not. transient ) .and. ( all( residual .lt. abort_residual ) ) ) then
        //         t_end = t
        //         do_output = .true.
        //         do_finalize = .true.
        //     end if

        //     if( ( i_output_data .gt. 0 ) .and. ( mod( iter, i_output_data ) .eq. 0 ) ) then
        //         do_output = .true.
        //     end if

        //     ! maximum iteration number reacher
        //     if( iter .ge. iter_restart + max_iter ) then
        //         t_end = t
        //         do_output = .true.
        //         do_finalize = .true.
        //     end if

        //     call print_residual( do_output_data .and. do_output )

        //     ! output
        //     if( do_output_data .and. do_output ) then
        //         call write_output( iter, t )
        //         do_output = .false.
        //     end if

        //     ! stop if required
        //     if( do_finalize ) exit
        // end do

        // if( iter .ge. max_iter ) &
        //     call add_warning( __LINE__, __FILE__, &
        //         'Iteration limit exceeded (' // set_string( max_iter ) // ')!' )

}

void print_residual_header()
{
        // def_out = get_def_out()

        // if( is_explicit ) then
        //     write( def_out, '(A9,1X,A12,1X,A12,1X,A1,1X,A1,1X,":")', advance='no' ) &
        //         'iter', 'time', 'dt', 'V', 'O'
        // else
        //     write( def_out, '(A9,1X,A12,1X,A12,1X,A1,1X,A1,1X,A6,1X,A6,1X,":")', advance='no' ) &
        //         'iter', 'time', 'dt', 'V', 'O', 'inner', 'lsoe'
        // end if
        // write( def_out, '(' // set_string( n_variables ) //'(1X,A12))' ) (variables(i)%name, i = 1, n_variables)

        // if( def_out .ne. _STDOUT_ ) call flush( def_out )

}

void print_residual( int do_output )
{
        // def_out = get_def_out()

        // output_str = ' '
        // if( do_output ) output_str = '*'

        // if( is_explicit ) then
        //     write( def_out, '(I9,1X,ES12.5,1X,ES12.5,1X,L1,1X,A1,1X,":")', advance='no' ) &
        //         iter, t, dt, is_viscous_dt, output_str
        // else
        //     write( def_out, '(I9,1X,ES12.5,1X,ES12.5,1X,L1,1X,A1,1X,I6,1X,I6,1X,":")', advance='no' ) &
        //         iter, t, dt, is_viscous_dt, output_str, n_iter_inner, n_iter_lsoe
        // end if
        // write( def_out, '(' // set_string( n_variables ) //'(1X,ES12.5))' ) residual

        // if( def_out .ne. _STDOUT_ ) call flush( def_out )
}