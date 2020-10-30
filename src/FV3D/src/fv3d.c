//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include "fv3d_private.h"
#include "mesh/mesh_module.h"
#include "equation/equation_module.h"
#include "fv/fv_module.h"
#include "analyze/analyze_module.h"
#include "output/output_module.h"
#include "restart/restart_module.h"
#include "timedisc/timedisc_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
string_t title = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void main_define();
void main_initialize();
void main_finalize();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
int main( int argc, string_t *argv )
{
    // define the program structure
    main_define();
    mesh_define();
    equation_define();
    fv_define();
    analyze_define();
    output_define();
    restart_define();
    timedisc_define();

    // call the global initialize routine
    global_initialize( argc, argv, 1, 1 );

    // calculation
    printf_r( "\n" );
    printf_r_block( '=', "Calculation" );
    timedisc();
    printf_r_emtpy_block( '=' );

    // end the program
    check_abort( 1 );
    return 1;
}

void main_define()
{
    register_initialize_routine( main_initialize );
    register_finalize_routine( main_finalize );

    string_t tmp = "untitled";
    set_parameter( "General/title", ParameterString, &tmp, "The project title", NULL, 0 );
}

void main_initialize()
{
    get_parameter( "General/title", ParameterString, &title );
}

void main_finalize()
{
    deallocate( title );
}