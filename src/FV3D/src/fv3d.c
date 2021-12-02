/*******************************************************************************
 * @file fv3d.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include "fv3d_private.h"
#include "mesh/mesh_module.h"
#include "navier_stokes/navier_stokes_module.h"
#include "fv/fv_module.h"
#include "analyze/analyze_module.h"
#include "output/output_module.h"
#include "restart/restart_module.h"
#include "timedisc/timedisc_module.h"

/*******************************************************************************
 * @brief Main function
 * @param argc
 * @param argv
 * @return int
 ******************************************************************************/
int main(int argc, string_t *argv)
{
    /* define the program structure */
    fv3d_define();
    mesh_define();
    navier_stokes_define();
    fv_define();
    analyze_define();
    output_define();
    timedisc_define();
    restart_define();

    /* call the global initialize routine */
    global_initialize(argc, argv, BTRU, BTRU, BTRU, BTRU);

    /* mesh info */
    PRINTF("\n");
    printf_r_sep_title('=', "Mesh");
    print_mesh_info();
    printf_r_sep('=');

    /* equation info */
    PRINTF("\n");
    printf_r_sep_title('=', "Variables");
    print_variables();
    printf_r_sep('=');

    /* calculation */
    PRINTF("\n");
    printf_r_sep_title('=', "Calculation");
    timedisc();
    printf_r_sep('=');

    /* end the program */
    check_abort(1);
    return 1;
}

/*******************************************************************************
 * @brief Define fv3d
 ******************************************************************************/
void fv3d_define()
{
    REGISTER_INITIALIZE_ROUTINE(fv3d_initialize);
    REGISTER_FINALIZE_ROUTINE(fv3d_finalize);

    string_t title = "untitled";
    SET_PARAMETER("General/title", StringParameter, &title,
                  "The project title", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize fv3d
 ******************************************************************************/
void fv3d_finalize()
{
    free_solver();
}

/*******************************************************************************
 * @brief Initialize fv3d
 ******************************************************************************/
void fv3d_initialize()
{
    string_t title = NULL;
    GET_PARAMETER("General/title", StringParameter, &title);
    init_solver(title);
    DEALLOCATE(title);
}
