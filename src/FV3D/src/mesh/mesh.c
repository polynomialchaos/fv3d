/*******************************************************************************
 * @file mesh.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#include <string.h>
#include "mesh_module.h"

string_t mesh_file = NULL;

/*******************************************************************************
 * @brief Define mesh
 ******************************************************************************/
void mesh_define()
{
    REGISTER_INITIALIZE_ROUTINE(mesh_initialize);
    REGISTER_FINALIZE_ROUTINE(mesh_finalize);

    string_t tmp = "untitled.mesh.h5";
    SET_PARAMETER("Mesh/mesh_file", StringParameter, &tmp,
                  "The mesh file", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize mesh
 ******************************************************************************/
void mesh_finalize()
{
    free_mesh();
    DEALLOCATE(mesh_file);
}

/*******************************************************************************
 * @brief Initialize mesh
 ******************************************************************************/
void mesh_initialize()
{
    GET_PARAMETER("Mesh/mesh_file", StringParameter, &mesh_file);
    init_mesh(mesh_file);
}