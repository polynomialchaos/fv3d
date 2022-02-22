/*******************************************************************************
 * @file solver.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "solver_private.h"

string_t solver_title = NULL;

/*******************************************************************************
 * @brief Free solver
 ******************************************************************************/
void free_solver()
{
    DEALLOCATE(solver_title);
}

/*******************************************************************************
 * @brief Return the simulation title
 ******************************************************************************/
cstring_t get_simulation_title()
{
    return solver_title;
}

/*******************************************************************************
 * @brief Initialize solver
 * @param title
 ******************************************************************************/
void init_solver(cstring_t title)
{
    solver_title = allocate_strcpy(title);
}
