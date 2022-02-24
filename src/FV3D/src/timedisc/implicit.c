/*******************************************************************************
 * @file implicit.c
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2022-02-22
 * @copyright Copyright (c) 2022 by Florian Eigentler.
 *  This work is licensed under terms of the MIT license (<LICENSE>).
 ******************************************************************************/
#include "timedisc/timedisc_module.h"

string_t implicit_scheme_name = NULL;
int max_iter_inner = 100;
string_t solver_name = NULL;
double tolerance_lsoe = 1e-12;
int max_iter_lsoe = 100;
int max_krylov_dims = 15;
int max_krylov_restarts = 2;

/*******************************************************************************
 * @brief Define implicit timedisc
 ******************************************************************************/
void implicit_define()
{
    BM_REGISTER_INITIALIZE_ROUTINE(implicit_initialize);
    BM_REGISTER_FINALIZE_ROUTINE(implicit_finalize);

    string_t tmp_opt[] = {"BDF-2", "Euler"};
    int tmp_opt_n = sizeof(tmp_opt) / sizeof(string_t);
    string_t tmp = tmp_opt[0];

    string_t tmp2_opt[] = {"BiCGStab", "GMRes"};
    int tmp2_opt_n = sizeof(tmp2_opt) / sizeof(string_t);
    string_t tmp2 = tmp2_opt[0];

    BM_SET_PARAMETER("TimeDisc/Implicit/scheme", StringParameter, &tmp,
                     "The implicit timestep scheme", &tmp_opt, tmp_opt_n);
    BM_SET_PARAMETER("TimeDisc/Implicit/max_iter_inner", DigitParameter,
                     &max_iter_inner,
                     "The maximum number of inner iterations", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/Implicit/solver", StringParameter, &tmp2,
                     "The linear system of equations solver",
                     &tmp2_opt, tmp2_opt_n);
    BM_SET_PARAMETER("TimeDisc/Implicit/tolerance_lsoe", NumberParameter,
                     &tolerance_lsoe,
                     "The linear solver tolerance", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/Implicit/max_iter_lsoe", DigitParameter,
                     &max_iter_lsoe,
                     "The linear solver maximum number of iterations", NULL, 0);
    BM_SET_PARAMETER("TimeDisc/Implicit/max_krylov_dims", DigitParameter,
                     &max_krylov_dims,
                     "The maximum Krylov space dimension in GMRes solver",
                     NULL, 0);
    BM_SET_PARAMETER("TimeDisc/Implicit/max_krylov_restarts", DigitParameter,
                     &max_krylov_restarts,
                     "The maximum restarts performed in GMRes solver", NULL, 0);
}

/*******************************************************************************
 * @brief Finalize implicit timedisc
 ******************************************************************************/
void implicit_finalize()
{
    BM_DEALLOCATE(implicit_scheme_name);
    BM_DEALLOCATE(solver_name);
    free_implicit();
}

/*******************************************************************************
 * @brief Initialize implicit timedisc
 ******************************************************************************/
void implicit_initialize()
{
    if (is_explicit() == BC_TRUE)
        return;

    BM_GET_PARAMETER("TimeDisc/Implicit/scheme", StringParameter,
                     &implicit_scheme_name);
    BM_GET_PARAMETER("TimeDisc/Implicit/max_iter_inner", DigitParameter,
                     &max_iter_inner);
    BM_GET_PARAMETER("TimeDisc/Implicit/solver", StringParameter, &solver_name);
    BM_GET_PARAMETER("TimeDisc/Implicit/tolerance_lsoe", NumberParameter,
                     &tolerance_lsoe);
    BM_GET_PARAMETER("TimeDisc/Implicit/max_iter_lsoe", DigitParameter,
                     &max_iter_lsoe);
    BM_GET_PARAMETER("TimeDisc/Implicit/max_krylov_dims", DigitParameter,
                     &max_krylov_dims);
    BM_GET_PARAMETER("TimeDisc/Implicit/max_krylov_restarts", DigitParameter,
                     &max_krylov_restarts);

    implicit_scheme_t implicit_scheme = -1;
    if (is_equal(implicit_scheme_name, "Euler"))
    {
        implicit_scheme = EulerImplicit;
    }
    else if (is_equal(implicit_scheme_name, "BDF-2"))
    {
        implicit_scheme = BDF2;
    }
    else
    {
        BM_CHECK_EXPRESSION(0);
    }

    implicit_solver_t implicit_solver = -1;
    if (is_equal(solver_name, "BiCGStab"))
    {
        implicit_solver = BiCGStab;
    }
    else if (is_equal(solver_name, "GMRes"))
    {
        implicit_solver = GMRes;
    }
    else
    {
        BM_CHECK_EXPRESSION(0);
    }

    init_implicit(implicit_scheme, implicit_solver, max_iter_inner,
                  tolerance_lsoe, max_iter_lsoe,
                  max_krylov_dims, max_krylov_restarts);
}