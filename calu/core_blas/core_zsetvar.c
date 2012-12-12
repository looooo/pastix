/**
 *
 * @file core_zsetvar.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 * CORE_zsetvar sets a single variable, x := alpha.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *         Scalar to set x to, passed by pointer so it can depend on runtime value.
 *
 * @param[out] x
 *         On exit, x = alpha.
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zsetvar = PCORE_zsetvar
#define CORE_zsetvar PCORE_zsetvar
#endif
void CORE_zsetvar(const PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *x)
{
    *x = *alpha;
}

/***************************************************************************//**
 *
 * QUARK_CORE_zsetvar sets a single variable, x := alpha.
 * Since x can be in the middle of a tile, we need to depend on the whole tile,
 * so add xlock argument.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *         Scalar to set x to, passed by pointer so it can depend on runtime value.
 *
 * @param[out] x
 *         On exit, x = alpha.
 *
 * @param[out] xlock
 *         Pointer to tile owning output variable x.
 *
 **/
void QUARK_CORE_zsetvar(Quark *quark, Quark_Task_Flags *task_flags,
                        const PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *x,
                        PLASMA_Complex64_t *xlock)
{
    DAG_SET_PROPERTIES("setvar", "orange");
    QUARK_Insert_Task(quark, CORE_zsetvar_quark, task_flags,
        sizeof(PLASMA_Complex64_t),  alpha,  INPUT,
        sizeof(PLASMA_Complex64_t),  x,              NODEP,  /* INOUT; see xlock */
        sizeof(PLASMA_Complex64_t),  xlock,          INOUT,  /* tile dependency containing x */
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zsetvar_quark = PCORE_zsetvar_quark
#define CORE_zsetvar_quark PCORE_zsetvar_quark
#endif
void CORE_zsetvar_quark(Quark *quark)
{
    const PLASMA_Complex64_t *alpha;
    PLASMA_Complex64_t *x;

    quark_unpack_args_2( quark, alpha, x );
    *x = *alpha;
}
