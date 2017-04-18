/**
 *
 * @file pastix_task_refine.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "bcsc.h"
#include "z_refine_functions.h"
#include "c_refine_functions.h"
#include "d_refine_functions.h"
#include "s_refine_functions.h"
#include "order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Select the refinement function to call depending on the matrix type and the
 * precision
 *
 *
 *******************************************************************************/
static void (*sopalinRefine[4][4])(pastix_data_t *pastix_data, void *x, void *b) =
{
    //  API_REFINE_GMRES
    {
        s_gmres_smp,
        d_gmres_smp,
        c_gmres_smp,
        z_gmres_smp
    },
    //  API_REFINE_PIVOT
    {
        s_pivot_smp,
        d_pivot_smp,
        c_pivot_smp,
        z_pivot_smp
    },
    //  API_REFINE_GRAD
    {
        s_grad_smp,
        d_grad_smp,
        c_grad_smp,
        z_grad_smp
    },
    //  API_REFINE_BICGSTAB
    {
        s_bicgstab_smp,
        d_bicgstab_smp,
        c_bicgstab_smp,
        z_bicgstab_smp
    }
};

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * @brief Perform iterative refinement.
 *
 * This routine is affected by the following parameters:
 *   IPARM_REFINEMENT, DPARM_EPSILON_REFINEMENT
 *
 * On exit, the following parameters are set:
 *   IPARM_ERROR_NUMBER
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] x
 *          The solution vector.
 *
 * @param[in,out] rhsnbr
 *          The number of right hand side members.
 *
 * @param[in] b
 *          The right hand side member.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 *
 *******************************************************************************/
int
pastix_task_refine( pastix_data_t *pastix_data,
                  void          *x,
                  pastix_int_t   rhsnbr,
                  void          *b )
{
    pastix_int_t  *iparm    = pastix_data->iparm;
    Order         *ordemesh = pastix_data->ordemesh;
    double timer;

    print_debug(DBG_STEP, "->pastix_task_refine\n");

    if (rhsnbr > 1)
    {
        errorPrintW("Refinement works only with 1 rhs, please call them one after the other.");
        rhsnbr = 1;
    }

    /* Prepare the refinment threshold, if not set by the user */
    if ( pastix_data->dparm[DPARM_EPSILON_REFINEMENT] < 0. ) {
        if ( (pastix_data->bcsc->flttype == PastixFloat) ||
             (pastix_data->bcsc->flttype == PastixComplex32) )
            pastix_data->dparm[DPARM_EPSILON_REFINEMENT] = 1e-6;
        else
            pastix_data->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
    }

    if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                         1,
                                         b,
                                         pastix_data->bcsc->gN,
                                         ordemesh->permtab ))
    {
        iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_BADPARAMETER;
        return PASTIX_ERR_BADPARAMETER;
    }

    if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                         1,
                                         x,
                                         pastix_data->bcsc->gN,
                                         ordemesh->permtab ))
    {
        iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_BADPARAMETER;
        return PASTIX_ERR_BADPARAMETER;
    }

    clockStart(timer);
    sopalinRefine[iparm[IPARM_REFINEMENT]][pastix_data->bcsc->flttype -2](pastix_data, x, b);
    clockStop(timer);
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
        pastix_print( 0, 0, OUT_TIME_REFINE, clockVal(timer) );
    }

    if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                         1,
                                         b,
                                         pastix_data->bcsc->gN,
                                         ordemesh->peritab ))
    {
        iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_BADPARAMETER;
        return PASTIX_ERR_BADPARAMETER;
    }

    if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                         1,
                                         x,
                                         pastix_data->bcsc->gN,
                                         ordemesh->peritab ))
    {
        iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_BADPARAMETER;
        return PASTIX_ERR_BADPARAMETER;
    }

    return PASTIX_SUCCESS;
}
