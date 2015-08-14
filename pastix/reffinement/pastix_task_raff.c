/**
 * 
 * @file pastix_task_raff.c
 * 
 * PaStiX reffinement functions implementations.
 *
 * PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 * LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
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
#include "z_raff_functions.h"
#include "order.h"

static void (*sopalinRaff[4][4])(pastix_data_t *pastix_data, void *x, void *b) = 
{
//  API_RAF_GMRES
    {
    s_gmres_smp,
    d_gmres_smp,
    c_gmres_smp,
    z_gmres_smp
    },
//  API_RAF_PIVOT
    {
//     s_pivot_smp,
//     d_pivot_smp,
//     c_pivot_smp,
//     z_pivot_smp
    NULL,
    NULL,
    NULL,
    NULL
    },
//  API_RAF_GRAD
    {
    s_grad_smp,
    d_grad_smp,
    c_grad_smp,
    z_grad_smp
    },
//  API_RAF_BICGSTAB
    {
//     s_bicgstab_smp,
//     d_bicgstab_smp,
//     c_bicgstab_smp,
//     z_bicgstab_smp
    NULL,
    NULL,
    NULL,
    NULL
    }
};

/**
 *******************************************************************************
 *
 * @ingroup pastix_raff
 *
 * pastix_task_raff - Reffinement task.
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
 *******************************************************************************/
void pastix_task_raff(pastix_data_t *pastix_data,
                      void          *x,
                      pastix_int_t   rhsnbr,
                      void          *b)
{
    pastix_int_t  * iparm    = pastix_data->iparm;
    double        * dparm    = pastix_data->dparm;
    SopalinParam  * sopar    = &(pastix_data->sopar);
    Order         * ordemesh = pastix_data->ordemesh;
    double timer;

    print_debug(DBG_STEP, "->pastix_task_raff\n");

    if (rhsnbr > 1)
    {
//         errorPrintW("Reffinement works only with 1 rhs, please call them one after the other.");
        rhsnbr = 1;
    }

        if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                             1,
                                             b,
                                             pastix_data->bcsc->gN,
                                             ordemesh->permtab ))
        {
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;
        }

        if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                             1,
                                             x,
                                             pastix_data->bcsc->gN,
                                             ordemesh->permtab ))
        {
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;
        }

//     sopar->itermax     = iparm[IPARM_ITERMAX];
//     sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
//     sopar->gmresim     = iparm[IPARM_GMRES_IM];
//     sopar->gN          = pastix_data->bcsc->gN;

    clockStart(timer);
    sopalinRaff[iparm[IPARM_REFINEMENT]][pastix_data->bcsc->flttype -2](pastix_data, x, b);
    clockStop(timer);
    pastix_print( 0, 0, OUT_TIME_RAFF, clockVal(timer) );

//     dparm[DPARM_RELATIVE_ERROR] = sopar->rberror;
//     iparm[IPARM_NBITER]         = sopar->itermax;

    if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                         1,
                                         b,
                                         pastix_data->bcsc->gN,
                                         ordemesh->peritab ))
    {
        iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
        return;
    }

    if( PASTIX_SUCCESS != bcscApplyPerm( pastix_data->bcsc,
                                         1,
                                         x,
                                         pastix_data->bcsc->gN,
                                         ordemesh->peritab ))
    {
        iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
        return;
    }

    /* Fin du reordering */

//     srafftime = (double)dparm[DPARM_RAFF_TIME];
//     MPI_Reduce(&srafftime,&rrafftime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

//     if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
//     {
//         fprintf(stdout, OUT_RAFF_ITER_NORM, (long)iparm[IPARM_NBITER], (double)dparm[DPARM_RELATIVE_ERROR]);
//         if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
//             if (dparm[DPARM_RELATIVE_ERROR] > 0)
//                 print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
//             if (dparm[DPARM_SCALED_RESIDUAL] > 0)
//                 print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
//         }
// 
//         fprintf(stdout, OUT_TIME_RAFF, rrafftime);
//     }
    iparm[IPARM_START_TASK]++;

    return;
}