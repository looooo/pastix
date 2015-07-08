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
#include "z_csc.h"
#include "bcsc.h"
#include "z_bcsc.h"
// #include "s_raff_functions.h"
// #include "d_raff_functions.h"
// #include "c_raff_functions.h"
#include "z_raff_functions.h"
// #include "s_csc_intern_updown.h"
// #include "d_csc_intern_updown.h"
// #include "c_csc_intern_updown.h"
#include "z_csc_intern_updown.h"

static void (*sopalinRaff[4][4])(SolverMatrix*, SopalinParam*) = 
{
//  API_RAF_GMRES
    {
        s_gmres_thread,
        d_gmres_thread,
        c_gmres_thread,
        z_gmres_thread},
//  API_RAF_PIVOT
    {
        s_pivot_thread,
        d_pivot_thread,
        c_pivot_thread,
        z_pivot_thread},
//  API_RAF_GRAD
    {
        s_grad_thread,
        d_grad_thread,
        c_grad_thread,
        z_grad_thread},
//  API_RAF_BICGSTAB
    {
        s_bicgstab_thread,
        d_bicgstab_thread,
        c_bicgstab_thread,
        z_bicgstab_thread}
};

static int (*bcscApplyPerm[4])(pastix_int_t, pastix_int_t, void*, pastix_int_t, pastix_int_t*) = 
{
    s_bcscApplyPerm,
    d_bcscApplyPerm,
    c_bcscApplyPerm,
    z_bcscApplyPerm
};

/*
 Function: pastix_task_raff

 Reffinement task

 Parameters:
 pastix_data - PaStiX data structure.
 pastix_comm - PaStiX MPI communicator.
 n           - Matrix size.
 b           - Right hand side.
 loc2glob    - local to global column number.
 */
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
 *          The PaStiX data structure.
 * 
 * @param[in] pastix_comm
 *          The PaStiX MPI communicator.
 * 
 * @param[in] n
 *          The matrix size.
 * 
 * @param[in,out] rhs
 *          The right hand side members.
 * 
 * @param[in] rhsnbr
 *          The number of right hand side members.
 * 
 * @param[in] loc2glob
 *          The local to global column number.
 *
 *******************************************************************************/
void pastix_task_raff(pastix_data_t *pastix_data,
                      MPI_Comm       pastix_comm,
                      pastix_int_t   n,
                      void          *b,
                      pastix_int_t   rhsnbr,
                      pastix_int_t  *loc2glob)
{
    pastix_int_t           * iparm    = pastix_data->iparm;
    double        * dparm    = pastix_data->dparm;
    SopalinParam  * sopar    = &(pastix_data->sopar);
    SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
    Order         * ordemesh = pastix_data->ordemesh;
    double          srafftime,rrafftime;
    pastix_int_t             procnum  = pastix_data->procnum;;
    void          * tmp;

    print_debug(DBG_STEP, "->pastix_task_raff\n");

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_onempi("%s", OUT_STEP_REFF);

    if (sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if (procnum == 0)
            errorPrintW("Refinment step incompatible with 2D distribution");
        return;
    }

    if (rhsnbr > 1)
    {
        errorPrintW("Reffinement works only with 1 rhs, please call them one after the other.");
        solvmatr->updovct.sm2xnbr = 1;
    }

    if (iparm[IPARM_ONLY_RAFF] == API_YES )
    {

        /* setting sopar->b for reffinement */
        if (sopar->b == NULL)
        {

        if( PASTIX_SUCCESS != bcscApplyPerm[iparm[IPARM_FLOAT]-2]( pastix_data->bcsc->gN,
                                                                   1,
                                                                   b,
                                                                   pastix_data->bcsc->gN,
                                                                   pastix_data->ordemesh->permtab )
        {
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;
        }
        
        MALLOC_INTERN(sopar->b, solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze * pastix_size_of( iparm[IPARM_FLOAT] ), char );
        memcpy(sopar->b, b, pastix_data->bcsc->gN * pastix_size_of( iparm[IPARM_FLOAT] ));
    }

    sopar->itermax     = iparm[IPARM_ITERMAX];
    sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
    sopar->gmresim     = iparm[IPARM_GMRES_IM];
    
    if(sopalinRaff[iparm[IPARM_REFINEMENT]][iparm[IPARM_FLOAT]-2](solvmatr, sopar) == NULL)
    {
        errorPrint("Refinement method and factorization type are incompatibles");
        iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
        return;
    }

    dparm[DPARM_RELATIVE_ERROR] = sopar->rberror;
    iparm[IPARM_NBITER]         = sopar->itermax;

    /* sopar->b was only needed for raff */
    memFree_null(sopar->b);

    /* b <- solution */
//         CscRhsUpdown[iparm[IPARM_FLOAT]-2](&(solvmatr->updovct),
//                                            solvmatr,
//                                            b, n, ordemesh->peritab,
//                                            iparm[IPARM_DOF_NBR],
//                                            iparm[IPARM_RHS_MAKING],
//                                            pastix_comm);
        memcpy(b, &(solvmatr->updovct.sm2xtab), pastix_data->bcsc->gN * pastix_size_of( iparm[IPARM_FLOAT] ));

        if( PASTIX_SUCCESS != bcscApplyPerm[iparm[IPARM_FLOAT]-2]( pastix_data->bcsc->gN,
                                                                   1,
                                                                   b,
                                                                   pastix_data->bcsc->gN,
                                                                   pastix_data->ordemesh->peritab )
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
