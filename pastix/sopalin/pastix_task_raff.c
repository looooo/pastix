/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
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
/*
 * File: pastix_task_raff.c
 *
 * PaStiX reffinement functions implementations.
 *
 */

#include "common.h"
#include "csc.h"
#include "bcsc.h"

static void (*sopalinRaff[4])(SolverMatrix*, SopalinParam*) = 
{
//  API_RAF_GMRES
    gmres_thread,
//  API_RAF_PIVOT
    pivot_thread,
//  API_RAF_GRAD
    grad_thread,
//  API_RAF_BICGSTAB
    bicgstab_thread
};

static void (*buildUpdoVect[4])(pastix_data_t*, pastix_int_t*, void*, MPI_Comm) = 
{
    s_buildUpdoVect,
    d_buildUpdoVect,
    c_buildUpdoVect,
    z_buildUpdoVect
};

static void (*CscRhsUpdown[4])(UpDownVector*, SolverMatrix*, void*, pastix_int_t,  pastix_int_t, int, int, MPI_Comm) = 
{
    s_CscRhsUpdown,
    d_CscRhsUpdown,
    c_CscRhsUpdown,
    z_CscRhsUpdown
};

static void (*CscRhsUpdown[4])(UpDownVector*, SolverMatrix*, void*, pastix_int_t, pastix_int_t,  pastix_int_t, int, int, MPI_Comm) = 
{
    s_CscdRhsUpdown,
    d_CscdRhsUpdown,
    c_CscdRhsUpdown,
    z_CscdRhsUpdown
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
#ifndef PASTIX_UPDO_ISEND
    if (THREAD_COMM_ON)
    {
        if (procnum == 0)
            errorPrintW("THREAD_COMM require -DPASTIX_UPDO_ISEND,"
                        " force API_THREAD_MULTIPLE");
        sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* PASTIX_UPDO_ISEND */
#ifndef STORAGE
    if (THREAD_COMM_ON)
    {
        if (procnum == 0)
            errorPrintW("THREAD_COMM require -DSTORAGE,"
                        " force API_THREAD_MULTIPLE");
        sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* STORAGE */

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
            switch(iparm[IPARM_FLOAT])
            {
            case API_REALSINGLE:
                MALLOC_INTERN(sopar->b,
                              solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                              pastix_float_t);
            break;
            case API_REALDOUBLE:
                MALLOC_INTERN(sopar->b,
                              solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                              pastix_double_t);
            break;
            case API_COMPLEXSINGLE:
                MALLOC_INTERN(sopar->b,
                              solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                              pastix_complex32_t);
            break;
            case API_COMPLEXDOUBLE:
            default:
                MALLOC_INTERN(sopar->b,
                              solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                              pastix_complex64_t);
                
            }
        }

        tmp = solvmatr->updovct.sm2xtab;
        solvmatr->updovct.sm2xtab = sopar->b;

        buildUpdoVect[iparm[IPARM_FLOAT]-2](pastix_data,
                                          loc2glob,
                                          b,
                                          pastix_comm);

        sopar->b = solvmatr->updovct.sm2xtab;
        solvmatr->updovct.sm2xtab = tmp;

    }

    sopar->itermax     = iparm[IPARM_ITERMAX];
    sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
#ifdef OOC
    if (iparm[IPARM_GMRES_IM] != 1)
    {
        iparm[IPARM_GMRES_IM] = 1;
        if (procnum == 0)
            errorPrintW("IPARM_GMRES_IM force to 1 when using OOC");
    }
#endif
    sopar->gmresim = iparm[IPARM_GMRES_IM];
    
    if(sopalinRaff[iparm[IPARM_REFINEMENT]](solvmatr, sopar) == NULL)
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
    if (iparm[IPARM_GRAPHDIST] == API_NO)
    {
        CscRhsUpdown[iparm[IPARM_FLOAT]-2](&(solvmatr->updovct),
                                           solvmatr,
                                           b, n, ordemesh->peritab,
                                           iparm[IPARM_DOF_NBR],
                                           iparm[IPARM_RHS_MAKING],
                                           pastix_comm);
    }
#ifdef PASTIX_DISTRIBUTED
    else
    {
        CscdRhsUpdown[iparm[IPARM_FLOAT]-2](&(solvmatr->updovct),
                                            solvmatr,
                                            b, n,
                                            pastix_data->glob2loc,
                                            ordemesh->peritab,
                                            iparm[IPARM_DOF_NBR], pastix_comm);
    }
#endif

    /* Fin du roerdering */

    srafftime = (double)dparm[DPARM_RAFF_TIME];
    MPI_Reduce(&srafftime,&rrafftime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

    if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
    {
        fprintf(stdout, OUT_RAFF_ITER_NORM, (long)iparm[IPARM_NBITER], (double)dparm[DPARM_RELATIVE_ERROR]);
        if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
            if (dparm[DPARM_RELATIVE_ERROR] > 0)
                print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
            if (dparm[DPARM_SCALED_RESIDUAL] > 0)
                print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
        }

        fprintf(stdout, OUT_TIME_RAFF, rrafftime);
    }
    iparm[IPARM_START_TASK]++;

    return;
} 
