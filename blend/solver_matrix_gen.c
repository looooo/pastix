/**
 *
 * @file solver_matrix_gen.c
 *
 * PaStiX solver structure generation function.
 *
 * @copyright 1998-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Tony Delarue
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @author Gregoire Pichon
 * @author Nolan Bredel
 * @author Alycia Lisito
 * @date 2024-07-05
 *
 **/
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>

#include "common.h"
#include "symbol/symbol.h"
#include "blend/solver.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "pastix/order.h"
#include "extendVector.h"
#include "simu.h"
#include "blendctrl.h"
#include "solver_matrix_gen_utils.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Initialize the solver matrix structure
 *
 * This function takes all the global preprocessing steps: the symbol matrix
 * and the result of the simulation step to generate the solver matrix that holds
 * only local information of each PaStiX process.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          On entry, the allocated pointer to a solver matrix structure.
 *          On exit, this structure holds alls the local information required to
 *          perform the numerical factorization.
 *
 * @param[in] symbmtx
 *          The global symbol matrix structure.
 *
 * @param[in] ordeptr
 *          The ordering structure.
 *
 * @param[in] simuctrl
 *          The information resulting from the simulation that will provide the
 *          data mapping, and the order of the task execution for the static
 *          scheduling.
 *
 * @param[in] ctrl
 *          The blend control structure that contains extra information
 *          computed during the analyze steps and the parameters of the analyze
 *          step.
 *
 * @param[in] comm
 *          TODO
 *
 * @param[in] isched
 *          TODO
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
int
solverMatrixGen( SolverMatrix          *solvmtx,
                 const symbol_matrix_t *symbmtx,
                 const pastix_order_t  *ordeptr,
                 const SimuCtrl        *simuctrl,
                 const BlendCtrl       *ctrl,
                 PASTIX_Comm            comm,
                 isched_t              *isched )
{
    pastix_int_t         i;
    pastix_int_t        *cblklocalnum;
    pastix_int_t        *bloklocalnum;
    pastix_int_t        *tasklocalnum;
    pastix_int_t        *browtmp;
    solver_cblk_recv_t **ftgttab = NULL;
    pastix_int_t        *faninnbr_tab = NULL;
    pastix_int_t         clustnbr = ctrl->clustnbr;
    (void)ordeptr;

    assert( symbmtx->baseval == 0 );

    solverInit( solvmtx );

#ifdef PASTIX_DYNSCHED
    solvmtx->btree = ctrl->btree;
#endif
    solvmtx->clustnum  = ctrl->clustnum;
    solvmtx->clustnbr  = clustnbr;
    solvmtx->procnbr   = ctrl->total_nbcores;
    solvmtx->thrdnbr   = ctrl->local_nbthrds;
    solvmtx->bublnbr   = ctrl->local_nbctxts;
    solvmtx->solv_comm = comm;

    /* Allocate the different local numbering arrays */
    MALLOC_INTERN( bloklocalnum, symbmtx->bloknbr,  pastix_int_t );
    MALLOC_INTERN( cblklocalnum, symbmtx->cblknbr,  pastix_int_t );
    MALLOC_INTERN( tasklocalnum, simuctrl->tasknbr, pastix_int_t );
    MALLOC_INTERN( ftgttab,      symbmtx->cblknbr,  solver_cblk_recv_t* );
    memset( ftgttab, 0, symbmtx->cblknbr * sizeof(solver_cblk_recv_t*) );

    if ( clustnbr > 1 ) {
        /*
         * k is the proc number
         * - faninnbr_tab[2*k]              corresponds to the first index of clbk fanin for k.
         * - faninnbr_tab[2*k+1]            corresponds to the first index of blok fanin for k.
         * - faninnbr_tab[2*clustnbr+2*k]   corresponds to the first index of clbk recv  for k.
         * - faninnbr_tab[2*clustnbr+2*k+1] corresponds to the first index of blok recv  for k.
         */
        MALLOC_INTERN( faninnbr_tab, clustnbr * 4, pastix_int_t );
        memset( faninnbr_tab, 0, clustnbr * 4 * sizeof( pastix_int_t ) );
    }

    /* Compute local indexes to compress the symbol information into solver */
    solvMatGen_fill_localnums( symbmtx, simuctrl, solvmtx,
                               cblklocalnum, bloklocalnum, tasklocalnum,
                               ftgttab, faninnbr_tab );

    solvmtx->cblkmin2d  = solvmtx->cblknbr;
    solvmtx->cblkschur  = solvmtx->cblknbr;
    solvmtx->gcblknbr   = symbmtx->cblknbr;
    solvmtx->gbloknbr   = symbmtx->bloknbr;

    /***************************************************************************
     * Fill in the local bloktab and cblktab
     */
    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN( solvmtx->cblktab,  solvmtx->cblknbr+1, SolverCblk   );
    MALLOC_INTERN( solvmtx->bloktab,  solvmtx->bloknbr+1, SolverBlok   );
    MALLOC_INTERN( solvmtx->browtab,  solvmtx->brownbr,   pastix_int_t );
    MALLOC_INTERN( browtmp,           symbmtx->browmax,   pastix_int_t );
    MALLOC_INTERN( solvmtx->gcbl2loc, symbmtx->cblknbr,   pastix_int_t );
    memset( solvmtx->gcbl2loc, 0xff,  symbmtx->cblknbr * sizeof(pastix_int_t) );
    {
        solver_cblk_recv_t *ftgtcblk;
        SolverCblk    *solvcblk = solvmtx->cblktab;
        SolverBlok    *solvblok = solvmtx->bloktab;
        symbol_cblk_t *symbcblk = symbmtx->cblktab;
        Cand          *candcblk = ctrl->candtab;
        pastix_int_t   nbcblk2d = 0;
        pastix_int_t   nbblok2d = 0;
        pastix_int_t   bcscidx  = 0; /* Index of the local classic cblk */
        pastix_int_t   sndeidx  = 0; /* Index of the current supernode in the original elimination tree */
        pastix_int_t   cblknum  = 0;
        pastix_int_t   brownum  = 0;
        pastix_int_t   coefnbr  = 0;
        pastix_int_t   nodenbr  = 0;

        for ( i = 0; i < symbmtx->cblknbr; i++, symbcblk++, candcblk++ ) {
            SolverBlok  *fblokptr;
            pastix_int_t nbcols;
            int recvcnt = 0;
            int tasks2D, flaglocal;

            flaglocal = ( simuctrl->cblktab[i].owned ) || ( ftgttab[i] != NULL );
            if ( !flaglocal ) {
                continue;
            }

            tasks2D = candcblk->cblktype & CBLK_TASKS_2D;

            if ( simuctrl->cblktab[i].owned ) {
                pastix_int_t cblksize, fcolnum, lcolnum;

                symbol_cblk_get_colnum( symbmtx, symbcblk, &fcolnum, &lcolnum );

                /*
                 * Register reception cblk
                 */
                ftgtcblk = ftgttab[i];
                while( ftgtcblk != NULL ) {
                    assert( (ftgtcblk->ownerid != -1) &&
                            (ftgtcblk->ownerid != ctrl->clustnum) );

                    /*
                     * solvcblk->cblktype set to CBLK_RECV temporarily for the blok
                     * faninnm initialisation.
                     */
                    solvcblk->cblktype = CBLK_RECV;
                    solvblok = solvMatGen_register_remote_cblk( solvmtx, symbmtx, ftgtcblk,
                                                                candcblk, cblklocalnum,
                                                                solvcblk, solvblok,
                                                                cblknum, brownum, i, faninnbr_tab );

                    /* Initialize missing fields and set to RECV */
                    solvcblk->brown2d   = brownum;
                    solvcblk->cblktype |= CBLK_RECV;
                    solvcblk->priority  = -1;

                    /* Update colmax is necessary */
                    solvmtx->colmax = pastix_imax( solvmtx->colmax, cblk_colnbr(solvcblk) );

                    /* Update the maximum reception buffer size */
                    cblksize        = cblk_colnbr(solvcblk) * solvcblk->stride;

                    if ( solvcblk->cblktype & CBLK_COMPRESSED ) {
                        cblksize += (solvblok - solvcblk->fblokptr);
                    }
                    solvmtx->maxrecv = pastix_imax( solvmtx->maxrecv, cblksize );

                    /* Update information about 1d/2d tasks */
                    solvMatGen_cblkIs2D( solvmtx, &nbcblk2d, &nbblok2d,
                                         (solvblok - solvcblk->fblokptr),
                                         tasks2D, cblknum );

                    /*
                     * lcolidx of the RECV block must be shifted is the
                     * reception does not start at the firs column
                     */
                    solvcblk->lcolidx = nodenbr + solvcblk->fcolnum - fcolnum;

                    recvcnt++;
                    solvmtx->recvcnt++;
                    cblknum++;

                    solvcblk->bcscnum = -( solvmtx->fanincnt + solvmtx->recvcnt );
                    solvcblk->gfaninnum = faninnbr_tab[2 * (clustnbr + ftgtcblk->ownerid)] + solvmtx->gcblknbr;
                    faninnbr_tab[2 * (clustnbr + ftgtcblk->ownerid)]++;
                    solvcblk++;

                    {
                        solver_cblk_recv_t *current = ftgtcblk;
                        ftgtcblk = ftgtcblk->next;
                        free( current );
                    }
                }
                ftgttab[i] = NULL;

                /*
                 * Register the local cblk
                 */
                solvblok = solvMatGen_register_local_cblk( symbmtx, candcblk, cblklocalnum,
                                                           solvcblk, solvblok,
                                                           cblknum, brownum, i, ctrl->clustnum );

                /* Store the information for the bcsc */
                solvcblk->bcscnum = bcscidx;
                bcscidx++;

                /* Store index for the RHS */
                solvcblk->lcolidx = nodenbr;
                assert( ((clustnbr >  1) && (nodenbr <= solvcblk->fcolnum)) ||
                        ((clustnbr == 1) && (nodenbr == solvcblk->fcolnum)) );

                /* Update local statistics */
                nbcols = cblk_colnbr( solvcblk );
                nodenbr += nbcols;
                coefnbr += solvcblk->stride * nbcols;
            }
            /*
             * Register a fanin
             */
            else {
                ftgtcblk = ftgttab[i];

                assert( ftgtcblk != NULL );
                assert( (ftgtcblk->ownerid != -1) &&
                        (ftgtcblk->ownerid != ctrl->clustnum) );

                /*
                    * solvcblk->cblktype set to CBLK_FANIN temporarily for the blok
                    * faninnm initialisation.
                    */
                solvcblk->cblktype = CBLK_FANIN;
                solvblok = solvMatGen_register_remote_cblk( solvmtx, symbmtx, ftgtcblk,
                                                            candcblk, cblklocalnum,
                                                            solvcblk, solvblok,
                                                            cblknum, brownum, i, faninnbr_tab );

                /* Set to fan-in */
                solvcblk->cblktype |= CBLK_FANIN;
                solvcblk->priority  = -1;
                solvmtx->fanincnt++;
                solvcblk->bcscnum = -( solvmtx->fanincnt + solvmtx->recvcnt );
                solvcblk->gfaninnum = faninnbr_tab[2 * ftgtcblk->ownerid] + solvmtx->gcblknbr;
                faninnbr_tab[2 * ftgtcblk->ownerid]++;

                nbcols = cblk_colnbr( solvcblk );

                /* Update colmax is necessary */
                solvmtx->colmax = pastix_imax( solvmtx->colmax, nbcols );

                free( ftgtcblk );
                ftgttab[i] = NULL;
            }

            /* Update the 1D/2D infos of the solvmtx through a cblk. */
            solvMatGen_cblkIs2D( solvmtx, &nbcblk2d, &nbblok2d,
                                 (solvblok - solvcblk->fblokptr), tasks2D, cblknum );

#if defined(PASTIX_WITH_MPI)
            if ( clustnbr > 1 ) {
#if defined(PASTIX_BLEND_FANIN_FR)
                if ( ( solvcblk->cblktype & CBLK_COMPRESSED ) &&
                     ( solvcblk->cblktype & ( CBLK_FANIN | CBLK_RECV ) ) ) {
                    solvcblk->cblktype &= ( ~CBLK_COMPRESSED );
                }
#endif
                /* No Schur complement in distributed for the moment */
                if( solvcblk->cblktype & CBLK_IN_SCHUR ) {
                    static int warning_schur = 1;
                    if ( warning_schur && (solvmtx->clustnum == 0) ) {
                        fprintf( stderr,
                                 "Warning: Schur complement support is not yet available with multiple MPI processes\n"
                                 "         It is thus disabled and the factorization will be fully performed\n" );
                        warning_schur = 0;
                    }
                    solvcblk->cblktype &= (~CBLK_IN_SCHUR);
                }
            }
#endif

            /* Store first local cblk in Schur */
            if ( (cblknum < solvmtx->cblkschur) &&
                 (solvcblk->cblktype & CBLK_IN_SCHUR) )
            {
                solvmtx->cblkschur = cblknum;
            }

            solvmtx->gcbl2loc[i] = cblknum;
            assert( cblknum == (solvcblk - solvmtx->cblktab) );

            /* Compute the original supernode in the nested dissection */
            sndeidx = solvMatGen_supernode_index( symbcblk, solvcblk, sndeidx, ordeptr );

            /*
             * Copy browtab information
             * In case of 2D tasks, we reorder the browtab to first store
             * the 1D contributions, and then the 2D updates.
             * This might also be used for low rank compression, to first
             * accumulate small dense contributions, and then, switch to a
             * low rank - low rank update scheme.
             */
            {
                pastix_int_t brownbr;
                brownbr = solvMatGen_reorder_browtab( symbmtx, symbcblk, solvmtx, solvcblk,
                                                      browtmp, cblklocalnum, bloklocalnum, brownum );

                /* Diagonal bloks of CBLK_RECV are added at the end of the browtab */
                while( recvcnt ) {
                    fblokptr = solvcblk[-recvcnt].fblokptr;
                    solvmtx->browtab[brownum + brownbr] = fblokptr - solvmtx->bloktab;
                    fblokptr->browind = brownum + brownbr;
                    brownbr++;

                    /* Supernode index is copied in too */
                    solvcblk[-recvcnt].sndeidx = solvcblk->sndeidx;
                    recvcnt--;
                }

                brownum += brownbr;

                assert( brownum <= solvmtx->brownbr );
                assert( solvcblk->brown2d >= solvcblk->brownum );
                assert( solvcblk->brown2d <= solvcblk->brownum + brownbr );
            }

            cblknum++;
            solvcblk++;
        }

        /*  Add a virtual cblk to avoid side effect in the loops on cblk bloks */
        if ( cblknum > 0 ) {
            solvMatGen_init_cblk( solvcblk, solvblok, candcblk, symbcblk,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  brownum, 0, -1, solvmtx->clustnum );
            solvcblk->lcolidx = nodenbr;
        }

        /*  Add a virtual blok to avoid side effect in the loops on cblk bloks */
        if ( solvmtx->bloknbr > 0 ) {
            solvMatGen_init_blok( solvblok, symbmtx->cblknbr + 1, symbmtx->cblknbr + 1,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  0, 0, 0 );
        }

        solvmtx->nodenbr = nodenbr;
        solvmtx->coefnbr = coefnbr;

        solvmtx->nb2dcblk = nbcblk2d;
        solvmtx->nb2dblok = nbblok2d;

        assert( solvmtx->cblkmax1d + 1 >= solvmtx->cblkmin2d );
        assert( solvmtx->brownbr  == brownum );
        assert( solvmtx->cblknbr  == cblknum );
        assert( solvmtx->faninnbr == solvmtx->fanincnt );
        assert( solvmtx->recvnbr  == solvmtx->recvcnt );
        assert( solvmtx->bloknbr  == solvblok - solvmtx->bloktab );
    }
    memFree_null( browtmp );

    /* Fill in tasktab */
    solvMatGen_fill_tasktab( solvmtx, isched, simuctrl,
                             tasklocalnum, cblklocalnum,
                             bloklocalnum, ctrl->clustnum, 0 );

    memFree_null(cblklocalnum);
    memFree_null(bloklocalnum);
    memFree_null(tasklocalnum);
    memFree_null(ftgttab);

    if ( clustnbr > 1 ) {
        memFree_null( faninnbr_tab );
    }

    /* Compute the maximum area of the temporary buffer */
    solvMatGen_max_buffers( solvmtx );
    solvMatGen_stats_last( solvmtx );

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Initialize the solver matrix structure in sequential
 *
 * This function takes all the global preprocessing steps: the symbol matrix,
 * and the result of the simulation step to generate the solver matrix for one
 * PaStiX process.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          On entry, the allocated pointer to a solver matrix structure.
 *          On exit, this structure holds alls the local information required to
 *          perform the numerical factorization.
 *
 * @param[in] symbmtx
 *          The global symbol matrix structure.
 *
 * @param[in] ordeptr
 *          The ordering structure.
 *
 * @param[in] simuctrl
 *          The information resulting from the simulation that will provide the
 *          data mapping, and the order of the task execution for the static
 *          scheduling.
 *
 * @param[in] ctrl
 *          The blend control structure that contains extra information
 *          computed during the analyze steps and the parameters of the analyze
 *          step.
 *
 * @param[in] comm
 *          TODO
 *
 * @param[in] isched
 *          TODO
 *
 * @param[in] is_dbg
 *          TODO
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
int
solverMatrixGenSeq( SolverMatrix          *solvmtx,
                    const symbol_matrix_t *symbmtx,
                    const pastix_order_t  *ordeptr,
                    const SimuCtrl        *simuctrl,
                    const BlendCtrl       *ctrl,
                    PASTIX_Comm            comm,
                    isched_t              *isched,
                    pastix_int_t           is_dbg )
{
    pastix_int_t  i;
    pastix_int_t *browtmp = 0;
    (void)ordeptr;

    assert( symbmtx->baseval == 0 );

    solverInit( solvmtx );

    solvmtx->clustnum  = ctrl->clustnum;
    solvmtx->clustnbr  = ctrl->clustnbr;
    solvmtx->procnbr   = ctrl->total_nbcores;
    solvmtx->thrdnbr   = ctrl->local_nbthrds;
    solvmtx->bublnbr   = ctrl->local_nbctxts;
    solvmtx->solv_comm = comm;

    /* Set values computed through solvMatGen_fill_localnum in distributed */
    solvmtx->tasknbr = simuctrl->tasknbr;
    solvmtx->cblknbr = symbmtx->cblknbr;
    solvmtx->bloknbr = symbmtx->bloknbr;
    solvmtx->brownbr = symbmtx->cblktab[ solvmtx->cblknbr ].brownum
                     - symbmtx->cblktab[0].brownum;

    solvmtx->cblkmin2d = solvmtx->cblknbr;
    solvmtx->cblkschur = solvmtx->cblknbr;
    solvmtx->gcblknbr  = symbmtx->cblknbr;

    /***************************************************************************
     * Fill in the local bloktab and cblktab
     */
    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN(solvmtx->cblktab, solvmtx->cblknbr+1, SolverCblk  );
    MALLOC_INTERN(solvmtx->bloktab, solvmtx->bloknbr+1, SolverBlok  );
    MALLOC_INTERN(solvmtx->browtab, solvmtx->brownbr,   pastix_int_t);
    MALLOC_INTERN(browtmp,          symbmtx->browmax,   pastix_int_t);
    {
        SolverCblk    *solvcblk = solvmtx->cblktab;
        SolverBlok    *solvblok = solvmtx->bloktab;
        symbol_cblk_t *symbcblk = symbmtx->cblktab;
        Cand          *candcblk = ctrl->candtab;
        pastix_int_t   nbcblk2d = 0;
        pastix_int_t   nbblok2d = 0;
        pastix_int_t   sndeidx  = 0; /* Index of the current supernode in the original elimination tree */
        pastix_int_t   cblknum  = 0;
        pastix_int_t   brownum  = 0;
        pastix_int_t   coefnbr  = 0;
        pastix_int_t   nodenbr  = 0;

        for(i=0; i<symbmtx->cblknbr; i++, symbcblk++, candcblk++)
        {
            pastix_int_t nbcols;
            int tasks2D = candcblk->cblktype & CBLK_TASKS_2D;

            /*
             * Register the local cblk
             */
            solvblok = solvMatGen_register_local_cblk( symbmtx, candcblk, NULL,
                                                       solvcblk, solvblok,
                                                       cblknum, brownum, i,
                                                       simuctrl->bloktab[ symbcblk->bloknum ].ownerclust );

            /* Store the information for the bcsc */
            solvcblk->bcscnum = i;
            solvcblk->lcolidx = solvcblk->fcolnum;

            /* Update local statistics */
            assert( nodenbr == solvcblk->fcolnum );
            nbcols = cblk_colnbr( solvcblk );
            nodenbr += nbcols;
            coefnbr += solvcblk->stride * nbcols;

            solvMatGen_cblkIs2D( solvmtx, &nbcblk2d, &nbblok2d,
                                 (solvblok - solvcblk->fblokptr), tasks2D, cblknum );

#if defined(PASTIX_WITH_MPI)
            if ( (solvcblk->cblktype & CBLK_COMPRESSED) &&
                 (solvcblk->cblktype & (CBLK_FANIN | CBLK_RECV)) )
            {
                solvcblk->cblktype &= (~CBLK_COMPRESSED);
            }
#endif

            /* Store first local cblk in Schur */
            if ( (cblknum < solvmtx->cblkschur) &&
                 (solvcblk->cblktype & CBLK_IN_SCHUR) )
            {
                solvmtx->cblkschur = cblknum;
            }

            /* Compute the original supernode in the nested dissection */
            sndeidx = solvMatGen_supernode_index( symbcblk, solvcblk, sndeidx, ordeptr );

            /*
             * Copy browtab information
             * In case of 2D tasks, we reorder the browtab to first store
             * the 1D contributions, and then the 2D updates.
             * This might also be used for low rank compression, to first
             * accumulate small dense contributions, and then, switch to a
             * low rank - low rank update scheme.
             */
            {
                pastix_int_t brownbr;
                brownbr = solvMatGen_reorder_browtab( symbmtx, symbcblk, solvmtx, solvcblk,
                                                      browtmp, NULL, NULL, brownum );

                brownum += brownbr;

                assert( brownum <= solvmtx->brownbr );
                assert( solvcblk->brown2d >= solvcblk->brownum );
                assert( solvcblk->brown2d <= solvcblk->brownum + brownbr );
            }
            cblknum++;
            solvcblk++;
        }

        /*  Add a virtual cblk to avoid side effect in the loops on cblk bloks */
        if ( cblknum > 0 ) {
            solvMatGen_init_cblk( solvcblk, solvblok, candcblk, symbcblk,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  symbcblk->brownum, 0, -1, ctrl->clustnum);
            solvcblk->lcolidx = solvcblk->fcolnum;
        }

        /*  Add a virtual blok to avoid side effect in the loops on cblk bloks */
        if ( solvmtx->bloknbr > 0 ) {
            solvMatGen_init_blok( solvblok, symbmtx->cblknbr + 1, symbmtx->cblknbr + 1,
                                  solvcblk[-1].lcolnum + 1, solvcblk[-1].lcolnum + 1,
                                  0, 0, 0 );
        }

        solvmtx->nodenbr = nodenbr;
        solvmtx->coefnbr = coefnbr;

        solvmtx->nb2dcblk = nbcblk2d;
        solvmtx->nb2dblok = nbblok2d;

        assert( solvmtx->cblkmax1d+1 >= solvmtx->cblkmin2d );
        assert( solvmtx->cblknbr == cblknum );
        assert( solvmtx->bloknbr == solvblok - solvmtx->bloktab );
    }
    memFree_null( browtmp );

    /*
     * Update browind fields
     */
    for(i=0; i<solvmtx->brownbr; i++)
    {
        pastix_int_t bloknum = solvmtx->browtab[i];
        solvmtx->bloktab[ bloknum ].browind = i;
    }

    /* Fill in tasktab */
    solvMatGen_fill_tasktab( solvmtx, isched, simuctrl,
                             NULL, NULL, NULL, ctrl->clustnum, is_dbg );

    /* Compute the maximum area of the temporary buffer */
    solvMatGen_max_buffers( solvmtx );
    solvMatGen_stats_last( solvmtx );

    return PASTIX_SUCCESS;
}
