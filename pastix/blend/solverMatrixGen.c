#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>

#include "common.h"
#include "dof.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "extendVector.h"
#include "elimin.h"
#include "cand.h"
#include "blendctrl.h"
#include "simu.h"
#include "solver_check.h"
#include "task.h"
#include "fanboth2.h"
#include "solverRealloc.h"
#include "solver_io.h"
#include "solverMatrixGen.h"

void build_smx(UpDownVector          *updovct,
               const SymbolMatrix    *symbptr,
               const SimuCtrl        *simuptr,
               const BlendCtrl *const ctrl,
               const Dof       *const dofptr)
{
    pastix_int_t i, j;
    pastix_int_t cursor = 0;
    pastix_int_t xcolnbr = 0;
    pastix_int_t Xnbr = 0;
    pastix_int_t localXnbr = 0;
    pastix_int_t delta;

    for(i=0;i<symbptr->cblknbr;i++)
    {
#ifdef DOF_CONSTANT
        delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;
        /*if(cbprtab[i] == ctrl->procnum)*/
        if(simuptr->bloktab[symbptr->cblktab[i].bloknum].ownerclust == ctrl->clustnum)
            localXnbr += delta;
        Xnbr += delta;
#else
        EXIT(MOD_BLEND,INTERNAL_ERR);
#endif
    }

    /** We build the whole second member **/
    for(i=0;i<symbptr->cblknbr;i++)
    {
        /** Compute xcolnbr **/
        delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;


        /** Add delta in the ODB variable **/
        xcolnbr = 0;
        for(j=symbptr->cblktab[i].bloknum+1;j<symbptr->cblktab[i+1].bloknum;j++)
        {
            xcolnbr += (symbptr->bloktab[j].lrownum - symbptr->bloktab[j].frownum + 1)*dofptr->noddval;
        }
        /** We only count non diagonal terms **/
        /** Now add the height of the cblk (-1 for diagonal term) in the DIAGONAL variables **/
        xcolnbr += delta-1;
    }

    /** Now we fill the local second member **/
    updovct->sm2xsze = localXnbr;
    updovct->sm2xnbr = 1;
    updovct->sm2xtab = NULL;

    /* Find the sm2xmax = cblk the broadest on other proc */
    updovct->sm2xmax = 0;
    for(i=0;i<symbptr->cblknbr;i++)
    {
        delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;
        if(updovct->sm2xmax < delta)
            updovct->sm2xmax = delta;
    }

    j = 0;
    for(i=0;i<symbptr->cblknbr;i++)
        /*if(cbprtab[i] == ctrl->procnum)*/
        if(simuptr->bloktab[symbptr->cblktab[i].bloknum].ownerclust == ctrl->clustnum)
        {
            delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;

            updovct->cblktab[j].sm2xind = cursor;
            j++;
            cursor += delta;
        }
}



pastix_int_t *
solverMatrixGen(const pastix_int_t clustnum,
                              SolverMatrix *solvmtx,
                              const SymbolMatrix *symbmtx,
                              const SimuCtrl * simuctrl,
                              const BlendCtrl * ctrl,
                              const Dof * dofptr)
{
    pastix_int_t            p, c;
    pastix_int_t            cursor, cursor2;
    pastix_int_t            flag;
    pastix_int_t            i, j, k;
    pastix_int_t            ftgtnum          = 0;
    pastix_int_t            coefnbr          = 0;
    pastix_int_t            nodenbr          = 0;
    pastix_int_t            odb_nbr          = 0;
    pastix_int_t            cblknum          = 0;
    pastix_int_t            bloknum          = 0;
    pastix_int_t            tasknum          = 0;
    pastix_int_t            indnbr           = 0;
    pastix_int_t          * cblklocalnum     = NULL;
    pastix_int_t          * bloklocalnum     = NULL;
    pastix_int_t          * tasklocalnum     = NULL;
    pastix_int_t          * ftgtlocalnum     = NULL;
    pastix_int_t          * bcofind          = NULL;
    pastix_int_t          * clust_mask       = NULL;
    pastix_int_t          * clust_first_cblk = NULL;
    pastix_int_t          * clust_highest    = NULL;
    pastix_int_t          * uprecvcblk       = NULL;
    pastix_int_t            flaglocal        = 0;

#ifdef PASTIX_DYNSCHED
    solvmtx->btree    = ctrl->btree;
#endif
    solvmtx->clustnum = ctrl->clustnum;
    solvmtx->clustnbr = ctrl->clustnbr;
    solvmtx->procnbr  = ctrl->total_nbcores;
    solvmtx->thrdnbr  = ctrl->local_nbthrds;
    solvmtx->bublnbr  = ctrl->local_nbctxts;
    solvmtx->ftgtcnt  = simuctrl->ftgtcnt;
#ifdef STARPU_GET_TASK_CTX
    solvmtx->starpu_subtree_nbr = symbmtx->starpu_subtree_nbr;
#endif

    if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
    {
        fprintf(stdout, "NUMBER of THREAD %ld \n", (long) solvmtx->thrdnbr );
        fprintf(stdout, "NUMBER of BUBBLE %ld \n", (long) solvmtx->bublnbr );
    }

    /* Copy the vector used to get a cluster number from a processor number */
    MALLOC_INTERN(solvmtx->proc2clust, solvmtx->procnbr, pastix_int_t);
    memcpy(solvmtx->proc2clust, ctrl->core2clust, sizeof(pastix_int_t)*solvmtx->procnbr);

    /** Be sure initialized **/
    solvmtx->coefmax = 0;

    /***************************************************************************
     * Compute local indices to compress the symbol information into solver
     */
    {
        pastix_int_t *localindex;

        MALLOC_INTERN(localindex, ctrl->clustnbr, pastix_int_t);
        memset( localindex, 0, ctrl->clustnbr * sizeof(pastix_int_t) );

        /* Compute local number of tasks on each cluster */
        MALLOC_INTERN(tasklocalnum, simuctrl->tasknbr, pastix_int_t);
        for(i=0; i<simuctrl->tasknbr; i++) {
            c = simuctrl->bloktab[simuctrl->tasktab[i].bloknum].ownerclust;

            tasklocalnum[i] = localindex[c];
            localindex[c]++;
        }
        solvmtx->tasknbr = localindex[clustnum];

        /* Compute the local numbering of the bloks and cblks on each cluster */
        MALLOC_INTERN(bloklocalnum, symbmtx->bloknbr, pastix_int_t);
        MALLOC_INTERN(cblklocalnum, symbmtx->cblknbr, pastix_int_t);

        memset( localindex, 0, ctrl->clustnbr * sizeof(pastix_int_t) );
        cblknum = 0;
        for(i=0; i<symbmtx->cblknbr; i++)
        {
            flaglocal = 0;
            for(j = symbmtx->cblktab[i].bloknum; j<symbmtx->cblktab[i+1].bloknum; j++)
            {
                c = simuctrl->bloktab[j].ownerclust;
                bloklocalnum[j] = localindex[c];
                localindex[c]++;

                if (c == clustnum)
                    flaglocal = 1;
            }

            if(flaglocal) {
                cblklocalnum[i] = cblknum;
                cblknum++;
            }
            else {
                cblklocalnum[i] = -1;
            }
        }
        solvmtx->bloknbr = localindex[clustnum];
        solvmtx->cblknbr = cblknum;

        memFree_null(localindex);
    }

    /***************************************************************************
     * Fill in bloktab and cblktab
     */

    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN(solvmtx->cblktab, solvmtx->cblknbr+1, SolverCblk);
    MALLOC_INTERN(solvmtx->bloktab, solvmtx->bloknbr,   SolverBlok);
    {
        SolverCblk *solvcblk = solvmtx->cblktab;
        SolverBlok *solvblok = solvmtx->bloktab;
        SymbolCblk *symbcblk = symbmtx->cblktab;
        SymbolBlok *symbblok = symbmtx->bloktab;
        SimuBlok   *simublok = simuctrl->bloktab;
        pastix_int_t blokamax = 0; /* Maximum area of a block in the global matrix */

        cblknum = 0;
        bloknum = 0;
        nodenbr = 0;
        coefnbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++, symbcblk++)
        {
            pastix_int_t fbloknum  = symbcblk[0].bloknum;
            pastix_int_t lbloknum  = symbcblk[1].bloknum;
            pastix_int_t stride    = 0;
            pastix_int_t nbcolumns = (symbcblk->lcolnum - symbcblk->fcolnum + 1) * dofptr->noddval;
            pastix_int_t nbrows;

            flaglocal = 0;
            cursor = bloknum;

            for( j=fbloknum; j<lbloknum; j++, symbblok++, simublok++ ) {
                nbrows = (symbblok->lrownum - symbblok->frownum + 1) * dofptr->noddval;

                blokamax = pastix_imax( blokamax, nbrows * nbcolumns );

                if(simublok->ownerclust == clustnum)
                {
                    flaglocal = 1;

                    /* Init the blok */
                    solvblok->frownum = symbblok->frownum * dofptr->noddval;
                    solvblok->lrownum = solvblok->frownum + nbrows - 1;
                    solvblok->cblknum = cblklocalnum[symbblok->cblknum];
                    //solvblok->levfval;
                    solvblok->coefind = stride;

                    stride += nbrows;
                    bloknum ++; solvblok++;
                }
            }
            if(flaglocal)
            {
                /* Init the cblk */
                solvcblk->fcolnum  = symbcblk->fcolnum * dofptr->noddval;
                solvcblk->lcolnum  = solvcblk->fcolnum + nbcolumns - 1;
                solvcblk->bloknum  = cursor;
                solvcblk->stride   = stride;
                solvcblk->procdiag = -1;
                solvcblk->coeftab  = NULL;
                solvcblk->ucoeftab = NULL;

                /* Extra statistic informations */
                nodenbr += nbcolumns;
                coefnbr += stride * nbcolumns;

                cblknum++; solvcblk++;
            }
        }

        /*  Add a virtual cblk to avoid side effect in the loops on cblk bloks */
        if (cblknum > 0)
        {
            solvcblk->fcolnum  = solvcblk->lcolnum + 1;
            solvcblk->lcolnum  = solvcblk->lcolnum + 1;
            solvcblk->bloknum  = bloknum;
            solvcblk->stride   = 0;
            solvcblk->procdiag = -1;
            solvcblk->coeftab  = NULL;
            solvcblk->ucoeftab = NULL;
        }

        solvmtx->nodenbr = nodenbr;
        solvmtx->coefnbr = coefnbr;
        solvmtx->arftmax = blokamax;

        assert( solvmtx->cblknbr == cblknum );
        assert( solvmtx->bloknbr == bloknum );
    }

    /***************************************************************************
     * Fill in tasktab
     */
    MALLOC_INTERN(solvmtx->tasktab, solvmtx->tasknbr+1, Task);
    {
        SimuTask    *simutask = simuctrl->tasktab;
        Task        *solvtask = solvmtx->tasktab;
        pastix_int_t nbftmax  = 0;

        tasknum = 0;
        ftgtnum = 0;
        indnbr  = 0;

        for(i=0; i<simuctrl->tasknbr; i++, simutask++)
        {
            if( simuctrl->bloktab[ simutask->bloknum ].ownerclust == clustnum )
            {
                assert( tasknum == tasklocalnum[i] );

                solvtask->taskid  = simutask->taskid;
                solvtask->prionum = simutask->prionum;
                solvtask->cblknum = cblklocalnum[ simutask->cblknum ];
                solvtask->bloknum = bloklocalnum[ simutask->bloknum ];
                solvtask->ftgtcnt = simutask->ftgtcnt;
                solvtask->ctrbcnt = simutask->ctrbcnt;
                solvtask->indnum  = indnbr;

                nbftmax = pastix_imax( nbftmax, solvtask->ftgtcnt );

                /*
                 * Count number of index needed in indtab:
                 *  => number of off-diagonal block below the block (including the block itself)
                 */
                odb_nbr = symbmtx->cblktab[ simutask->cblknum + 1 ].bloknum - simutask->bloknum - 1;

                switch(solvtask->taskid)
                {
                case COMP_1D:
                    indnbr += (odb_nbr*(odb_nbr+1))/2;
                    break;
                default:
                    fprintf(stderr, "solverMatrixgen: Error no task type \n");
                    EXIT(MOD_BLEND,INTERNAL_ERR);
                }

                tasknum++; solvtask++;
            }
        }
        assert(tasknum == solvmtx->tasknbr);

        /* One more task to avoid side effect */
        solvtask->taskid  = -1;
        solvtask->prionum = -1;
        solvtask->cblknum = solvmtx->cblknbr+1;
        solvtask->bloknum = solvmtx->bloknbr+1;
        solvtask->ftgtcnt = 0;
        solvtask->ctrbcnt = 0;
        solvtask->indnum  = indnbr;

        /* Store the final indnbr */
        solvmtx->indnbr  = indnbr;
        solvmtx->nbftmax = nbftmax;
    }

    /***************************************************************************
     * Fill in the ttsktab arrays (one per thread)
     *
     * TODO: This would definitely be better if each thread was initializing
     * it's own list on its own memory node.
     */
    {
        SimuProc *simuproc = &(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum]);

        /* Number of processor in this cluster */
        k = solvmtx->bublnbr;
        MALLOC_INTERN(solvmtx->ttsknbr, k, pastix_int_t  );
        MALLOC_INTERN(solvmtx->ttsktab, k, pastix_int_t* );

        for(p=0; p<k; p++, simuproc++)
        {
            pastix_int_t priomin = INTVALMAX;
            pastix_int_t priomax = 0;
            pastix_int_t ttsknbr = extendint_Size( simuproc->tasktab );
            pastix_int_t j, jloc;

            solvmtx->ttsknbr[p] = ttsknbr;

            if(ttsknbr > 0) {
                MALLOC_INTERN(solvmtx->ttsktab[p], ttsknbr, pastix_int_t);
            }
            else {
                solvmtx->ttsktab[p] = NULL;
            }

            for(i=0; i<ttsknbr; i++)
            {
                j    = extendint_Read(simuproc->tasktab, i);
                jloc = tasklocalnum[j];
                solvmtx->ttsktab[p][i] = jloc;

#if (defined PASTIX_DYNSCHED) || (defined TRACE_SOPALIN)
                solvmtx->tasktab[jloc].threadid = p;
#endif
                priomax = pastix_imax( solvmtx->tasktab[jloc].prionum, priomax );
                priomin = pastix_imin( solvmtx->tasktab[jloc].prionum, priomin );
            }

#ifdef PASTIX_DYNSCHED
            solvmtx->btree->nodetab[p].priomin = priomin;
            solvmtx->btree->nodetab[p].priomax = priomax;
#endif
        }
    }

    /***************************************************************************
     * Fill in ftgttab
     */
    {
        solvmtx->ftgtnbr = 0;
        solvmtx->ftgttab = NULL;

        /* Compute local number of outgoing contributions */
        for(c=0; c<ctrl->clustnbr; c++) {
            if(c == clustnum) {
                assert( extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c])) == 0 );
                continue;
            }
            solvmtx->ftgtnbr += extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c]));
        }

        if(solvmtx->ftgtnbr > 0) {
            SimuCluster *simuclust = &(simuctrl->clustab[clustnum]);
            FanInTarget *solvftgt;
            pastix_int_t ftgtnbr;

            MALLOC_INTERN(solvmtx->ftgttab, solvmtx->ftgtnbr, FanInTarget);

            /* Allocate array to store local indices */
            ftgtnbr = simuctrl->bloktab[symbmtx->bloknbr].ftgtnum;
            MALLOC_INTERN(ftgtlocalnum, ftgtnbr, pastix_int_t);
            memset(ftgtlocalnum, -1, ftgtnbr * sizeof(pastix_int_t));

            cursor = 0;
            solvftgt = solvmtx->ftgttab;

            for(c=0; c<ctrl->clustnbr; c++)
            {
                ftgtnbr = extendint_Size(&(simuclust->ftgtsend[c]));
                for(i=0; i<ftgtnbr; i++)
                {
                    ftgtnum = extendint_Read(&(simuclust->ftgtsend[c]), i);
                    ftgtlocalnum[ftgtnum] = cursor;

                    /* Copy information computed during simulation */
                    memcpy(solvftgt->infotab, simuctrl->ftgttab[ftgtnum].ftgt.infotab, MAXINFO*sizeof(pastix_int_t));

                    /* Update with Degre of freedoms */
                    solvftgt->infotab[FTGT_FCOLNUM] *= dofptr->noddval;
                    solvftgt->infotab[FTGT_LCOLNUM] *= dofptr->noddval;
                    solvftgt->infotab[FTGT_LCOLNUM] += dofptr->noddval - 1;
                    solvftgt->infotab[FTGT_FROWNUM] *= dofptr->noddval;
                    solvftgt->infotab[FTGT_LROWNUM] *= dofptr->noddval;
                    solvftgt->infotab[FTGT_LROWNUM] += dofptr->noddval - 1;

                    /* Convert to local numbering */
                    solvftgt->infotab[FTGT_TASKDST] = tasklocalnum[solvftgt->infotab[FTGT_TASKDST]];
                    solvftgt->infotab[FTGT_BLOKDST] = bloklocalnum[solvftgt->infotab[FTGT_BLOKDST]];

                    /* Restore ctrbcnt (modified durind simulation) */
                    solvftgt->infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[cursor].infotab[FTGT_CTRBNBR];
                    solvftgt->coeftab = NULL;

                    cursor++; solvftgt++;
                }
            }
        }
    }


    /***************************************************************************
     * Fill in indtab
     */
    {
        solvmtx->indtab = NULL;
        if (solvmtx->indnbr) {
            MALLOC_INTERN(solvmtx->indtab, solvmtx->indnbr, pastix_int_t);
        }

        indnbr = 0;
        for(i=0; i<simuctrl->tasknbr; i++)
        {
            pastix_int_t bloknum = simuctrl->tasktab[i].bloknum;
            pastix_int_t cblknum = simuctrl->tasktab[i].cblknum;

            if(simuctrl->bloktab[bloknum].ownerclust != clustnum)
                continue;

            assert(indnbr == solvmtx->tasktab[tasklocalnum[i]].indnum);
            assert(bloklocalnum[simuctrl->tasktab[i].bloknum] == solvmtx->tasktab[tasklocalnum[i]].bloknum);
            assert(cblklocalnum[simuctrl->tasktab[i].cblknum] == solvmtx->tasktab[tasklocalnum[i]].cblknum);

            switch(simuctrl->tasktab[i].taskid) {
            case COMP_1D:
            {
                pastix_int_t fbloknum = symbmtx->cblktab[cblknum  ].bloknum+1;
                pastix_int_t lbloknum = symbmtx->cblktab[cblknum+1].bloknum;

                /* For each couple (bloknum,j)\ j>=bloknum of off-diagonal block, check where goes the contribution */
                for(bloknum=fbloknum; bloknum<lbloknum; bloknum++)
                {
                    pastix_int_t firstbloknum = 0;
                    pastix_int_t facebloknum  = 0;

                    for(j=bloknum; j<lbloknum; j++)
                    {
                        facebloknum = symbolGetFacingBloknum(symbmtx, bloknum, j, firstbloknum, ctrl->ricar);

                        if(facebloknum >= 0) {
                            firstbloknum = facebloknum;

                            if(simuctrl->bloktab[facebloknum].ownerclust != clustnum)
                            {
                                solvmtx->indtab[indnbr] = ftgtlocalnum[CLUST2INDEX(facebloknum, clustnum)];
                                assert(solvmtx->indtab[indnbr] < solvmtx->ftgtnbr);
                            }
                            else
                            {
                                solvmtx->indtab[indnbr] = - tasklocalnum[simuctrl->bloktab[facebloknum].tasknum];
                                assert( (- solvmtx->indtab[indnbr]) < solvmtx->tasknbr );
                            }
                        }
                        else {
                            /* The facing block does not exist */
                            solvmtx->indtab[indnbr] =  solvmtx->ftgtnbr+1;
                        }
                        indnbr++;
                    }
                }
                break;
            }
            default:
                fprintf(stderr, "Error in solverMatrixgen: taskid unknown\n");
                EXIT(MOD_BLEND,INTERNAL_ERR);
            }
        }
        assert(indnbr == solvmtx->indnbr);
    }

    /***************************************************************************
     * Compute the maximum area of the temporary buffer used during computation
     *
     * It is either:
     *    - The panel of a diagonal block used in hetrf/sytrf factorizations of
     *      width MAXSIZEOFBLOCKS = 64
     *    - The area of the GEMM computation in a compacted update when done on
     *      CPUs
     *
     * Rk: This loop is not merged with the main block loop, since strides have
     * to be peviously computed.
     */
    {
        SolverCblk *solvcblk = solvmtx->cblktab;
        SolverBlok *solvblok = solvmtx->bloktab;
        pastix_int_t gemmmax = 0;
        pastix_int_t diagmax = 0;
        pastix_int_t gemmarea;
        pastix_int_t diagarea;

        /* Let's keep the block dimensions to print statistics informations */
        pastix_int_t maxg_m = 0;
        pastix_int_t maxg_n = 0;
        pastix_int_t maxd_m = 0;
        pastix_int_t maxd_n = 0;

        for(i=0;i<solvmtx->cblknbr;i++, solvcblk++)
        {
            pastix_int_t fbloknum = solvcblk[0].bloknum;
            pastix_int_t lbloknum = solvcblk[1].bloknum;
            pastix_int_t m = solvcblk->stride;
            pastix_int_t n = solvblok->lrownum - solvblok->frownum + 1;

            /* Temporary buffer for factorization is required only if the block
             * is larger than the blocking size */
            diagarea = n * pastix_imin( 64, pastix_imax( 0, n - 64) );
            if ( diagarea > diagmax ) {
                diagmax = diagarea;
                maxd_m = n;
                maxd_n = pastix_imin( 64, pastix_imax( 0, n - 64) );
            }

            m -= n;

            /*
             * Compute the surface of the panel for LDLt factorization
             * This could be cut down if we know at analyse time which operation
             * will be performed.
             */
            diagarea = m * n;
            if ( diagarea > diagmax ) {
                diagmax = diagarea;
                maxd_m = m;
                maxd_n = n;
            }

            solvblok++;

            /* Area of GEMM updates */
            for( j=fbloknum+1; j<lbloknum; j++, solvblok++ ) {
                n = solvblok->lrownum - solvblok->frownum + 1;

                gemmarea = m * n;
                if ( gemmarea > gemmmax ) {
                    gemmmax = gemmarea;
                    maxg_m = m;
                    maxg_n = n;
                }

                m-= n;
            }
        }

        solvmtx->coefmax = pastix_imax( gemmmax, diagmax );
        fprintf(stderr,
                "Coefmax: diagonal %ld (%ld x %ld)\n"
                "         update   %ld (%ld x %ld)\n",
                diagmax, maxd_m, maxd_n,
                gemmmax, maxg_m, maxg_n );
    }

    /****************************************/
    /** Compute the information for the    **/
    /** Forward and Backward triangular    **/
    /** Solution                           **/
    /****************************************/

    /* Pour l'instant uniquement si on est en 1d */
    if (ctrl->level2D == 0)
    {

        /** The initial symbol matrix is not expanded **/
        nodenbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
            nodenbr += (symbmtx->cblktab[i].lcolnum-symbmtx->cblktab[i].fcolnum+1)*dofptr->noddval;
        solvmtx->updovct.gnodenbr= nodenbr;
        /*fprintf(stderr," GNODENBR %ld \n", (long)solvmtx->updovct.gnodenbr);*/

        /** Build the browtabs for each diagonal block **/
        MALLOC_INTERN(solvmtx->updovct.cblktab, solvmtx->cblknbr,UpDownCblk);
        cursor = 0;
        MALLOC_INTERN(clust_mask,       ctrl->clustnbr, pastix_int_t);
        MALLOC_INTERN(clust_first_cblk, ctrl->clustnbr, pastix_int_t);
        MALLOC_INTERN(clust_highest,    ctrl->clustnbr, pastix_int_t);

        solvmtx->updovct.downmsgnbr = 0;

        for(i=0;i<symbmtx->cblknbr;i++)
        {
            pastix_int_t brownbr;

            /*if(cbprtab[i] == clustnum)*/
            if(simuctrl->bloktab[symbmtx->cblktab[i].bloknum].ownerclust == clustnum)
            {
                /*** Compute the list of clusters in the BROW (each cluster is list one time) ***/
                bzero(clust_mask, sizeof(pastix_int_t)*ctrl->clustnbr);
                bzero(clust_first_cblk, sizeof(pastix_int_t)*ctrl->clustnbr);
                for(j=0;j<ctrl->clustnbr;j++)
                    /*clust_highest[j] = - simuctrl->cblknbr;*/
                    clust_highest[j] = -1;


                brownbr = 0;
                for(j=0; j<ctrl->egraph->verttab[i].innbr;j++)
                {
                    pastix_int_t cluster;
                    pastix_int_t cblk;
                    cluster = simuctrl->bloktab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]].ownerclust;
                    cblk = ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]];
#ifdef DEBUG_M
                    ASSERT( ctrl->candtab[cblk].treelevel   <= 0,   MOD_BLEND);
                    ASSERT( ctrl->candtab[cblk].distrib     == D1,  MOD_BLEND);
                    ASSERT( simuctrl->tasktab[cblk].cblknum == cblk,MOD_BLEND);
#endif

                    /*if( ctrl->candtab[cblk].treelevel >= clust_highest[cluster] )*/
                    if( simuctrl->tasktab[cblk].prionum >= clust_highest[cluster] )
                    {
                        /*clust_highest[cluster] = ctrl->candtab[cblk].treelevel;*/
                        clust_highest[cluster] = simuctrl->tasktab[cblk].prionum;
                        clust_first_cblk[cluster] = cblk;
                    }

                    if(clust_mask[cluster] == 0)
                    {
                        clust_mask[cluster] = 1;
                        brownbr++;
                    }
                }

                solvmtx->updovct.cblktab[cursor].browprocnbr = brownbr;
                if(solvmtx->updovct.cblktab[cursor].browprocnbr>0)
                {
                    MALLOC_INTERN(solvmtx->updovct.cblktab[cursor].browproctab,
                                  solvmtx->updovct.cblktab[cursor].browprocnbr,
                                  pastix_int_t);
                }
                else
                    solvmtx->updovct.cblktab[cursor].browproctab = NULL;

                if(clust_mask[clustnum] == 1)
                    solvmtx->updovct.cblktab[cursor].msgnbr = brownbr-1;
                else
                    solvmtx->updovct.cblktab[cursor].msgnbr = brownbr;
                solvmtx->updovct.downmsgnbr += solvmtx->updovct.cblktab[cursor].msgnbr;

                /*** Alloc the vector that will contain the global cblknum with the max priority for each processor in browproctab  ***/
                if(solvmtx->updovct.cblktab[cursor].browprocnbr>0)
                {
                    MALLOC_INTERN(solvmtx->updovct.cblktab[cursor].browcblktab,
                                  solvmtx->updovct.cblktab[cursor].browprocnbr,
                                  pastix_int_t);
                }
                else
                    solvmtx->updovct.cblktab[cursor].browcblktab = NULL;

                brownbr = 0;
                for(j=0;j<ctrl->clustnbr;j++)
                    if(clust_mask[j] == 1)
                    {
                        solvmtx->updovct.cblktab[cursor].browproctab[brownbr]   = j;
                        solvmtx->updovct.cblktab[cursor].browcblktab[brownbr++] = clust_first_cblk[j];
                    }

                solvmtx->updovct.cblktab[cursor].ctrbnbr = (pastix_int_t)ctrl->egraph->verttab[i].innbr;
                cursor++;
            }
        }

        /********************************************************************/
        /*** Find the list of local blocks in front of a diagonal blocks ****/
        /********************************************************************/
        cursor  = 0;
        cursor2 = 0;
        MALLOC_INTERN(solvmtx->updovct.gcblk2list, symbmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            solvmtx->updovct.gcblk2list[i] = -1;
            flag = 0;
            for(j=0; j<ctrl->egraph->verttab[i].innbr;j++)
                if( simuctrl->bloktab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]].ownerclust == clustnum)
                {
                    if(flag == 0)
                    {
                        flag = 1;
                        solvmtx->updovct.gcblk2list[i] = cursor;
                        cursor++;
                    }
                    cursor2++;
                }
        }
        solvmtx->updovct.gcblk2listnbr = symbmtx->cblknbr;

        MALLOC_INTERN(solvmtx->updovct.listptr, cursor+1, pastix_int_t);
        solvmtx->updovct.listptrnbr = cursor+1;
        MALLOC_INTERN(solvmtx->updovct.listblok, cursor2, pastix_int_t);
        MALLOC_INTERN(solvmtx->updovct.listcblk, cursor2, pastix_int_t);

        cursor  = 0;
        cursor2 = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            flag = 0;
            for(j=0; j<ctrl->egraph->verttab[i].innbr;j++)
                if( simuctrl->bloktab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]].ownerclust == clustnum)
                {
                    if(flag == 0)
                    {
                        solvmtx->updovct.listptr[cursor2] = cursor;
                        cursor2++;
                        flag = 1;
                    }
#ifdef OOC
                    {
                        pastix_int_t tmp1,tmp2, tmp3;
                        pastix_int_t iter;
                        tmp1 = bloklocalnum[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j] ];
                        tmp2 = cblklocalnum[ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]]];
                        for (iter = solvmtx->updovct.listptr[cursor2-1]; iter < cursor; iter ++)
                            if (solvmtx->tasktab[tmp2].prionum <
                                solvmtx->tasktab[solvmtx->updovct.listcblk[iter]].prionum )
                            {
                                /* No problem with using solvmtx->updovct.listcblk[iter]
                                 * during first loop, cursor = 0, we don't use it,
                                 * and we set first.
                                 */
                                tmp3 = solvmtx->updovct.listcblk[iter];
                                solvmtx->updovct.listcblk[iter] = tmp2;
                                tmp2 = tmp3;

                                tmp3 = solvmtx->updovct.listblok[iter];
                                solvmtx->updovct.listblok[iter] = tmp1;
                                tmp1 = tmp3;
                            }
                        solvmtx->updovct.listblok[cursor] = tmp1;
                        solvmtx->updovct.listcblk[cursor] = tmp2;

                    }
#else
                    solvmtx->updovct.listblok[cursor] = bloklocalnum[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j] ];
                    solvmtx->updovct.listcblk[cursor] = cblklocalnum[ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]]];
#endif
                    cursor++;
                }

        }
        solvmtx->updovct.listptr[cursor2] = cursor;
        solvmtx->updovct.listnbr = cursor;

        solvmtx->updovct.loc2globnbr = solvmtx->cblknbr;
        MALLOC_INTERN(solvmtx->updovct.loc2glob, solvmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
            if(cblklocalnum[i] >= 0)
                solvmtx->updovct.loc2glob[cblklocalnum[i]] = i;

        memFree_null(clust_mask);
        memFree_null(clust_first_cblk);
        memFree_null(clust_highest);

        /***** Fill lblk2gcblk ******/
        solvmtx->updovct.gcblknbr = symbmtx->cblknbr;
        MALLOC_INTERN(solvmtx->updovct.lblk2gcblk, symbmtx->bloknbr, pastix_int_t);
        for(i=0;i<symbmtx->bloknbr;i++)
            if(simuctrl->bloktab[i].ownerclust == clustnum)
                solvmtx->updovct.lblk2gcblk[bloklocalnum[i]] = symbmtx->bloktab[i].cblknum;

        /* Calcul du nombre de messages a recevoir lors de la remontée */
        MALLOC_INTERN(uprecvcblk, symbmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
            uprecvcblk[i] = 0;
        for (i=0; i<solvmtx->bublnbr; i++)
            for (j=0; j < solvmtx->ttsknbr[i]; j++)
            {
                cblknum = solvmtx->tasktab[solvmtx->ttsktab[i][j]].cblknum;
                for (k =  solvmtx->cblktab[cblknum+1].bloknum-1;
                     k >= solvmtx->cblktab[cblknum].bloknum+1; k--)
                    /* if the contribution is not local */
                    if (solvmtx->bloktab[k].cblknum <= 0)
                        uprecvcblk[solvmtx->updovct.lblk2gcblk[k]] = 1;
            }
        solvmtx->updovct.upmsgnbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
            solvmtx->updovct.upmsgnbr += uprecvcblk[i];
        memFree_null(uprecvcblk);

        /*********************************/
        /*     Temporaire                */
        /*  Pour tester descente remonte */
        /*********************************/
        build_smx(&(solvmtx->updovct), symbmtx, simuctrl, ctrl, dofptr);

    }
    /*********************** END TRIANGULAR INFO BUILDING ******************************************/

    memFree_null(cblklocalnum);
    memFree_null(bloklocalnum);
    memFree_null(tasklocalnum);

    if(ftgtlocalnum != NULL)
        memFree_null(ftgtlocalnum);

    return bcofind;
}
