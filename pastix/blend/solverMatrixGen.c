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

/*#define DEBUG_PRIO*/

void working_array_boundaries(SolverMatrix *solvmtx, MPI_Comm pastix_comm);

void printSymbolMatrix(FILE *file, SymbolMatrix *symbptr)
{
    pastix_int_t i, j;
    for(i=0;i<symbptr->cblknbr;i++)
    {
        fprintf(file, "CBLK %ld [%ld : %ld ] \n",(long)i, (long)symbptr->cblktab[i].fcolnum, (long)symbptr->cblktab[i].lcolnum);
        for(j=symbptr->cblktab[i].bloknum;j<symbptr->cblktab[i+1].bloknum;j++)
            fprintf(file, "--BLOK %ld [%ld : %ld ]\n", (long)j, (long)symbptr->bloktab[j].frownum, (long)symbptr->bloktab[j].lrownum);
        fprintf(file, "\n");
    }
}

void build_smx(UpDownVector          *updovct,
               const SymbolMatrix          *symbptr,
               pastix_int_t                   *blprtab,
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
        if(ctrl->core2clust[blprtab[symbptr->cblktab[i].bloknum]] == ctrl->clustnum)
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
        if(ctrl->core2clust[blprtab[symbptr->cblktab[i].bloknum]] == ctrl->clustnum)
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
    pastix_int_t            i, j, jloc, k;
    pastix_int_t            ftgtnum          = 0;
    pastix_int_t            coefnbr          = 0;
    pastix_int_t            coefind          = 0;
    pastix_int_t            nodenbr          = 0;
    pastix_int_t            odb_nbr          = 0;
    pastix_int_t            cblknum          = 0;
    pastix_int_t            bloknum          = 0;
    pastix_int_t            tasknum          = 0;
    pastix_int_t            min              = INTVALMAX;
    pastix_int_t            facebloknum      = 0;
    pastix_int_t            delta            = 0;
    pastix_int_t            indnbr           = 0;
    pastix_int_t            btagnbr          = 0;
    pastix_int_t            bcofnbr          = 0;
    pastix_int_t          * blprtab          = simuctrl->blprtab;
    pastix_int_t          * cblklocalnum     = NULL;
    pastix_int_t          * bloklocalnum     = NULL;
    pastix_int_t          * tasklocalnum     = NULL;
    pastix_int_t          * btaglocalnum     = NULL;
    pastix_int_t          * ftgtlocalnum     = NULL;
    pastix_int_t          * bcofind          = NULL;
    pastix_int_t          * proc2clust       = NULL;
    pastix_int_t          * clust_mask       = NULL;
    pastix_int_t          * clust_first_cblk = NULL;
    pastix_int_t          * clust_highest    = NULL;
    pastix_int_t          * uprecvcblk       = NULL;
    pastix_int_t          * localindex       = NULL;
    pastix_int_t            flaglocal        = 0;
    pastix_int_t            min_procdst;
    pastix_int_t            min_task;
    pastix_int_t            maxprionum2;

    /** Set procnum and procnbr **/
    /*solvmtx->procnum = procnum;
     solvmtx->procnbr = ctrl->procnbr;*/
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

    /* Initialize pointer */
    proc2clust = ctrl->core2clust;

    /** Be sure initialized **/
    solvmtx->cpftmax = 0;
    solvmtx->bpftmax = 0;
    solvmtx->coefmax = 0;

    MALLOC_INTERN(localindex, ctrl->clustnbr, pastix_int_t);
    memset( localindex, 0, ctrl->clustnbr * sizeof(pastix_int_t) );

    /* Compute local number of tasks on each cluster */
    MALLOC_INTERN(tasklocalnum, simuctrl->tasknbr, pastix_int_t);
    for(i=0; i<simuctrl->tasknbr; i++) {
        c = proc2clust[blprtab[simuctrl->tasktab[i].bloknum]];

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
            c = proc2clust[blprtab[j]];
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

    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN(solvmtx->cblktab, solvmtx->cblknbr+1, SolverCblk);
    MALLOC_INTERN(solvmtx->bloktab, solvmtx->bloknbr,   SolverBlok);

    /* Fill in bloktab and cblktab */
    {
        cblknum = 0;
        bloknum = 0;
        nodenbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            flaglocal = 0;
            cursor = bloknum;

            for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
                if(proc2clust[blprtab[j]] == clustnum)
                {
                    flaglocal = 1;
                    solvmtx->bloktab[bloknum].frownum = symbmtx->bloktab[j].frownum * dofptr->noddval;
                    solvmtx->bloktab[bloknum].lrownum = symbmtx->bloktab[j].lrownum * dofptr->noddval + dofptr->noddval-1;
                    solvmtx->bloktab[bloknum].cblknum = cblklocalnum[symbmtx->bloktab[j].cblknum];
                    bloknum ++;
                }
            if(flaglocal)
            {
                solvmtx->cblktab[cblknum].fcolnum = symbmtx->cblktab[i].fcolnum * dofptr->noddval;
                solvmtx->cblktab[cblknum].lcolnum = symbmtx->cblktab[i].lcolnum * dofptr->noddval + dofptr->noddval-1;
                solvmtx->cblktab[cblknum].bloknum = cursor;
                solvmtx->cblktab[cblknum].coeftab = NULL;
                solvmtx->cblktab[cblknum].ucoeftab = NULL;
                nodenbr += symbmtx->cblktab[i].lcolnum - symbmtx->cblktab[i].fcolnum + 1;
                cblknum++;
            }

        }
        solvmtx->nodenbr = nodenbr;
#ifdef DEBUG_BLEND
        ASSERT(solvmtx->cblknbr == cblknum,MOD_BLEND);
        if(solvmtx->bloknbr != bloknum)
            fprintf(stderr, "bloknbr %ld bloknum %ld \n", (long)solvmtx->bloknbr, (long)bloknum);
        ASSERT(solvmtx->bloknbr == bloknum,MOD_BLEND);
#endif

        if (cblknum > 0)
        {
            /*  virtual cblk to avoid side effect in the loops on cblk bloks */
            solvmtx->cblktab[cblknum].fcolnum = solvmtx->cblktab[cblknum-1].lcolnum+1;
            solvmtx->cblktab[cblknum].lcolnum = solvmtx->cblktab[cblknum-1].lcolnum+1;
            solvmtx->cblktab[cblknum].bloknum = bloknum;
            solvmtx->cblktab[cblknum].coeftab = NULL;
            solvmtx->cblktab[cblknum].ucoeftab = NULL;
        }
    }
    /****************************/
    /*      Fill tasktab        */
    /****************************/
    MALLOC_INTERN(solvmtx->tasktab, solvmtx->tasknbr+1, Task);
    /** Initialize the tasknext field **/
    /** We need that to know if a task has already been chained **/
    for(i=0;i<solvmtx->tasknbr+1;i++)
        solvmtx->tasktab[i].tasknext = -1;


    tasknum = 0;
    ftgtnum = 0;
    indnbr  = 0;
    maxprionum2 = 0;

    for(i=0;i<simuctrl->tasknbr;i++)
    {

        maxprionum2 = 0;
        if(proc2clust[blprtab[simuctrl->tasktab[i].bloknum]] == clustnum)
        {
            ASSERTDBG(tasknum == tasklocalnum[i], MOD_BLEND);

            /*fprintf(stdout ,"i %ld indnbr %ld \n", (long)i, (long)indnbr);*/
            solvmtx->tasktab[tasknum].taskid  = simuctrl->tasktab[i].taskid;
            solvmtx->tasktab[tasknum].prionum = simuctrl->tasktab[i].prionum;
            solvmtx->tasktab[tasknum].cblknum = cblklocalnum[simuctrl->tasktab[i].cblknum];

            /**** Set the prionum for DR : the locally subtree are ordered in natural cblk order ****/
            if(ctrl->candtab[simuctrl->tasktab[i].cblknum].fcandnum == ctrl->candtab[simuctrl->tasktab[i].cblknum].lcandnum)
            {
                solvmtx->tasktab[tasknum].prionum2 = cblklocalnum[simuctrl->tasktab[i].cblknum];
                if(solvmtx->tasktab[tasknum].prionum2 > maxprionum2)
                    maxprionum2 = solvmtx->tasktab[tasknum].prionum2;
            }
            else
            {
                solvmtx->tasktab[tasknum].prionum2 = -1; /** The prionum2 for these task is set outside this loop **/
            }

            if(solvmtx->tasktab[tasknum].taskid == COMP_1D)
                solvmtx->tasktab[tasknum].bloknum = solvmtx->cblktab[solvmtx->tasktab[tasknum].cblknum].bloknum;
            else
                solvmtx->tasktab[tasknum].bloknum = bloklocalnum[simuctrl->tasktab[i].bloknum];

            solvmtx->tasktab[tasknum].ctrbcnt = simuctrl->tasktab[i].ctrbcnt;
            solvmtx->tasktab[tasknum].ftgtcnt = simuctrl->tasktab[i].ftgtcnt;
            solvmtx->tasktab[tasknum].btagptr = NULL;
            solvmtx->tasktab[tasknum].indnum  = indnbr;

            /*fprintf(stdout, "task %ld, taskid %ld cblk %ld  prionum %ld ctrbnbr %ld \n",(long)tasknum, (long)solvmtx->tasktab[tasknum].taskid, (long)solvmtx->tasktab[tasknum].cblknum, (long)solvmtx->tasktab[tasknum].prionum, (long)solvmtx->tasktab[tasknum].ctrbcnt);*/

            /** Count number of index needed in indtab **/
            /* number odb below the block (included the block) */
            odb_nbr = symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum - simuctrl->tasktab[i].bloknum -1;
#ifdef DEBUG_BLEND
            ASSERT(simuctrl->tasktab[i].taskid == solvmtx->tasktab[tasknum].taskid,MOD_BLEND);
#endif
            switch(simuctrl->tasktab[i].taskid)
            {
            case COMP_1D:
                indnbr += (odb_nbr*(odb_nbr+1))/2;
                break;
            case DIAG:
                indnbr++;
                /** Count number of BlockTarget and BlockCoeff **/

                flaglocal = 0;
                for(c  = ctrl->candtab[simuctrl->tasktab[i].cblknum].fccandnum;
                    c <= ctrl->candtab[simuctrl->tasktab[i].cblknum].lccandnum; c++)
                    for(j = simuctrl->tasktab[i].bloknum+1;
                        j < symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum; j++)
                        if(proc2clust[blprtab[j]]==c)
                        {
                            flaglocal = 1;
                            btagnbr++;
                            break;
                        }
                if(flaglocal)
                    bcofnbr++;
                break;
            case E1:
                indnbr++;
                /** Count number of BlockTarget and BlockCoeff **/
                for(c  = ctrl->candtab[simuctrl->tasktab[i].cblknum].fccandnum;
                    c <= ctrl->candtab[simuctrl->tasktab[i].cblknum].lccandnum; c++)
                    for(j = simuctrl->tasktab[i].bloknum;
                        j < symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum; j++)
                        if(proc2clust[blprtab[j]]==c)
                        {
                            btagnbr++;
                            break;
                        }
                bcofnbr++;
                break;
            case E2:
                indnbr++;
                break;
            default:
                fprintf(stderr, "solverMatrixgen: Error no task type \n");
                EXIT(MOD_BLEND,INTERNAL_ERR);
            }

            if(simuctrl->tasktab[i].taskid == E1 || simuctrl->tasktab[i].taskid == E2)
            {
                /*************************************************************/
                /** Chain (in a cycle) all the tasks that use the same btag **/
                /*************************************************************/
                if((solvmtx->tasktab[tasknum].tasknext == -1) &&  /** This task has not already been chained **/
                   (simuctrl->tasktab[i].tasknext      != -1))    /** This task is not the last diag block   **/
                {
                    ASSERTDBG(tasklocalnum[i] == tasknum, MOD_BLEND);
                    ASSERTDBG(simuctrl->tasktab[i].taskid == solvmtx->tasktab[tasknum].taskid, MOD_BLEND);

                    cursor = simuctrl->tasktab[i].tasknext;
                    while(cursor != i)
                    {
                        ASSERTDBG(simuctrl->tasktab[i].taskid == simuctrl->tasktab[cursor].taskid, MOD_BLEND);

                        /*fprintf(stdout, "i %d id %d cursor %d \n", i, simuctrl->tasktab[i].taskid, cursor);*/
                        if(proc2clust[ blprtab[simuctrl->tasktab[cursor].bloknum] ] == clustnum)
                            break;
                        cursor = simuctrl->tasktab[cursor].tasknext;

                    }
                    solvmtx->tasktab[tasknum].tasknext = tasklocalnum[cursor];
                }
            }
            else
            {
                ASSERTDBG(solvmtx->tasktab[tasknum].tasknext == -1 && simuctrl->tasktab[i].tasknext == -1,MOD_BLEND);
            }
            solvmtx->tasktab[tasknum].taskmstr = -1;
            tasknum++;
        }
    }
    /** One more to avoid side effect **/
    solvmtx->tasktab[tasknum].indnum = indnbr;

    /** For coherence in Malt (it uses takstab[i+1].tasknum to count access **/
    /*if(solvmtx->tasktab[tasknum-1].taskid == COMP_1D)
     solvmtx->tasktab[tasknum-1].indnum = indnbr-1;*/
    ASSERTDBG(tasknum == solvmtx->tasknbr,MOD_BLEND);

    /********************************************************************/
    /* Set the prionum2 for task that have several candidate thread     */
    /********************************************************************/
    for(i=0;i<solvmtx->tasknbr;i++)
        if(solvmtx->tasktab[i].prionum2 < 0)
            solvmtx->tasktab[i].prionum2 = maxprionum2 + solvmtx->tasktab[i].prionum;


    /********************************************/
    /* Fill the processor tasktab indice vector */
    /********************************************/
    /* Number of processor in this cluster */
    k = solvmtx->bublnbr;
    MALLOC_INTERN(solvmtx->ttsknbr, k, pastix_int_t);

    for(p = 0;p<k;p++)
    {
        solvmtx->ttsknbr[p] = extendint_Size(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum + p].tasktab);
        print_debug(DBG_BUBBLES, "La bulle %d execute %d taches\n", (int)p, (int)solvmtx->ttsknbr[p]);
    }

    MALLOC_INTERN(solvmtx->ttsktab, k, pastix_int_t *);
    for(p = 0;p<k;p++)
    {
#ifdef PASTIX_DYNSCHED
        pastix_int_t min = INTVALMAX;
        pastix_int_t max = 0;
#endif

        if(solvmtx->ttsknbr[p] > 0)
        {
            MALLOC_INTERN(solvmtx->ttsktab[p], solvmtx->ttsknbr[p], pastix_int_t);
        }
        else
            solvmtx->ttsktab[p] = NULL;

        for(i=0;i<solvmtx->ttsknbr[p];i++)
        {
            j    = extendint_Read(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum + p].tasktab, i);
            jloc = tasklocalnum[j];
            solvmtx->ttsktab[p][i] = jloc;

#if (defined PASTIX_DYNSCHED) || (defined TRACE_SOPALIN)
            solvmtx->tasktab[jloc].threadid = p;
#endif
#ifdef TRACE_SOPALIN
            solvmtx->tasktab[jloc].fcandnum = ctrl->candtab[simuctrl->tasktab[j].cblknum].fcandnum;
            solvmtx->tasktab[jloc].lcandnum = ctrl->candtab[simuctrl->tasktab[j].cblknum].lcandnum;
            solvmtx->tasktab[jloc].id       = simuctrl->tasktab[j].cblknum;
#endif

#ifdef PASTIX_DYNSCHED
            if ( solvmtx->tasktab[jloc].prionum > max )
                max = solvmtx->tasktab[jloc].prionum;
            if ( solvmtx->tasktab[jloc].prionum < min )
                min = solvmtx->tasktab[jloc].prionum;
#endif
        }

#ifdef PASTIX_DYNSCHED
        solvmtx->btree->nodetab[p].priomin = min;
        solvmtx->btree->nodetab[p].priomax = max;
#endif
    }

#ifdef PRINT_ORDOTASK
    {
        FILE *ordofile;
        char  ordofilename[250];
        sprintf(ordofilename, "Ordo.%02d", clustnum);
        ordofile = fopen(ordofilename, "w");
        for(p = 0;p<k;p++)
            for(i=0;i<solvmtx->ttsknbr[p];i++)
            {
                j = extendint_Read(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum + p].tasktab, i);
                fprintf(ordofile, "%ld %ld\n", j, tasklocalnum[j]);
            }
        fclose(ordofile);
    }
#endif

    /*******************/
    /** Fill ftgttab  **/
    /*******************/

    solvmtx->ftgtnbr = 0;
    for(c=0;c < ctrl->clustnbr;c++)
    {
        if(c != clustnum)
            solvmtx->ftgtnbr += extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c]));
    }

    if(solvmtx->ftgtnbr > 0)
    {
        MALLOC_INTERN(solvmtx->ftgttab, solvmtx->ftgtnbr, FanInTarget);
    }
    else
        solvmtx->ftgttab = NULL;
    MALLOC_INTERN(ftgtlocalnum, simuctrl->bloktab[symbmtx->bloknbr].ftgtnum, pastix_int_t);
    memset(ftgtlocalnum, -1, sizeof(pastix_int_t)*simuctrl->bloktab[symbmtx->bloknbr].ftgtnum);
    cursor = 0;
    for(c=0;c<ctrl->clustnbr;c++)
        for(i=0;i<extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c]));i++)
        {
            ftgtnum = extendint_Read(&(simuctrl->clustab[clustnum].ftgtsend[c]), i);
            ftgtlocalnum[ftgtnum] = cursor;
            memcpy(solvmtx->ftgttab[cursor].infotab, simuctrl->ftgttab[ftgtnum].ftgt.infotab, MAXINFO*sizeof(pastix_int_t));


#ifdef DOF_CONSTANT
            solvmtx->ftgttab[cursor].infotab[FTGT_FCOLNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LCOLNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LCOLNUM] += dofptr->noddval - 1;
            solvmtx->ftgttab[cursor].infotab[FTGT_FROWNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LROWNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LROWNUM] += dofptr->noddval - 1;
#endif

            solvmtx->ftgttab[cursor].infotab[FTGT_TASKDST] = tasklocalnum[solvmtx->ftgttab[cursor].infotab[FTGT_TASKDST]];

            solvmtx->ftgttab[cursor].infotab[FTGT_BLOKDST] = bloklocalnum[solvmtx->ftgttab[cursor].infotab[FTGT_BLOKDST]];
            /* Reinit ftgt ctrbcnt */
            solvmtx->ftgttab[cursor].infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[cursor].infotab[FTGT_CTRBNBR];
            /** Allocation for FanInTarget not assured by blend (done by sopalin)**/
            solvmtx->ftgttab[cursor].coeftab = NULL;

            /*if(p == 2)
             fprintf(stdout, "Ftgt %ld prio %ld to task %ld \n", (long)cursor, (long)solvmtx->ftgttab[cursor].infotab[FTGT_PRIONUM],
             (long)solvmtx->ftgttab[cursor].infotab[FTGT_TASKDST]);*/


            cursor++;
        }


    /*********************/
    /*    Fill indtab    */
    /*********************/
    solvmtx->btagnbr  = btagnbr;
    solvmtx->indnbr   = indnbr;
    solvmtx->bcofnbr  = bcofnbr;
    solvmtx->indtab   = NULL;
    solvmtx->btagtab  = NULL;
    bcofind           = NULL;
    solvmtx->bcoftab  = NULL;
    if (indnbr)
        MALLOC_INTERN(solvmtx->indtab, indnbr, pastix_int_t);
    if (btagnbr)
    {
        MALLOC_INTERN(solvmtx->btagtab, btagnbr, BlockTarget);
        MALLOC_INTERN(bcofind,          btagnbr, pastix_int_t);
        memset(bcofind, -1, btagnbr*sizeof(pastix_int_t));
    }
    if (bcofnbr)
        MALLOC_INTERN(solvmtx->bcoftab, bcofnbr, BlockCoeff);
    tasknum     = 0;
    indnbr      = 0;
    btagnbr     = 0;
    bcofnbr     = 0;
    for(i=0;i<simuctrl->tasknbr;i++)
    {
        if(proc2clust[blprtab[simuctrl->tasktab[i].bloknum]] == clustnum)
        {
            ASSERTDBG(tasklocalnum[i] == tasknum,MOD_BLEND);
            ASSERTDBG(indnbr == solvmtx->tasktab[tasklocalnum[i]].indnum, MOD_BLEND);
            ASSERTDBG(bloklocalnum[simuctrl->tasktab[i].bloknum] == solvmtx->tasktab[tasklocalnum[i]].bloknum,MOD_BLEND);
            ASSERTDBG(cblklocalnum[simuctrl->tasktab[i].cblknum] == solvmtx->tasktab[tasklocalnum[i]].cblknum,MOD_BLEND);

            switch(simuctrl->tasktab[i].taskid)
            {
            case COMP_1D:
                for(bloknum = symbmtx->cblktab[simuctrl->tasktab[i].cblknum].bloknum+1;
                    bloknum < symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum; bloknum++)
                {
                    facebloknum = 0;
                    for(j=bloknum;j<symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum;j++)
                    {
                        facebloknum = symbolGetFacingBloknum(symbmtx, bloknum, j, facebloknum, ctrl->ricar);
                        /*#ifdef NAPA*/
                        if(facebloknum >= 0)
                        {
                            /*#endif*/
                            if(proc2clust[blprtab[facebloknum]]!=clustnum)
                            {
                                solvmtx->indtab[indnbr] = ftgtlocalnum[CLUST2INDEX(facebloknum, clustnum)];
#ifdef DEBUG_PRIO
                                solvmtx->ftgttab[solvmtx->indtab[indnbr]].infotab[FTGT_PRIONUM]
                                    = solvmtx->tasktab[tasklocalnum[i]].prionum;
                                /*fprintf(stdout, "SOLV Task1D %ld FTGT %ld  taskprio %ld ftgtprio %ld \n", (long)tasklocalnum[i], (long)solvmtx->indtab[indnbr] , (long)solvmtx->tasktab[tasklocalnum[i]].prionum, (long)solvmtx->ftgttab[solvmtx->indtab[indnbr]].infotab[FTGT_PRIONUM]);*/
#endif

                                ASSERTDBG(solvmtx->indtab[indnbr] < solvmtx->ftgtnbr,MOD_BLEND);
                            }
                            else
                            {
                                solvmtx->indtab[indnbr] = -tasklocalnum[simuctrl->bloktab[facebloknum].tasknum];

#ifdef DEBUG_BLEND
                                if(!(-solvmtx->indtab[indnbr] < solvmtx->tasknbr))
                                    fprintf(stderr, "facetasknum %ld tasklocalnum %ld \n",
                                            (long)simuctrl->bloktab[facebloknum].tasknum, (long)-solvmtx->indtab[indnbr]);
                                ASSERT(-solvmtx->indtab[indnbr] < solvmtx->tasknbr,MOD_BLEND);
#endif
                            }
                            indnbr++;
                            ASSERTDBG(indnbr <= solvmtx->indnbr,MOD_BLEND);
                            /*#ifdef NAPA*/
                        }
                        else
                        { /** THE FACING BLOCK DO NOT EXIST **/
                            solvmtx->indtab[indnbr] =  solvmtx->ftgtnbr+1;
                            indnbr++;
                        }
                        /*#endif*/
                    }
                }
                break;
            case DIAG:
                /* We just need the index of one of the bloktarget */
                if(i == simuctrl->tasknbr - 1)
                    solvmtx->indtab[indnbr] = -1;
                else
                {
                    solvmtx->indtab[indnbr] = btagnbr;
                    solvmtx->bcoftab[bcofnbr].sendcnt = 0;
                    /*fprintf(stdout, "TASK %ld, btagnum %ld indnum %ld taskid %ld \n",(long)tasklocalnum[i], (long)btagnbr, (long)indnbr, (long)solvmtx->tasktab[tasklocalnum[i]].taskid);*/
                    flaglocal = 0;
                    for(c =  ctrl->candtab[simuctrl->tasktab[i].cblknum].fccandnum;
                        c <= ctrl->candtab[simuctrl->tasktab[i].cblknum].lccandnum; c++)
                        for(bloknum = simuctrl->tasktab[i].bloknum+1;
                            bloknum < symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum; bloknum++)
                            if(proc2clust[blprtab[bloknum]]==c)
                            {
                                flaglocal = 1;
                                /*solvmtx->btagtab[btagnbr].bcofnum  = bcofnbr;*/

                                bcofind[btagnbr] = bcofnbr;

                                /**solvmtx->btagtab[btagnbr].infotab[BTAG_TASKDST] = tasklocalnum[simuctrl->bloktab[bloknum].tasknum];**/
                                /** Find the lowest priority on the block receiver of the btag **/
                                cursor      = simuctrl->bloktab[bloknum].tasknum;
                                min         = simuctrl->tasktab[cursor].prionum;
                                min_procdst = blprtab[bloknum];
                                min_task    = tasklocalnum[cursor];

                                solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT] = 1;

                                while(simuctrl->tasktab[cursor].tasknext != simuctrl->bloktab[bloknum].tasknum)
                                {
                                    ASSERTDBG(simuctrl->tasktab[cursor].taskid == E1,MOD_BLEND);

                                    cursor = simuctrl->tasktab[cursor].tasknext;
                                    if(proc2clust[blprtab[simuctrl->tasktab[cursor].bloknum]] == c)
                                    {
                                        solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT]++;
                                        if(simuctrl->tasktab[cursor].prionum < min)
                                        {
                                            min         = simuctrl->tasktab[cursor].prionum;
                                            min_procdst = blprtab[simuctrl->tasktab[cursor].bloknum];
                                            min_task    = tasklocalnum[cursor];
                                        }
                                    }
                                }

                                /*solvmtx->btagtab[btagnbr].infotab[BTAG_PROCDST] = c;*/
                                solvmtx->btagtab[btagnbr].infotab[BTAG_TASKDST] = min_task;
                                solvmtx->btagtab[btagnbr].infotab[BTAG_PROCDST] = min_procdst;
                                solvmtx->btagtab[btagnbr].infotab[BTAG_PRIONUM] = min;

#ifdef DEBUG_BLEND
                                {
                                    pastix_int_t taskcnt;
                                    taskcnt = 1;
                                    for(j=bloknum+1; j<symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum;j++)
                                        if(proc2clust[blprtab[j]]==c)
                                            taskcnt++;
                                    ASSERT(taskcnt == solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT],MOD_BLEND);
                                }


                                {
                                    pastix_int_t tnbr;
                                    tnbr = 1;
                                    cursor = simuctrl->bloktab[bloknum].tasknum;
                                    while(simuctrl->tasktab[cursor].tasknext != simuctrl->bloktab[bloknum].tasknum)
                                    {
                                        cursor = simuctrl->tasktab[cursor].tasknext;
                                        if(proc2clust[blprtab[simuctrl->tasktab[cursor].bloknum]] == c)
                                            tnbr++;
                                    }
                                    if(tnbr != solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT])
                                        fprintf(stderr, "tnbr = %ld taskcnt = %ld \n", (long)tnbr, (long)solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT]);

                                    ASSERT(tnbr == solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT],MOD_BLEND);
                                }
#endif


                                btagnbr++;
                                solvmtx->bcoftab[bcofnbr].sendcnt++;
                                break;
                            }
                    if(flaglocal)
                    {
                        solvmtx->bcoftab[bcofnbr].infotab[BCOF_FCOLNUM] = solvmtx->bloktab[bloklocalnum[simuctrl->tasktab[i].bloknum]].frownum;
                        solvmtx->bcoftab[bcofnbr].infotab[BCOF_LCOLNUM] = solvmtx->bloktab[bloklocalnum[simuctrl->tasktab[i].bloknum]].lrownum;
                        solvmtx->bcoftab[bcofnbr].infotab[BCOF_FROWNUM] = solvmtx->bloktab[bloklocalnum[simuctrl->tasktab[i].bloknum]].frownum;
                        solvmtx->bcoftab[bcofnbr].infotab[BCOF_LROWNUM] = solvmtx->bloktab[bloklocalnum[simuctrl->tasktab[i].bloknum]].lrownum;
                        solvmtx->bcoftab[bcofnbr].coeftab = NULL;
                        bcofnbr++;
                    }
                }
                indnbr++;
                break;

            case E1:
                /* We just need the index of one of the bloktarget */
                solvmtx->indtab[indnbr] = btagnbr;
                solvmtx->bcoftab[bcofnbr].sendcnt = 0;

                /*
                 * If diag block is local taskmstr is the diag task
                 * else it's the first local E1 task
                 */
                {
                    bloknum = symbmtx->cblktab[simuctrl->tasktab[i].cblknum].bloknum;
                    if(proc2clust[blprtab[bloknum]] == clustnum)
                    {
                        solvmtx->tasktab[tasknum].taskmstr = tasklocalnum[simuctrl->bloktab[bloknum].tasknum];
                    }
                    else
                    {
                        cursor   = tasknum;
                        min      = solvmtx->tasktab[cursor].prionum;
                        min_task = tasknum;
                        while(solvmtx->tasktab[cursor].tasknext != tasknum)
                        {
                            cursor = solvmtx->tasktab[cursor].tasknext;
                            ASSERTDBG(solvmtx->tasktab[cursor].taskid == E1,MOD_BLEND);

                            if(solvmtx->tasktab[cursor].prionum < min)
                            {
                                min      = solvmtx->tasktab[cursor].prionum;
                                min_task = cursor;
                            }
                        }
                        solvmtx->tasktab[tasknum].taskmstr = min_task;
                    }
                }

                for(c=ctrl->candtab[simuctrl->tasktab[i].cblknum].fccandnum; c<=ctrl->candtab[simuctrl->tasktab[i].cblknum].lccandnum;c++)
                    for(bloknum=simuctrl->tasktab[i].bloknum;bloknum<symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum;bloknum++)
                        if(proc2clust[blprtab[bloknum]]==c)
                        {
#ifdef DEBUG_BLEND
                            pastix_int_t debug = 1;
#endif
                            bcofind[btagnbr] = bcofnbr;

                            cursor = i+bloknum-simuctrl->tasktab[i].bloknum+1;

                            ASSERTDBG(simuctrl->tasktab[cursor].bloknum == bloknum,MOD_BLEND);

                            /* solvmtx->btagtab[btagnbr].infotab[BTAG_TASKDST] = tasklocalnum[cursor];*/

                            /** Find the lowest priority on the block receiver of the btag **/
                            min_task    = tasklocalnum[cursor];
                            min         = simuctrl->tasktab[cursor].prionum;
                            min_procdst = blprtab[simuctrl->tasktab[cursor].bloknum];
                            cursor      = simuctrl->tasktab[cursor].tasknext;

                            solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT] = 1;
                            while(cursor!= i+bloknum-simuctrl->tasktab[i].bloknum+1)
                            {
                                if(proc2clust[blprtab[simuctrl->tasktab[cursor].bloknum]] == c)
                                {
                                    solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT]++;
                                    ASSERTDBG(simuctrl->tasktab[cursor].taskid == E2,MOD_BLEND);
#ifdef DEBUG_BLEND
                                    debug++;
#endif
                                    if(simuctrl->tasktab[cursor].prionum < min)
                                    {
                                        min = simuctrl->tasktab[cursor].prionum;
                                        min_procdst = blprtab[simuctrl->tasktab[cursor].bloknum];
                                        min_task = tasklocalnum[cursor];
                                    }
                                }
                                cursor = simuctrl->tasktab[cursor].tasknext;
                            }
                            /*solvmtx->btagtab[btagnbr].infotab[BTAG_PROCDST] = c;*/
                            solvmtx->btagtab[btagnbr].infotab[BTAG_TASKDST] = min_task;
                            solvmtx->btagtab[btagnbr].infotab[BTAG_PROCDST] = min_procdst;
                            solvmtx->btagtab[btagnbr].infotab[BTAG_PRIONUM] = min;

#ifdef DEBUG_BLEND
                            {
                                pastix_int_t taskcnt;
                                taskcnt = 1;
                                for(j=bloknum+1; j<symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum;j++)
                                    if(proc2clust[blprtab[j]]==c)
                                        taskcnt++;
                                ASSERT(taskcnt == solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT],MOD_BLEND);
                            }
                            /*solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT]  = 1;
                             for(j=bloknum+1; j<symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum;j++)
                             if(proc2clust[blprtab[j]]==c)
                             solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT]++;*/


                            ASSERT(proc2clust[ min_procdst ] == c ,MOD_BLEND);
                            if(debug != solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT])
                                fprintf(stdout, " debug %ld taskcnt %ld \n", (long)debug, (long)solvmtx->btagtab[btagnbr].infotab[BTAG_TASKCNT]);
#endif
                            btagnbr++;
                            solvmtx->bcoftab[bcofnbr].sendcnt++;
                            break;
                        }

                solvmtx->bcoftab[bcofnbr].infotab[BCOF_FCOLNUM] = solvmtx->cblktab[cblklocalnum[simuctrl->tasktab[i].cblknum]].fcolnum;
                solvmtx->bcoftab[bcofnbr].infotab[BCOF_LCOLNUM] = solvmtx->cblktab[cblklocalnum[simuctrl->tasktab[i].cblknum]].lcolnum;
                solvmtx->bcoftab[bcofnbr].infotab[BCOF_FROWNUM] = solvmtx->bloktab[bloklocalnum[simuctrl->tasktab[i].bloknum]].frownum;
                solvmtx->bcoftab[bcofnbr].infotab[BCOF_LROWNUM] = solvmtx->bloktab[bloklocalnum[simuctrl->tasktab[i].bloknum]].lrownum;
                solvmtx->bcoftab[bcofnbr].coeftab = NULL;
                bcofnbr++;

                indnbr++;
                break;

            case E2:

                /*
                 * If diag block is local taskmstr is the diag task
                 * else it's the first local E1 task
                 */
            {
                bloknum = simuctrl->tasktab[i].bloknum2;
                if(proc2clust[blprtab[bloknum]] == clustnum)
                {
                    solvmtx->tasktab[tasknum].taskmstr = tasklocalnum[simuctrl->bloktab[bloknum].tasknum];
                }
                else
                {
                    cursor   = tasknum;
                    min      = solvmtx->tasktab[cursor].prionum;
                    min_task = tasknum;
                    while(solvmtx->tasktab[cursor].tasknext != tasknum)
                    {
                        cursor = solvmtx->tasktab[cursor].tasknext;
                        ASSERTDBG(solvmtx->tasktab[cursor].taskid == E2,MOD_BLEND);

                        if(solvmtx->tasktab[cursor].prionum < min)
                        {
                            min      = solvmtx->tasktab[cursor].prionum;
                            min_task = cursor;
                        }
                    }
                    solvmtx->tasktab[tasknum].taskmstr = min_task;
                }
            }

            facebloknum = simuctrl->tasktab[i].facebloknum;

            ASSERT(indnbr == solvmtx->tasktab[tasklocalnum[i]].indnum,MOD_BLEND);

            if(proc2clust[blprtab[facebloknum]]!=clustnum)
            {
                solvmtx->indtab[solvmtx->tasktab[tasklocalnum[i]].indnum]
                    = ftgtlocalnum[CLUST2INDEX(facebloknum, clustnum)];
#ifdef DEBUG_PRIO
                solvmtx->ftgttab[solvmtx->indtab[indnbr]].infotab[FTGT_PRIONUM]
                    = solvmtx->tasktab[tasklocalnum[i]].prionum;
                /*fprintf(stdout, "SOLV Task2D %ld FTGT %ld  taskprio %ld ftgtprio %ld \n", (long)tasklocalnum[i], (long)solvmtx->indtab[indnbr], (long)solvmtx->tasktab[tasklocalnum[i]].prionum, (long)solvmtx->ftgttab[solvmtx->indtab[indnbr]].infotab[FTGT_PRIONUM]);*/
#endif
            }
            else
            {
                solvmtx->indtab[solvmtx->tasktab[tasklocalnum[i]].indnum]
                    = -tasklocalnum[simuctrl->bloktab[facebloknum].tasknum];
            }
            indnbr++;
            break;

            default:
                fprintf(stderr, "Error in solverMatrixgen \n");
                EXIT(MOD_BLEND,INTERNAL_ERR);
            }
            tasknum++;
        }

    }
    ASSERTDBG(indnbr  == solvmtx->indnbr,  MOD_BLEND);
    ASSERTDBG(bcofnbr == solvmtx->bcofnbr, MOD_BLEND);
    ASSERTDBG(btagnbr == solvmtx->btagnbr, MOD_BLEND);

    /*
     *  Initialize data for 2D distribution
     */
    solvmtx->btgsnbr = 0;
    solvmtx->btgrnbr = 0;
    if (ctrl->level2D != 0)
    {
        /* link local btag (blend ???) (attention pb en SMP ???) */
        /* and compute the number of block to send */
        for (i=0; i<solvmtx->btagnbr; i++)
        {
            if (proc2clust[solvmtx->btagtab[i].infotab[BTAG_PROCDST]] != clustnum)
                solvmtx->btgsnbr++;
        }

        /*       for (i=0;i<SOLV_FTGTNBR;i++) */
        /*        ASSERTDBG((FANIN_FCOLNUM(i)!=0) || (FANIN_LCOLNUM(i)!=0),MOD_SOPALIN); */
    }

    /********************/
    /** Fill Solver    **/
    /** cblk and blok  **/
    /********************/
    cblknum = 0;
    bloknum = 0;
    coefnbr = 0;
    coefind = 0;
    for(i=0;i<solvmtx->cblknbr;i++)
    {
        coefind = 0;
        solvmtx->cblktab[i].stride   = 0;
        solvmtx->cblktab[i].procdiag = -1;
        for(j=solvmtx->cblktab[i].bloknum;j<solvmtx->cblktab[i+1].bloknum;j++)
        {
            /* Solvmtx is already expanded in dll */
            delta =  solvmtx->bloktab[j].lrownum - solvmtx->bloktab[j].frownum +1;
            solvmtx->cblktab[i].stride += delta;
            solvmtx->bloktab[j].coefind = coefind;
            coefind += delta;
        }
        coefnbr += (solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum +1) * solvmtx->cblktab[i].stride;
        coefind = coefnbr;
    }
    solvmtx->coefnbr = coefnbr;

    /*****************************************/
    /**  Find coefmax, cpftmax and bpftmax  **/
    /*****************************************/
    /** Find bpftmax **/
    /* bpftmax is the number of coef of the largest block target in reception */
    solvmtx->bpftmax = 0;
    /* the largest block target is the largest local block */


    /** Find cpftmax **/
    /* cpftmax is the number of coef of the largest fan in target in reception */
    solvmtx->cpftmax = 0;

    {
        /***** Find coefmax *****
         * coefmax is the number of coef of the largest temporary block used
         * to compute contribution block on the CPUs.
         * It can be seen as the maximum surface of the
         * C matrix in the GEMM operations.
         */
        pastix_int_t max_m = 0;
        pastix_int_t max_n = 0;

        solvmtx->coefmax = 0;

        if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_YES) {
            for(i=0;i<solvmtx->tasknbr;i++) {
                if(solvmtx->tasktab[i].taskid == COMP_1D) {
                    delta = 0;
                    assert(solvmtx->tasktab[i].cblknum >= 0);
                    for(j=solvmtx->tasktab[i].bloknum;
                        j<solvmtx->cblktab[solvmtx->tasktab[i].cblknum+1].bloknum;
                        j++) {
                        coefind = solvmtx->bloktab[j].lrownum - solvmtx->bloktab[j].frownum+1;
#ifdef PASTIX_ESC
                        while(((j+1) < solvmtx->cblktab[solvmtx->tasktab[i].cblknum+1].bloknum)
                              && (solvmtx->bloktab[j].cblknum == solvmtx->bloktab[j+1].cblknum)) {
                            j++;
                            coefind += solvmtx->bloktab[j].lrownum - solvmtx->bloktab[j].frownum+1;
                        }
#endif
                        if(coefind > delta)
                            delta = coefind;
                    }
                    coefnbr = solvmtx->cblktab[solvmtx->tasktab[i].cblknum].stride * delta;
                    if(coefnbr > solvmtx->coefmax) {
                        solvmtx->coefmax = coefnbr;
                        max_m = solvmtx->cblktab[solvmtx->tasktab[i].cblknum].stride;
                        max_n = delta;
                    }
                }
            }
            /* Maximum size of diagonal blocks */
            for(i=0;i<solvmtx->cblknbr;i++) {
                delta = solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum+1;
                if(delta*delta > solvmtx->coefmax) {
                    solvmtx->coefmax = delta*delta;
                    max_m = delta;
                    max_n = delta;
                }
            }

            fprintf(stderr, "Actual coefmax = %ld (%ld x %ld)\n",
                    (long)solvmtx->coefmax, (long)max_m, (long)max_n );
        }

        /* First compute the maximum size of contribution block */
        solvmtx->coefmax = 0;
        {
            pastix_int_t itercblk;
            pastix_int_t m, n;
            for (itercblk = 0; itercblk < solvmtx->cblknbr; itercblk++) {
                pastix_int_t stride = solvmtx->cblktab[ itercblk ].stride;
                for(i=solvmtx->cblktab[ itercblk ].bloknum+1;
                    i<solvmtx->cblktab[ itercblk + 1].bloknum;i++) {
                    m = stride - solvmtx->bloktab[i].coefind;
                    n = solvmtx->bloktab[i].lrownum - solvmtx->bloktab[i].frownum+1;
#ifdef PASTIX_ESC
                    while(((i+1) < solvmtx->cblktab[itercblk+1].bloknum)
                          && (solvmtx->bloktab[i].cblknum == solvmtx->bloktab[i+1].cblknum)) {
                        i++;
                        n += solvmtx->bloktab[n].lrownum - solvmtx->bloktab[n].frownum+1;
                    }
#endif
                    delta = m * n;
                    if(delta > solvmtx->coefmax) {
                        solvmtx->coefmax = delta;
                        max_m = m;
                        max_n = n;
                    }
                }
                /* kernel_trsm require COLNBR * (stride - COLNBR + 1) in LDLt */
                /* horizontal dimension */
                n = solvmtx->cblktab[itercblk].lcolnum -
                    solvmtx->cblktab[itercblk].fcolnum + 1;
                /* vertical dimension */
                m = stride - n + 1;
                delta = m * n;
                if(delta > solvmtx->coefmax) {
                    solvmtx->coefmax = delta;
                    max_m = m;
                    max_n = n;
                }
            }
        }

        if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_YES) {
            fprintf(stderr, "New suggested coefmax = %ld (%ld x %ld)\n",
                    (long)solvmtx->coefmax, (long)max_m, (long)max_n );
            /* First compute the maximum size of contribution block */
            {
                pastix_int_t max = 0;
                pastix_int_t n;
                for(i=0;i<solvmtx->cblknbr-1;i++) {
                    n = solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum+1;
                    delta = n * 64;
                    if(delta > max) {
                        max = delta;
                        max_m = n;
                        max_n = 64;
                    }
                }
                fprintf(stderr, "Max diagblock coefmax without shur = %ld (%ld x %ld)\n",
                        (long)max, (long)max_m, (long)max_n );

                n = solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum+1;
                delta = n * 64;
                if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_YES)
                    fprintf(stderr, "Max diagblock on shur = %ld (%ld x %ld)\n",
                            (long)delta, (long)n, (long)64 );
            }
        }
    }

    /** Find the cpftmax **/
    /* OIMBE on peut trouver le bon : flemmard */
    solvmtx->cpftmax = 0;
    for(i=0;i<simuctrl->ftgtnbr;i++)
    {
        if((simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR]>0))
            /*&& (proc2clust[simuctrl->ftgttab[i].ftgt.infotab[FTGT_PROCDST]] == clustnum))*/
        {
            coefnbr = (simuctrl->ftgttab[i].ftgt.infotab[FTGT_LCOLNUM] - simuctrl->ftgttab[i].ftgt.infotab[FTGT_FCOLNUM] + 1)*dofptr->noddval * (simuctrl->ftgttab[i].ftgt.infotab[FTGT_LROWNUM] - simuctrl->ftgttab[i].ftgt.infotab[FTGT_FROWNUM] + 1)*dofptr->noddval;
            if(coefnbr > solvmtx->cpftmax)
                solvmtx->cpftmax = coefnbr;
        }
    }

    /** Find the bpftmax **/
    solvmtx->bpftmax = 0;
    for(i=0;i<simuctrl->tasknbr;i++)
    {
        if(simuctrl->tasktab[i].taskid == E1 || simuctrl->tasktab[i].taskid == E2)
            if(proc2clust[simuctrl->blprtab[simuctrl->tasktab[i].bloknum]] == clustnum)
            {
                delta = (symbmtx->cblktab[simuctrl->tasktab[i].cblknum].lcolnum - symbmtx->cblktab[simuctrl->tasktab[i].cblknum].fcolnum+1)
                    * dofptr->noddval;
                coefnbr = delta * (symbmtx->bloktab[simuctrl->tasktab[i].bloknum2].lrownum - symbmtx->bloktab[simuctrl->tasktab[i].bloknum2].frownum+1)
                    *dofptr->noddval;
                if(coefnbr > solvmtx->bpftmax)
                    solvmtx->bpftmax = coefnbr;

            }
    }



    /** Find the nbftmax **/
    solvmtx->nbftmax = 0;
    for(i=0;i<simuctrl->tasknbr;i++)
        if(simuctrl->tasktab[i].ftgtcnt> solvmtx->nbftmax)
            solvmtx->nbftmax = simuctrl->tasktab[i].ftgtcnt;

    /** Find the area max **/

    /** We search the biggest sum of AUB that can be send to the local cluster from a cblk on another cluster **/
    /* @@@@@@@@@@@@ TODO  ****/

    /** We search the biggest cblk that can receive some ftgt **/
    solvmtx->arftmax = 0;
    for(i=0;i<simuctrl->cblknbr;i++)
    {
        pastix_int_t size = 0;
        delta = (symbmtx->cblktab[i].lcolnum -  symbmtx->cblktab[i].fcolnum + 1)*dofptr->noddval;
        if(ctrl->candtab[i].fccandnum != ctrl->candtab[i].lccandnum)
        {
            for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
                size += delta * (symbmtx->bloktab[j].lrownum-symbmtx->bloktab[j].frownum+1) * dofptr->noddval;

            if(size> solvmtx->arftmax)
                solvmtx->arftmax = size;
        }
    }


    if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, "COEFMAX %ld CPFTMAX %ld BPFTMAX %ld NBFTMAX %ld ARFTMAX %ld \n", (long)solvmtx->coefmax, (long)solvmtx->cpftmax,
                (long)solvmtx->bpftmax, (long)solvmtx->nbftmax, (long)solvmtx->arftmax);


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
            if(proc2clust[blprtab[symbmtx->cblktab[i].bloknum]] == clustnum)
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
                    cluster = proc2clust[ blprtab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]] ];
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
                if( proc2clust[ blprtab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]] ] == clustnum)
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
                if( proc2clust[ blprtab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]] ] == clustnum)
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
            if(proc2clust[blprtab[i]] == clustnum)
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
        build_smx(&(solvmtx->updovct), symbmtx, blprtab, ctrl, dofptr);

    }
    /*********************** END TRIANGULAR INFO BUILDING ******************************************/

    memFree_null(cblklocalnum);
    memFree_null(bloklocalnum);
    memFree_null(tasklocalnum);

    if(btaglocalnum != NULL)
        memFree_null(btaglocalnum);

    if(ftgtlocalnum != NULL)
        memFree_null(ftgtlocalnum);

#ifdef DEBUG_BLEND
    for(i=0;i<solvmtx->btagnbr;i++)
        ASSERT(bcofind[i]>=0,MOD_BLEND);
#endif

    return bcofind;
}
