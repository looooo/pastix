/**
 *
 * @file sequential_zpotrf.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "sopalin_data.h"

int
core_zpotrfsp1d( SolverMatrix *solvmtx,
                 SolverCblk   *cblk,
                 double        criteria );

void
pastix_static_zpotrf( sopalin_data_t *sopalin_data )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk *cblk;
    pastix_int_t  i, ii;
    pastix_int_t tasknbr, *tasktab;
    Task *t;
    FILE *stream = fopen("facto.txt", "w");

    tasknbr = datacode->ttsknbr[0];
    tasktab = datacode->ttsktab[0];

    for (ii=0; ii<tasknbr; ii++){
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Compute task */
        switch( t->taskid )
        {
        case COMP_1D:
                    /* /\* Wait for contributions *\/ */
                    /* API_CALL(z_wait_contrib_comp_1d)(sopalin_data, me, i); */

                    /* z_ooc_wait_task(sopalin_data,i, me); */
                    /* z_ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(i),me); */

                    /* /\* trace_begin_task1(thread_data->tracefile, *\/ */
                    /* /\*                   SOPALIN_CLOCK_TRACE, *\/ */
                    /* /\*                   SOLV_PROCNUM, me, *\/ */
                    /* /\*                   STATE_COMP1D, *\/ */
                    /* /\*                   TASK_PROC( i ), *\/ */
                    /* /\*                   SOLV_TASKTAB[i], *\/ */
                    /* /\*                   stolen ); *\/ */

            /* Compute */
            core_zpotrfsp1d( datacode, cblk, 0.002397);
            fprintf(stream, "i=%ld, task=%ld, cblk=%ld\n",
                    (long)ii, (long)i, (long)t->cblknum);
            coeftab_zdumpcblk( datacode, cblk, stream );
            break;

        default:
            errorPrint("Taskid unknown for task %ld\n", (long)i);
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
    }
    fclose(stream);
}
