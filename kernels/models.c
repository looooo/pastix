/**
 *
 * @file models.c
 *
 *  PaStiX modelling routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2017-02-09
 *
 **/
#include "common.h"
#include "solver.h"
#include "models.h"

pastix_model_entry_t *model_entries = NULL;                /*< Array to all entries                 */
volatile int32_t      model_entries_nbr = -1;              /*< Index of the last entry in the array */
int32_t               model_size = 0;

void
modelInit( const SolverMatrix *solvmtx )
{
    pastix_int_t cblknbr   = solvmtx->cblknbr;
    pastix_int_t cblkmin2d = solvmtx->cblkmin2d;
    pastix_int_t total_number_of_tasks = 0;
    pastix_int_t nbfact, nbtrsm, nbgemm;
    pastix_int_t cblknum;
    SolverCblk *cblk;

    /* Factorization kernels */
    nbfact = cblknbr;

    /* TRSM kernels */
    nbtrsm = cblkmin2d + (cblknbr - cblkmin2d) * solvmtx->cblkmaxblk;
    if ( solvmtx->factotype == PastixFactLU ) {
        nbtrsm *= 2;
    }

    /* GEMM kernels */
    nbgemm = solvmtx->bloknbr - cblknbr;
    if ( solvmtx->factotype == PastixFactLU ) {
        nbgemm *= 2;
    }

    cblk = solvmtx->cblktab+cblkmin2d;
    for(cblknum = cblkmin2d; cblknum < cblknbr; cblknum++, cblk++ ) {
        pastix_int_t nbodb = (cblk[1].fblokptr - cblk[0].fblokptr) - 1;

        if ( solvmtx->factotype == PastixFactLU ) {
            nbgemm += nbodb * nbodb;
        }
        else {
            nbgemm += (nbodb * (nbodb-1)) / 2;
        }
    }

    total_number_of_tasks = nbfact + nbtrsm + nbgemm;
    model_entries = malloc( total_number_of_tasks * sizeof(pastix_model_entry_t) );
    model_size = total_number_of_tasks;

    return;
}

void
modelDumpAndExit( char **directory )
{
    pastix_model_entry_t *entry = model_entries;
    pastix_int_t i;
    FILE *f;

    f = fopen( "model.csv", "a" );
    if ( f == NULL ) {
        return;
    }

    for(i=0; i <= model_entries_nbr; i++, entry++ ) {
        fprintf( f, "%d;%d;%d;%d;%e\n",
                 entry->ktype, entry->m, entry->n, entry->k, entry->time );
    }

    fclose( f );
    free( model_entries );
    model_entries = NULL;
    model_entries_nbr = -1;
    model_size = 0;

    (void)directory;
}
