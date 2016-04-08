/**
 *
 * @file sopalin_zgetrf.c
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
#include "isched.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include <dague.h>
#include <dague/data.h>
#include <dague/data_distribution.h>

int dsparse_zgetrf_sp( dague_context_t *dague,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data );
#endif

void
sequential_zgetrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    double              threshold = sopalin_data->diagthreshold;
    pastix_complex64_t *work;
    pastix_int_t  i;
    (void)pastix_data;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    cblk = datacode->cblktab;
    current_cblk = 0;

    /* Compress information coming from A */
    Clock timer;
    clockStart(timer);

    pastix_int_t splitsize     = pastix_data->iparm[IPARM_COMPRESS_SIZE];
    pastix_int_t compress_cblk = splitsize;
    pastix_int_t compress_blok = 10;
    char  *tolerance = getenv("TOLERANCE");
    double tol = atof(tolerance);

    for (i=0; i<datacode->cblknbr; i++, cblk++) {
        pastix_int_t dima, dimb;
        SolverBlok *fblok, *lblok;
        pastix_complex64_t *U, *L;

        /* Horizontal dimension */
        dima = cblk->lcolnum - cblk->fcolnum + 1;
        /* Vertical dimension */
        dimb = cblk->stride  - dima;

        fblok = cblk->fblokptr;   /* this diagonal block */
        lblok = cblk[1].fblokptr; /* the next diagonal block */

        U = cblk->ucoeftab;
        L = cblk->lcoeftab;

        if ( fblok+1 < lblok ){
            pastix_int_t total = dimb;
            SolverBlok *blok   = fblok;
            pastix_int_t bloksize;

            while (total != 0){
                pastix_int_t rankU = -1;
                pastix_int_t rankL = -1;

                blok++;
                bloksize = blok->lrownum - blok->frownum + 1;

                /* TODO: Need to factorize thos two call and reduce the v space to (dima * rkm) */
                if (dima > compress_cblk && bloksize > compress_blok){
                    pastix_complex64_t *uU, *vU;
                    pastix_complex64_t *data;
                    data = U + dima + dimb - total;

                    rankU = pastix_imin( bloksize, dima );
                    uU = malloc( bloksize * rankU * sizeof(pastix_complex64_t));
                    vU = malloc( dima     * dima  * sizeof(pastix_complex64_t));

                    rankU = core_z_compress_LR(tol, bloksize, dima,
                                               data, cblk->stride,
                                               uU, bloksize,
                                               vU, dima);

                    /* Set blok to zero has it will be used has a local memory */
                    /* This is needed when using an extra surface to aggregate some contributions */
                    if (rankU != -1){
                        pastix_int_t i;
                        for (i=0; i<dima; i++){
                            memset(data + cblk->stride * i, 0, bloksize * sizeof(pastix_complex64_t));
                        }
                        blok->coefU_u_LR = uU;
                        blok->coefU_v_LR = vU;
                    }
                }


                if (dima > compress_cblk && bloksize > compress_blok){
                    pastix_complex64_t *uL, *vL;
                    pastix_complex64_t *data;
                    data = L + dima + dimb - total;

                    rankL = pastix_imin( bloksize, dima );
                    uL = malloc( bloksize * rankL * sizeof(pastix_complex64_t));
                    vL = malloc( dima     * dima  * sizeof(pastix_complex64_t));
                    rankL = core_z_compress_LR(tol, bloksize, dima,
                                               data, cblk->stride,
                                               uL, bloksize,
                                               vL, dima);
                    /* Set blok to zero has it will be used has a local memory */
                    /* This is needed when using an extra surface to aggregate some contributions */
                    if (rankL != -1){
                        pastix_int_t i;
                        for (i=0; i<dima; i++){
                            memset(data + cblk->stride * i, 0, bloksize * sizeof(pastix_complex64_t));
                        }
                        blok->coefL_u_LR = uL;
                        blok->coefL_v_LR = vL;
                    }
                }

                blok->rankU = rankU;
                blok->rankL = rankL;
                total -= bloksize;
            }
        }
    }

    clockStop(timer);
    time_comp += clockVal(timer);

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        core_zgetrfsp1d_LR( datacode, cblk, threshold, work );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "getrf_L.txt" );
#endif

    memFree_null( work );
}

void
thread_pzgetrf( int rank, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Wait */
        do {} while( cblk->ctrbcnt );

        /* Compute */
        core_zgetrfsp1d( datacode, cblk, sopalin_data->diagthreshold, work );
    }

#if defined(PASTIX_DEBUG_FACTO) && 0
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
    if (rank == 0) {
        coeftab_zdump( datacode, "getrf_L.txt" );
    }
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
#endif

    memFree_null( work );
}

void
thread_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzgetrf, sopalin_data );
}

#if defined(PASTIX_WITH_PARSEC)
void
parsec_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    dague_context_t *ctx;

    /* Start PaRSEC */
    if (pastix_data->parsec == NULL) {
        int argc = 0;
        pastix_parsec_init( pastix_data, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    /* Run the facto */
    dsparse_zgetrf_sp( ctx, sopalin_data->solvmtx->parsec_desc, sopalin_data );
}
#endif

static void (*zgetrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zgetrf,
    thread_zgetrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zgetrf,
#else
    NULL,
#endif
    NULL
};

void
sopalin_zgetrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zgetrf)(pastix_data_t *, sopalin_data_t *) = zgetrf_table[ sched ];

    if (zgetrf == NULL) {
        zgetrf = thread_zgetrf;
    }
    zgetrf( pastix_data, sopalin_data );
}
