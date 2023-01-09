/**
 *
 * @file z_bvec_tests.c
 *
 * Tests for performance of gemv in parallel and sequential versions.
 * Times is computed on the average of 50 calls to gemv function.
 * Size of the matrix is given by the size of a Laplacian.
 *
 * @copyright 2015-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @date 2022-10-06
 *
 * @precisions normal z -> c d s
 *
 **/
#include <pastix.h>
#include <stdlib.h>
#include <string.h>
#include "refinement/z_refine_functions.h"
#include "common/common.h"
#include "common/flops.h"
#include "bcsc/bcsc.h"
#include "bcsc/bvec.h"
#include "bcsc/bcsc_z.h"
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include <time.h>
#include <pastix/order.h>
#include <cblas.h>
#include <lapacke.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include "z_tests.h"

struct z_argument_init_s
{
    unsigned long long int seed;
    pastix_int_t           gm;
    pastix_int_t           lm;
    pastix_int_t           n;
    pastix_int_t           lda;
    pastix_complex64_t    *A;
    SolverMatrix          *solvmtx;
};

static inline void
z_init_seq( struct z_argument_init_s *arg )
{
    pastix_complex64_t *A       = arg->A;
    SolverMatrix       *solvmtx = arg->solvmtx;

    if ( solvmtx != NULL ) {
        SolverCblk         *cblk    = solvmtx->cblktab;
        pastix_int_t        ii;

        for( ii=0; ii<solvmtx->cblknbr; ii++, cblk++ )
        {
            if ( cblk->cblktype & (CBLK_FANIN | CBLK_RECV) ) {
                continue;
            }

            core_zplrnt( cblk_colnbr( cblk ), arg->n, A + cblk->lcolidx, arg->lda,
                         arg->gm, cblk->fcolnum, 0, arg->seed );
        }
    }
    else {
        core_zplrnt( arg->lm, arg->n, A, arg->lda,
                     arg->gm, 0, 0, arg->seed );
    }
}

static inline void
z_init_smp( isched_thread_t *ctx,
            void            *args )
{
    struct z_argument_init_s *arg     = (struct z_argument_init_s*)args;
    pastix_complex64_t       *A       = arg->A;
    SolverMatrix             *solvmtx = arg->solvmtx;
    pastix_int_t              rank;

    rank = ctx->rank;

    if ( solvmtx != NULL ) {
        pastix_int_t  ii, tid, tasknbr;
        pastix_int_t *tasktab;
        SolverMatrix *solvmtx = arg->solvmtx;
        SolverCblk   *cblk;
        Task         *task;

        tasknbr = solvmtx->ttsknbr[rank];
        tasktab = solvmtx->ttsktab[rank];

        for (ii=0; ii<tasknbr; ii++)
        {
            tid  = tasktab[ii];
            task = solvmtx->tasktab + tid;
            cblk = solvmtx->cblktab + task->cblknum;

            assert( !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) );

            core_zplrnt( cblk_colnbr( cblk ), arg->n, A + cblk->lcolidx, arg->lda,
                         arg->gm, cblk->fcolnum, 0, arg->seed );
        }
    }
    else {
        pastix_int_t begin, end, size, nn;

        size = ctx->global_ctx->world_size;
        nn   = arg->lm / size;

        begin = nn * rank;
        if ( rank == (size - 1) ) {
            end = arg->lm;
        }
        else {
            end = nn * (rank + 1);
        }

        assert( begin >= 0 );
        assert( end <= arg->lm );
        assert( begin < end );

        core_zplrnt( (end-begin), arg->n, A + begin, arg->lda,
                     arg->gm, begin, 0, arg->seed );
    }
}

void
z_init( pastix_data_t      *pastix_data,
        unsigned long long  seed,
        pastix_int_t        n,
        pastix_complex64_t *A,
        pastix_int_t        lda )
{
    pastix_scheduler_t       sched = pastix_data->iparm[IPARM_SCHEDULER];
    pastix_bcsc_t           *bcsc  = pastix_data->bcsc;
    struct z_argument_init_s args  = {
        .seed    = seed,
        .gm      = bcsc ? bcsc->gN : lda,
        .lm      = bcsc ? bcsc->n  : lda,
        .n       = n,
        .lda     = lda,
        .A       = A,
        .solvmtx = pastix_data->solvmatr,
    };

    if ( sched == PastixSchedSequential ) {
        z_init_seq( &args );
    }
    else {
        isched_parallel_call( pastix_data->isched, z_init_smp, &args );
    }
}

int
z_bvec_gemv_check( pastix_data_t *pastix_data,
                   int check, int m, int n )
{
    int                 i;
    int                 rc = 0;
    pastix_complex64_t *A, *x, *y, *y0;
    pastix_complex64_t  alpha = 3.5;
    pastix_complex64_t  beta = -4.8;
    struct z_solver     solver;

    z_refine_init( &solver, pastix_data );

    A = malloc( sizeof(pastix_complex64_t) * (size_t)m * (size_t)n );
    x = malloc( sizeof(pastix_complex64_t) * (size_t)n );
    y = malloc( sizeof(pastix_complex64_t) * (size_t)m );

    /**
     * generate matrice of size 'm * n'
     */
    z_init( pastix_data, 3839, n, A, m );

    /**
     * generate vectors of size 'n' and 'm',
     * y is distributed, x is small and always replicated on all nodes
     */
    z_init( pastix_data, 2308, 1, y, m );
    core_zplrnt( n, 1, x, n, n, 0, 0, 95753 );

    if ( check ) {
        double normA, normX, normY, normR, result;
        double eps = LAPACKE_dlamch_work( 'e' );

        normA = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', m, n, A, m );
        normX = LAPACKE_zlange( LAPACK_COL_MAJOR, '1', n, 1, x, n );
        normY = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', m, 1, y, m );

        y0 = malloc( sizeof(pastix_complex64_t) * m );
        solver.copy( pastix_data, m, y, y0 );

        solver.gemv( pastix_data, m, n, 3.5, A, m, x, -4.8, y );

        cblas_zgemv( CblasColMajor, CblasNoTrans, m, n,
                     CBLAS_SADDR(alpha), A, m, x, 1,
                     CBLAS_SADDR(beta), y0, 1 );

        core_zgeadd( PastixNoTrans, m, 1, -1., y, m, 1., y0, m );

        normR = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', m, 1, y0, m );

        result = normR / ((normA + normX + normY) * n * eps);
        if (  isnan(result) || isinf(result) || (result > 10.0) ) {
            rc = 1;
        }

        free( y0 );
    }
    else {
        double t, flops;
        Clock timer;

        timer = clockGetLocal();
        for ( i = 0; i < 50 ; ++i) {
            solver.gemv( pastix_data, m, n, alpha, A, m, x, beta, y );
        }
        timer = clockGetLocal() - timer;

        t = clockVal(timer) / 50.;
        flops = FLOPS_ZGEMV( m, n ) / t;
        printf("    Time for zgemv ( %ld x %ld ) : %e s ( %8.2g %cFlop/s)\n",
               (long)m, (long)n, t,
               pastix_print_value( flops ), pastix_print_unit( flops ) );
    }
    free( A );
    free( x );
    free( y );

    return rc;
}

int
z_bvec_check( pastix_data_t *pastix_data )
{
    pastix_bcsc_t      *bcsc = pastix_data->bcsc;
    pastix_int_t        lm, gm, n = 1;
    int                 lrc, grc, rc = 0;
    pastix_complex64_t *x, *y;
    pastix_complex64_t *x0, *y0;
    pastix_complex64_t  alpha = 3.5;
    struct z_solver     solver;
    unsigned long long  seedX = 3927;
    unsigned long long  seedY = 9846;
    double              eps, result;

    eps = LAPACKE_dlamch_work('e');

    z_refine_init( &solver, pastix_data );

    gm = bcsc->gN;
    lm = bcsc->n;

    x = malloc( sizeof(pastix_complex64_t) * lm * n );
    y = malloc( sizeof(pastix_complex64_t) * lm * n );
    memset( x, 0xff, sizeof(pastix_complex64_t) * lm * n );
    memset( y, 0xff, sizeof(pastix_complex64_t) * lm * n );

    /*
     * generate vectors
     */
    z_init( pastix_data, seedX, n, x, lm );

    /*
     * Check the copy
     */
    {
        if ( pastix_data->procnum == 0 ) {
            printf(" == Check copy function : ");
        }
        solver.copy( pastix_data, lm, x, y );

        lrc = z_bvec_compare( pastix_data, lm, n, x, lm, y, lm );
        MPI_Allreduce( &lrc, &grc, 1, MPI_INT, MPI_SUM, pastix_data->inter_node_comm );

        if ( pastix_data->procnum == 0 ) {
            printf( grc == 0 ? "SUCCESS\n" : "FAILURE\n" );
        }
        rc += grc;
    }

    /* Generate a backup of y to compute through cblas in order to compare with */
    y0 = malloc( sizeof(pastix_complex64_t) * lm * n );
    memset( y0, 0xff, sizeof(pastix_complex64_t) * lm * n );
    z_init( pastix_data, seedY, n, y,  lm );
    z_init( pastix_data, seedY, n, y0, lm );

    /*
     * Check the axpy
     */
    {
        if ( pastix_data->procnum == 0 ) {
            printf(" == Check axpy function : ");
        }
        solver.axpy( pastix_data, lm, alpha, x, y );

        cblas_zaxpy( lm * n, CBLAS_SADDR(alpha), x, 1, y0, 1 );

        lrc = z_bvec_compare( pastix_data, lm, n, y, lm, y0, lm );
        MPI_Allreduce( &lrc, &grc, 1, MPI_INT, MPI_SUM, pastix_data->inter_node_comm );

        if ( pastix_data->procnum == 0 ) {
            printf( grc == 0 ? "SUCCESS\n" : "FAILURE\n" );
        }
        rc += grc;
    }

    /* Restore y and y0 */
    z_init( pastix_data, seedY, n, y,  lm );
    z_init( pastix_data, seedY, n, y0, lm );

    /*
     * Check the scal
     */
    {
        if ( pastix_data->procnum == 0 ) {
            printf(" == Check scal function : ");
        }

        solver.scal( pastix_data, lm, alpha, y );

        cblas_zscal( lm * n, CBLAS_SADDR(alpha), y0, 1 );

        lrc = z_bvec_compare( pastix_data, lm, n, y, lm, y0, lm );
        MPI_Allreduce( &lrc, &grc, 1, MPI_INT, MPI_SUM, pastix_data->inter_node_comm );

        if ( pastix_data->procnum == 0 ) {
            printf( grc == 0 ? "SUCCESS\n" : "FAILURE\n" );
        }
        rc += grc;
    }

    free( y0 );

    /* Restore x0, y, and y0 */
    z_init( pastix_data, seedY, n, y, lm );
    if ( pastix_data->procnum == 0 ) {
        x0 = malloc( sizeof(pastix_complex64_t) * gm );
        y0 = malloc( sizeof(pastix_complex64_t) * gm );
        core_zplrnt( gm, n, x0, gm, gm, 0, 0, seedX );
        core_zplrnt( gm, n, y0, gm, gm, 0, 0, seedY );
    }

    /*
     * Check the dot
     */
    {
        pastix_complex64_t ldotc, gdotc;

        if ( pastix_data->procnum == 0 ) {
            printf(" == Check dot function : ");
        }

        /* Restore y */
        ldotc = solver.dot( pastix_data, lm, x, y );

        if ( pastix_data->procnum == 0 ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
            cblas_zdotc_sub( gm, x0, 1, y0, 1, &gdotc );
#else
            gdotc = cblas_zdotc( gm, x0, 1, y0, 1 );
#endif
        }

        MPI_Bcast( &gdotc, 1, MPI_DOUBLE_COMPLEX, 0, pastix_data->inter_node_comm );
        result = cabs(ldotc - gdotc) / (gm * eps);
        if (  isinf(cabs(ldotc)) || isinf(cabs(ldotc)) ||
              isnan(result) || isinf(result) || (result > 10.0) )
        {
            lrc = 1;
        }
        else {
            lrc = 0;
        }

        MPI_Allreduce( &lrc, &grc, 1, MPI_INT, MPI_SUM, pastix_data->inter_node_comm );
        if ( pastix_data->procnum == 0 ) {
            if ( grc == 0 ) {
                printf( "SUCCESS\n" );
            }
            else {
#if defined(PRECISION_z) || defined(PRECISION_c)
                printf( "FAILURE\n local dot = %e + i * %e, global dot = %e + i * %e ( %e )\n",
                        creal(ldotc), cimag(ldotc), creal(gdotc), cimag(gdotc), result );
#else
                printf( "FAILURE\n local dot = %e, global dot = %e ( %e )\n",
                        ldotc, gdotc, result );
#endif
            }
        }
        rc += grc;
    }

    /*
     * Check the norm
     */
    {
        double lnorm, gnorm;

        if ( pastix_data->procnum == 0 ) {
            printf(" == Check norm function : ");
        }

        lnorm = solver.norm( pastix_data, lm, x );

        if ( pastix_data->procnum == 0 ) {
            gnorm = cblas_dznrm2( gm, x0, 1 );
        }

        MPI_Bcast( &gnorm, 1, MPI_DOUBLE, 0, pastix_data->inter_node_comm );
        result = fabs(lnorm - gnorm) / (gm * eps);
        if (  isinf(lnorm) || isinf(lnorm) ||
              isnan(result) || isinf(result) || (result > 10.0) ) {
            lrc = 1;
        }
        else {
            lrc = 0;
        }

        MPI_Allreduce( &lrc, &grc, 1, MPI_INT, MPI_SUM, pastix_data->inter_node_comm );
        if ( pastix_data->procnum == 0 ) {
            if ( grc == 0 ) {
                printf( "SUCCESS\n" );
            }
            else {
                printf( "FAILURE\n local norm = %e, global norm = %e ( %e )\n",
                        lnorm, gnorm, result );
            }
        }
        rc += grc;
    }

    if ( pastix_data->procnum == 0 ) {
        free( x0 );
        free( y0 );
    }

    free( x );
    free( y );

    return rc;
}

int
z_bvec_time( pastix_data_t *pastix_data )
{
    pastix_bcsc_t      *bcsc = pastix_data->bcsc;
    Clock               timer;
    pastix_int_t        lm, gm, n = 1;
    int                 i;
    int                 rc = 0;
    pastix_complex64_t *x, *y;
    pastix_complex64_t  alpha = 3.5;
    struct z_solver     solver;

    z_refine_init( &solver, pastix_data );

    gm = bcsc->gN;
    lm = bcsc->n;

    x = malloc( sizeof(pastix_complex64_t) * lm * n );
    y = malloc( sizeof(pastix_complex64_t) * lm * n );

    /**
     * generate vectors
     */
    z_init( pastix_data, 3087, n, y, lm );
    z_init( pastix_data, 8998, n, x, lm );

    if ( pastix_data->procnum == 0 ) {
        printf("========== copy function time ==========\n");
    }

    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.copy( pastix_data, lm, x, y );
    }
    timer = clockGetLocal() - timer;

    if ( pastix_data->procnum == 0 ) {
        printf("    Time for copy (N= %ld)         : %e \n", (long)gm, clockVal(timer)/50.);
        printf("========== axpy function time ==========\n");
    }

    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.axpy( pastix_data, lm, alpha, x, y );
    }
    timer = clockGetLocal() - timer;

    if ( pastix_data->procnum == 0 ) {
        printf("    Time for axpy (N= %ld)         : %e \n", (long)gm, clockVal(timer)/50.);
        printf("========== dot function time ==========\n");
    }

    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.dot( pastix_data, lm, x, y );
    }
    timer = clockGetLocal() - timer;

    if ( pastix_data->procnum == 0 ) {
        printf("    Time for dot  (N= %ld)         : %e \n", (long)gm, clockVal(timer)/50.);
        printf("========== norm function time ==========\n");
    }

    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.norm( pastix_data, lm, x );
    }
    timer = clockGetLocal() - timer;

    if ( pastix_data->procnum == 0 ) {
        printf("    Time for norm (N= %ld)         : %e \n", (long)gm, clockVal(timer)/50.);
        printf("========== scal function time ==========\n");
    }

    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.scal( pastix_data, lm, alpha, x );
    }
    timer = clockGetLocal() - timer;

    if ( pastix_data->procnum == 0 ) {
        printf("    Time for scal (N= %ld)         : %e \n", (long)gm, clockVal(timer)/50.);
    }

    free( x );
    free( y );

    return rc;
}

int
z_bcsc_spmv_time( pastix_data_t    *pastix_data,
                  const spmatrix_t *spm,
                  pastix_int_t      nrhs )
{
    int                 i;
    int                 rc = 0;
    pastix_complex64_t *x, *y;
    pastix_complex64_t  alpha = 3.5;
    pastix_complex64_t  beta = -4.8;
    struct z_solver     solver;
    double t, flops;
    Clock timer;

    z_refine_init( &solver, pastix_data );

    x = malloc( sizeof(pastix_complex64_t) * spm->nexp * nrhs );
    y = malloc( sizeof(pastix_complex64_t) * spm->nexp * nrhs );

    /**
     * generate matrices
     */
    z_init( pastix_data, 9794, nrhs, x, spm->nexp );
    z_init( pastix_data, 7839, nrhs, y, spm->nexp );

    timer = clockGetLocal();
    for ( i = 0; i < 50 ; ++i) {
        solver.spmv( pastix_data, PastixNoTrans, alpha, x, beta, y );
    }
    timer = clockGetLocal() - timer;

    t = clockVal(timer) / 50.;
    flops = FLOPS_ZGEMM( spm->nnzexp, 1, 1 ) / t;
    printf("    Time for zspmv ( n=%ld; nnz=%ld ) : %e s ( %8.2g %cFlop/s)\n",
           (long)spm->nexp, (long)spm->nnzexp, t,
           pastix_print_value( flops ), pastix_print_unit( flops ) );

    free( x );
    free( y );

    return rc;
}

int
z_bvec_compare( pastix_data_t            *pastix_data,
                pastix_int_t              m,
                pastix_int_t              n,
                const pastix_complex64_t *A,
                pastix_int_t              lda,
                const pastix_complex64_t *B,
                pastix_int_t              ldb )
{
    pastix_int_t i, j;
    int rc_loc  = 0;
    int rc_glob = 0;

    for ( j = 0; j < n; j++ ) {
        for ( i = 0; i < m; i++, A++, B++ ) {
            if ( *A != *B ) {
                rc_loc++;
                goto reduce;
            }
        }

        A += lda - m;
        B += ldb - m;
    }

  reduce:
    MPI_Allreduce( &rc_loc, &rc_glob, 1, MPI_INT, MPI_SUM, pastix_data->inter_node_comm );

    (void)pastix_data;
    return rc_glob;
}

int
z_bvec_applyorder_check ( pastix_data_t *pastix_data,
                          spmatrix_t    *spm,
                          pastix_int_t   nrhs )
{
    pastix_complex64_t *b_in, *b_out;
    pastix_rhs_t        Bp;
    size_t              size;
    pastix_int_t        m, ldb;
    int                 rc = 0;
    pastix_int_t        sum_m, my_sum;
    pastix_int_t        clustnum = spm->clustnum;
    pastix_int_t        clustnbr = spm->clustnbr;

    /**
     * Generates the b and x vector such that A * x = b
     * Computes the norms of the initial vectors if checking purpose.
     */
    m     = spm->nexp;
    ldb   = m + 10;
    size  = ldb * nrhs;
    b_in  = malloc( size * sizeof(pastix_complex64_t) );
    b_out = malloc( size * sizeof(pastix_complex64_t) );
    memset( b_in,  0xff, sizeof(pastix_complex64_t) * size );
    memset( b_out, 0xff, sizeof(pastix_complex64_t) * size );

    sum_m  = 0;
    my_sum = 0;

    if ( clustnum > 0 ) {
        MPI_Recv( &my_sum, 1, PASTIX_MPI_INT, clustnum - 1, PastixTagAmount, spm->comm, MPI_STATUSES_IGNORE );
    }
    if ( clustnum < ( clustnbr - 1 ) ) {
        sum_m = my_sum + m;
        MPI_Send( &sum_m, 1, PASTIX_MPI_INT, clustnum + 1, PastixTagAmount, spm->comm );
    }

    /* Generates a replicated B */
    if ( spm->loc2glob == NULL ) {
        core_zplrnt( m, nrhs, b_in, ldb, m, 0, 0, 4978 );
    }
    else {
        core_zplrnt( m, nrhs, b_in, ldb, spm->gNexp, my_sum, 0, 4978 );
    }
    memcpy( b_out, b_in, size * sizeof(pastix_complex64_t) );

    /* Forces the datatypes. */
    spm->flttype = SpmComplex64;
    assert( pastix_data->solvmatr != NULL );
    pastix_data->solvmatr->flttype = SpmComplex64;
    assert( pastix_data->bcsc != NULL );
    pastix_data->bcsc->flttype = SpmComplex64;

    /* Make sure internal peritab and spm are reset */
    pastix_data->csc = spm;
    if ( pastix_data->ordemesh->peritab_exp ) {
        memFree_null( pastix_data->ordemesh->peritab_exp );
    }

    /* Compute P * b_in */
    pastixRhsInit( &Bp );
    pastix_subtask_applyorder( pastix_data, PastixDirForward, m, nrhs, b_out, ldb, Bp );

    /* Compute b_out = P^t * (P * b_in) */
    pastix_subtask_applyorder( pastix_data, PastixDirBackward, m, nrhs, b_out, ldb, Bp );
    pastixRhsFinalize( Bp );

    /* Checks the result. */
    rc = z_bvec_compare( pastix_data, m, nrhs, b_in, ldb, b_out, ldb );

    free( b_in );
    free( b_out );

    (void)sum_m;
    return rc;
}
