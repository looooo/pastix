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
    pastix_int_t        n;
    pastix_complex64_t *x;
};

static inline void
z_init_smp( isched_thread_t *ctx,
            void            *args )
{
    struct z_argument_init_s *arg = (struct z_argument_init_s*)args;
    pastix_int_t              n   = arg->n;
    pastix_complex64_t       *x   = arg->x;
    pastix_int_t              begin, end, rank, size, nn;

    size  = ctx->global_ctx->world_size;
    rank  = ctx->rank;
    nn    = n / size;

    begin = nn * rank;
    if ( rank == (size - 1) ) {
        end = n;
    }
    else {
        end = nn * (rank + 1);
    }

    assert( begin >= 0 );
    assert( end <= n );
    assert( begin < end );

    core_zplrnt( (end-begin), 1, x + begin, n, n, begin, 0, 7213 );
}

static inline void
z_init( pastix_data_t      *pastix_data,
        pastix_int_t        n,
        pastix_complex64_t *x )
{
    pastix_scheduler_t sched = pastix_data->iparm[IPARM_SCHEDULER];

    if ( sched == PastixSchedSequential ) {
        core_zplrnt( n, 1, x, n, n, 0, 0, 7213 );
    }
    else {
        struct z_argument_init_s args = { n, x };
        isched_parallel_call ( pastix_data->isched, z_init_smp, &args );
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
     * generate matrice of size 'n * nrhs'
     */
    z_init( pastix_data, m * n, A );

    /**
     * generate vectors of size 'nrhs'
     */
    z_init( pastix_data, m, y );
    z_init( pastix_data, n, x );

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
z_bvec_check( pastix_data_t *pastix_data,
              pastix_int_t m )
{
    Clock timer;
    int                 i;
    int                 rc = 0;
    pastix_complex64_t *x, *y;
    pastix_complex64_t  alpha = 3.5;
    struct z_solver     solver;

    z_refine_init( &solver, pastix_data );

    x = malloc( sizeof(pastix_complex64_t) * m );
    y = malloc( sizeof(pastix_complex64_t) * m );

    /**
     * generate vectors of size 'nrhs'
     */
    z_init( pastix_data, m, y );
    z_init( pastix_data, m, x );

    printf("========== copy function time ==========\n");
    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.copy( pastix_data, m, x, y );
    }
    timer = clockGetLocal() - timer;
    printf("    Time for copy (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== axpy function time ==========\n");
    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.axpy( pastix_data, m, alpha, x, y );
    }
    timer = clockGetLocal() - timer;
    printf("    Time for axpy (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== dot function time ==========\n");
    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.dot( pastix_data, m, x, y );
    }
    timer = clockGetLocal() - timer;
    printf("    Time for dot  (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== norm function time ==========\n");
    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.norm( pastix_data, m, x );
    }
    timer = clockGetLocal() - timer;
    printf("    Time for norm (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== scal function time ==========\n");
    timer = clockGetLocal();
    for( i=0; i<50; ++i ) {
        solver.scal( pastix_data, m, alpha, x );
    }
    timer = clockGetLocal() - timer;
    printf("    Time for scal (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

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
    z_init( pastix_data, spm->nexp * nrhs, x );
    z_init( pastix_data, spm->nexp * nrhs, y );

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

    /**
     * Generates the b and x vector such that A * x = b
     * Computes the norms of the initial vectors if checking purpose.
     */
    m     = spm->nexp;
    ldb   = m + 10;
    size  = ldb * nrhs;
    b_in  = malloc( size * sizeof(pastix_complex64_t) );
    b_out = malloc( size * sizeof(pastix_complex64_t) );

    /* Generate B */
    z_init( pastix_data, size, b_in );
    memcpy( b_out, b_in, size * sizeof(pastix_complex64_t) );

    /* Forces the datatypes. */
    spm->flttype = SpmComplex64;
    if ( pastix_data->solvmatr != NULL ) {
        pastix_data->solvmatr->flttype = SpmComplex64;
    }

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

    return rc;
}
