/**
 *
 * @file pzgeqp3.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define PRECISION_z

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)

static const PLASMA_Complex64_t zone  = (PLASMA_Complex64_t)  1.;
static const PLASMA_Complex64_t mzone = (PLASMA_Complex64_t) -1.;
static const PLASMA_Complex64_t zzero = (PLASMA_Complex64_t)  0.;


/***************************************************************************//**
 *  Parallel tile QR factorization with column pivoting - dynamic scheduling
 **/
void plasma_pzgeqp3_quark( PLASMA_desc A, int *jpvt, PLASMA_Complex64_t *tau,
                           PLASMA_Complex64_t *work, double *rwork,
                           PLASMA_sequence *sequence, PLASMA_request *request )
{
    PLASMA_Complex64_t *Aij, *Ajj, *Aik, *Ajk, *Fk;
    int jj, j, nj, jk, ldf, i2, ii, k, j2, kk, lda, mbi, nbj, nbk;
    int n;
    int tempnn;

    /* TODO: must be static to exist after return. Should be put into workspace. */
    static PLASMA_Complex64_t beta;

    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    PLASMA_Complex64_t *F = work;
    ldf = A.n;
    PLASMA_Complex64_t *aux = &work[A.nb*ldf];
    double *norms1 = rwork;
    double *norms2 = rwork + A.n;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for( j = 0; j < A.n; ++j ) {
        jpvt[j] = j;
    }

    /* Compute the initial norm of each column of A */
    for( n = 0; n < A.nt; n++ ) {
        tempnn = n == A.nt-1 ? A.n - n * A.nb : A.nb;
        QUARK_CORE_zgeqp3_norms(
            plasma->quark, &task_flags,
            plasma_desc_submatrix(A, 0, n*A.nb, A.m, tempnn),
            norms1 + n*A.nb,
            norms2 + n*A.nb );
    }

    for( jj = 0; jj < A.nt; ++jj ) {
        Ajj = A(jj,jj);

        j  = jj*A.nb;
        nj = A.n - j;
        nbj = min( A.nb, nj );

        /* F = 0, so gemv can accumulate into F */
        for( kk = jj; kk < A.nt; ++kk ) {
            nbk = min( A.nb, A.n - kk*A.nb );
            Fk = &F[(kk-jj)*A.nb];
            QUARK_CORE_zlaset( plasma->quark, &task_flags,
                PlasmaUpperLower, nbk, nbj,
                zzero, zzero, Fk, ldf );
        }

        for( k = 0; k < nbj; ++k ) {
            jk = j + k;

            QUARK_CORE_zgeqp3_pivot( plasma->quark, &task_flags,
                A, F, ldf, jj, k, jpvt, norms1, norms2 );

            /* update current column */
            if ( k > 0 ) {
                /* A[jk:m,jk] -= A[jk:m,j:jk] * F[k,0:k]^H */
                i2 = k;
                for( ii = jj; ii < A.mt; ++ii ) {
                    mbi = min( A.mb, A.m - ii*A.mb );
                    lda = BLKLDD( A, ii );
                    Aij = A(ii,jj);
                    /* really gemv, but F needs to be conjugated, so use gemm */
                    QUARK_CORE_zgemm_tile( plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaConjTrans, mbi-i2, 1, k, A.nb,
                        &mzone, &Aij[i2],         lda,
                                &F[k],            ldf,
                        &zone,  &Aij[i2 + k*lda], lda,
                        Aij, F, Aij );  /* tile dependencies */
                    i2 = 0;
                }
            }

            /* Householder */
            QUARK_CORE_zgeqp3_larfg( plasma->quark, &task_flags,
                A, jj, jj, k, k, &tau[jk], &beta );

            /* compute F */
            if ( k < nj-1 ) {
                /* F[k+1:nj,k] = tau[jk] * A[jk:m,jk+1:n]^T * A[jk:m,jk] */
                /* assumes F = 0, above */
                i2 = k;
                for( ii = jj; ii < A.mt; ++ii ) {
                    mbi = min( A.mb, A.m - ii*A.mb );
                    lda = BLKLDD( A, ii );
                    Aij = A(ii,jj);
                    j2  = k+1;
                    for( kk = jj; kk < A.nt; ++kk ) {
                        nbk = min( A.nb, A.n - kk*A.nb );
                        Aik = A(ii,kk);
                        Fk  = &F[(kk-jj)*A.nb];
                        QUARK_CORE_zgemv_tile( plasma->quark, &task_flags,
                            PlasmaConjTrans, mbi-i2, nbk-j2,
                            &tau[jk], &Aik[i2 + j2*lda], lda,
                                      &Aij[i2 +  k*lda], 1,
                            &zone,    & Fk[j2 +  k*ldf], 1,
                            Aik, Aij, Fk );  /* tile dependencies */
                        j2 = 0;
                    }
                    i2 = 0;
                }
            }

            /* incremental update of F */
            if ( k > 0 ) {
                /* aux = - A[jk:m,j:jk]^T * A[jk:m,jk] */
                i2 = k;
                for( ii = jj; ii < A.mt; ++ii ) {
                    mbi = min( A.mb, A.m - ii*A.mb );
                    lda = BLKLDD( A, ii );
                    Aij = A(ii,jj);
                    if ( ii == jj ) {
                        /* first gemv, beta=0; later beta=1 */
                        QUARK_CORE_zgemv_tile( plasma->quark, &task_flags,
                            PlasmaConjTrans, mbi-i2, k,
                            &mzone, &Aij[i2], lda,
                                    &Aij[i2 + k*lda], 1,
                            &zzero, aux, 1,
                            Aij, Aij, aux );  /* tile dependencies */
                    }
                    else {
                        QUARK_CORE_zgemv_tile( plasma->quark, &task_flags,
                            PlasmaConjTrans, mbi-i2, k,
                            &mzone, &Aij[i2], lda,
                                    &Aij[i2 + k*lda], 1,
                            &zone,  aux, 1,
                            Aij, Aij, aux );  /* tile dependencies */
                    }
                    i2 = 0;
                }
                /* F[0:nj,k] += tau[jk] * F[0:nj,0:k] * aux */
                /* (compared to lapack, minus moved above into aux) */
                for( kk = jj; kk < A.nt; ++kk ) {
                    nbk = min( A.nb, A.n - kk*A.nb );
                    Fk = &F[(kk-jj)*A.nb];
                    QUARK_CORE_zgemv_tile( plasma->quark, &task_flags,
                        PlasmaNoTrans, nbk, k,
                        &tau[jk], Fk, ldf,
                                  aux, 1,
                        &zone,    &Fk[k*ldf], 1,
                        Fk, aux, Fk );  /* tile dependencies */
                }
            }

            /* update pivot row and norms */
            /* A[jk,jk+1:n] -= A[jk,j:jk+1] * F[k+1:nj,0:k+1]^T */
            /* norms1[jk+1:n] =  sqrt( norms1[jk+1:n]**2 - A[jk,jk+1:n]**2 ) */
            j2 = k+1;
            lda = BLKLDD( A, jj );
            for( kk = jj; kk < A.nt; ++kk ) {
                nbk = min( A.nb, A.n - kk*A.nb );
                Ajk = A(jj,kk);
                Fk  = &F[(kk-jj)*A.nb];
                QUARK_CORE_zgeqp3_update( plasma->quark, &task_flags,
                    Ajj, lda, Ajk, lda, Fk, ldf, k, j2, nbk,
                    &norms1[kk*A.nb], &norms2[kk*A.nb],
                    sequence, request );
                j2 = 0;
            }

            /* save beta (from zlarfg) to diagonal */
            QUARK_CORE_zsetvar( plasma->quark, &task_flags,
                &beta, &Ajj[k + k*lda], Ajj );
        }

        /* trailing matrix update */
        for( ii = jj+1; ii < A.mt; ++ii ) {
            mbi = min( A.mb, A.m - ii*A.mb );
            lda = BLKLDD( A, ii );
            for( kk = jj+1; kk < A.nt; ++kk ) {
                nbk = min( A.nb, A.n - kk*A.nb );

                Aik = A(ii,kk);          /* Aik  mbi x nbk */
                Aij = A(ii,jj);          /* Aik  mbi x nbj */
                Fk  = &F[(kk-jj)*A.nb];  /* Fk^T nbj x nbk */
                QUARK_CORE_zgemm( plasma->quark, &task_flags,
                    PlasmaNoTrans, PlasmaConjTrans, mbi, nbk, nbj, A.nb,
                    mzone, Aij, lda,
                           Fk,  ldf,
                    zone,  Aik, lda );
            }
        }
    }
}
