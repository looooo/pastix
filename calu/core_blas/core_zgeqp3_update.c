/**
 *
 * @file core_zgeqp3_update.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

#define A(m,n) BLKADDR( A, PLASMA_Complex64_t, m, n )

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 * CORE_zgeqp3_update updates row k of one tile of A
 * and subtracts that row from the column norms.
 *
 *******************************************************************************
 *
 * @param[in] Ajj
 *         Diagonal tile (jj,jj) of A.
 *
 * @param[in] lda1
 *         Leading dimension of Ajj.
 *
 * @param[in,out] Ajk
 *         Tile (jj,kk) of A, kk >= jj.
 *         On exit, updates row k (i.e., as if Q was applied to trailing matrix).
 *
 * @param[in] lda2
 *         Leading dimension of Ajk.
 *
 * @param[in] Fk
 *         Tile kk of F.
 *
 * @param[in] ldf
 *         Leading dimension of Fk.
 *
 * @param[in] k
 *         Row within tile to update, 0 <= k < nb.
 *
 * @param[in] j
 *         Column to start updating; for diagonal tile, j=k+1, else j=0.
 *
 * @param[in] nb
 *         Number of columns in kk-th block-column of A.
 *
 * @param[in,out] norms1
 *         kk-th block of partial column norms vector, dimension nb.
 *         On exit, norms1 := norms1 - Ajk[k,:].
 *
 * @param[in] norms2
 *         kk-th block of original column norms vector, dimension nb.
 *
 * @param[out] INFO
 *         See returned value.
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval >0 if i, a numerical instability was detected at column i.
 *                    This case is not currently handled, but in the future
 *                    it will be handled internally.
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_update = PCORE_zgeqp3_update
#define CORE_zgeqp3_update PCORE_zgeqp3_update
#endif
void CORE_zgeqp3_update( const PLASMA_Complex64_t *Ajj, int lda1,
                         PLASMA_Complex64_t       *Ajk, int lda2,
                         const PLASMA_Complex64_t *Fk,  int ldf,
                         int k, int j, int nb,
                         double *norms1, const double *norms2,
                         int *info )
{
    int j2;
    double temp, temp2;
    double tol3z = sqrt( LAPACKE_dlamch_work('e'));
    const PLASMA_Complex64_t zone  =  1.0;
    const PLASMA_Complex64_t mzone = -1.0;
    
    *info = 0;
    
    /* update row k of A -- this is vector*matrix */
    /* Ajk[k,j:nb] -= Ajj[k,0:k+1] * Fk[j:nb,0:k+1].T */
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasConjTrans, 1, nb-j, k+1,
                 CBLAS_SADDR(mzone), &Ajj[k],          lda1,
                                     &Fk [j],          ldf,
                 CBLAS_SADDR(zone),  &Ajk[k + j*lda2], lda2 );
    
    for( j2 = j; j2 < nb; ++j2 ) {
        if ( norms1[j2] != 0. ) {
            /* NOTE: The following lines follow from the analysis in Lapack Working Note 176. */
            temp = cabs( Ajk[k + j2*lda2] ) / norms1[j2];
            temp = max( 0., (1. + temp)*(1. - temp) );
            temp2 = norms1[j2] / norms2[j2];
            temp2 = temp * temp2*temp2;
            if( temp2 <= tol3z ) {
                /* todo handle numerical instability. current just report that it occurred */
                *info = j2;
            }
            else {
                norms1[j2] = norms1[j2]*sqrt( temp );
            }
        }
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgeqp3_update( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_Complex64_t *Ajj, int lda1,
                               PLASMA_Complex64_t *Ajk, int lda2,
                               PLASMA_Complex64_t *Fk,  int ldf,
                               int k, int j, int nb,
                               double *norms1, const double *norms2,
                               PLASMA_sequence *sequence, PLASMA_request *request )
{
    DAG_SET_PROPERTIES("update", "magenta");
    QUARK_Insert_Task(quark, CORE_zgeqp3_update_quark, task_flags,
        sizeof(PLASMA_Complex64_t)*nb*nb,  Ajj,             INPUT,
        sizeof(int),                       &lda1,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,  Ajk,                     INOUT,
        sizeof(int),                       &lda2,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,  Fk,              INPUT,
        sizeof(int),                       &ldf,    VALUE,
        sizeof(int),                       &k,      VALUE,
        sizeof(int),                       &j,      VALUE,
        sizeof(int),                       &nb,     VALUE,
        sizeof(double)*nb,                 norms1,                  INOUT,
        sizeof(double)*nb,                 norms2,                  NODEP,  /* INOUT, but implied by norms1 */
        sizeof(PLASMA_sequence*),          &sequence, VALUE,
        sizeof(PLASMA_request*),           &request,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_update_quark = PCORE_zgeqp3_update_quark
#define CORE_zgeqp3_update_quark PCORE_zgeqp3_update_quark
#endif
void CORE_zgeqp3_update_quark( Quark *quark )
{
    const PLASMA_Complex64_t *Ajj, *Fk;
    PLASMA_Complex64_t *Ajk;
    int lda1, lda2, ldf, k, j, nb;
    double *norms1;
    const double *norms2;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    
    int info;
    quark_unpack_args_13( quark, Ajj, lda1, Ajk, lda2, Fk, ldf, k, j, nb, norms1, norms2, sequence, request );
    CORE_zgeqp3_update(          Ajj, lda1, Ajk, lda2, Fk, ldf, k, j, nb, norms1, norms2, &info );
    if (info != PLASMA_SUCCESS) {
        plasma_sequence_flush(quark, sequence, request, info);
    }
}
