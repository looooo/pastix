/**
 *
 * @file core_zgbrce.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Azzam Haidar
 * @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

#define A(_m, _n)  (PLASMA_Complex64_t *)plasma_geteltaddr(A, ((_m)-1), ((_n)-1), eltsize)
#define V(_m)      &(V[(_m)-1])
#define TAU(_m)    &(TAU[(_m)-1])

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zgbrce is a kernel that will operate on a region (triangle) of data bounded
 *  by st and ed. This kernel apply a right update, create a new nnz, then it
 *  eliminate it, and move to the next right update, create a new nnz, eliminate it
 *  and so on until finishing. When this is done, it take advantage that data are on
 *  cache and will apply the left on the remaining part of this region
 *  that has not been updated by the left yet.
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower:
 *         @arg PlasmaUpper:
 *
 * @param[in] N
 *          The order of the matrix A.
 *
 * @param[in, out] A
 *          A pointer to the descriptor of the matrix A.
 *
 * @param[out] V
 *          PLASMA_Complex64_t array, dimension (N).
 *          The scalar elementary reflectors are written in this
 *          array. So it is used as a workspace for V at each step
 *          of the bulge chasing algorithm.
 *
 * @param[out] TAU
 *          PLASMA_Complex64_t array, dimension (N).
 *          The scalar factors of the elementary reflectors are written
 *          in thisarray. So it is used as a workspace for TAU at each step
 *          of the bulge chasing algorithm.
 *
 * @param[in] st
 *          A pointer to the start index where this kernel will operate.
 *
 * @param[in] ed
 *          A pointer to the end index where this kernel will operate.
 *
 * @param[in] eltsize
 *          PLASMA internal value which refer to the size of the precision.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
/***************************************************************************//**
 *  TYPE 1-BDL Householder
 *  add -1 because of C
 ******************************************************************************/
int
CORE_zgbrce(PLASMA_enum uplo, int N,
            PLASMA_desc *A,
            PLASMA_Complex64_t *V,
            PLASMA_Complex64_t *TAU,
            int st,
            int ed,
            int eltsize)
{
    int    NB, J1, J2, J3, KDM2, len, pt;
    int    len1, len2, t1ed, t2st;
    int    i;
    static PLASMA_Complex64_t zzero = 0.0;
    PLASMA_desc vA=*A;

    /* Check input arguments */
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (ed <= st) {
        coreblas_error(6, "Illegal value of st and ed (internal)");
        return -6;
    }

    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    NB = A->mb;
    KDM2 = A->mb-2;
    if( uplo == PlasmaLower ){
        /* ========================
         *       LOWER CASE
         * ========================*/
        for (i = ed; i >= st+1 ; i--){
            /* apply Householder from the right. and create newnnz outside the band if J3 < N */
            J1  = ed+1;
            J2  = min((i+1+KDM2), N);
            J3  = min((J2+1), N);
            len = J3-J1+1;
            if(J3>J2)*A(J3,(i-1))=zzero;/* could be removed because A is supposed to be band.*/

            t1ed  = (J3/NB)*NB;
            t2st  = max(t1ed+1,J1);
            len1  = t1ed-J1+1;
            len2  = J3-t2st+1;
            if(len1>0)CORE_zlarfx2(PlasmaRight, len1, conj(*V(i)), conj(*TAU(i)), A(J1,  i-1), ELTLDD(vA, J1)  , A(J1  , i), ELTLDD(vA, J1)  );
            if(len2>0)CORE_zlarfx2(PlasmaRight, len2, conj(*V(i)), conj(*TAU(i)), A(t2st,i-1), ELTLDD(vA, t2st), A(t2st, i), ELTLDD(vA, t2st));
            len    = J3-J2;
            if(len>0){
                /* generate Householder to annihilate a(j+kd,j-1) within the band */
                *V(J3)         = *A(J3,(i-1));
                *A(J3,(i-1))   = 0.0;
                LAPACKE_zlarfg_work( 2, A(J2,(i-1)), V(J3), 1, TAU(J3));
            }
        }
        /* APPLY LEFT ON THE REMAINING ELEMENT OF KERNEL 2 */
        for (i = ed; i >= st+1 ; i--){
            J2  = min((i+1+KDM2), N);
            J3  = min((J2+1), N);
            len    = J3-J2;
            if(len>0){
                pt    = J2;
                J1    = i;
                J2    = min(ed,N);
                t1ed  = (J2/NB)*NB;
                t2st  = max(t1ed+1,J1);
                len1  = t1ed-J1+1;
                len2  = J2-t2st+1;
                if(len1>0)CORE_zlarfx2(PlasmaLeft, len1 , *V(J3), conj(*TAU(J3)), A(pt, i   ), ELTLDD(vA, pt),  A((pt+1),  i  ), ELTLDD(vA, pt+1) );
                if(len2>0)CORE_zlarfx2(PlasmaLeft, len2 , *V(J3), conj(*TAU(J3)), A(pt, t2st), ELTLDD(vA, pt),  A((pt+1), t2st), ELTLDD(vA, pt+1) );
            }
        }
    } else {
        /* ========================
         *       UPPER CASE
         * ========================*/
        for (i = ed; i >= st+1 ; i--){
            /* apply Householder from the right. and create newnnz outside the band if J3 < N */
            J1  = ed+1;
            J2  = min((i+1+KDM2), N);
            J3  = min((J2+1), N);
            len = J3-J1+1;
            if(J3>J2)*A((i-1), J3)=zzero;

            t1ed  = (J3/NB)*NB;
            t2st  = max(t1ed+1,J1);
            len1  = t1ed-J1+1;
            len2  = J3-t2st+1;
            if(len1>0)CORE_zlarfx2(PlasmaLeft, len1 , *V(i), conj(*TAU(i)), A(i-1, J1  ), ELTLDD(vA, i-1),  A(i,  J1 ), ELTLDD(vA, i) );
            if(len2>0)CORE_zlarfx2(PlasmaLeft, len2 , *V(i), conj(*TAU(i)), A(i-1, t2st), ELTLDD(vA, i-1),  A(i, t2st), ELTLDD(vA, i) );
            /* if nonzero element a(j+kd,j-1) has been created outside the band (if index < N) then eliminate it. */
            len    = J3-J2;
            if(len>0){
                /* generate Householder to annihilate a(j+kd,j-1) within the band */
                *V(J3)         = *A(i-1, J3);
                *A(i-1, J3)  = 0.0;
                LAPACKE_zlarfg_work( 2, A(i-1, J2), V(J3), 1, TAU(J3));
            }
        }
        /* APPLY RIGHT ON THE REMAINING ELEMENT OF KERNEL 2 */
        for (i = ed; i >= st+1 ; i--){
            /* find if there was a nnz created. if yes apply right else nothing to be done. */
            J2  = min((i+1+KDM2), N);
            J3  = min((J2+1), N);
            len    = J3-J2;
            if(len>0){
                pt    = J2;
                J1    = i;
                J2    = min(ed,N);
                t1ed  = (J2/NB)*NB;
                t2st  = max(t1ed+1,J1);
                len1  = t1ed-J1+1;
                len2  = J2-t2st+1;
                if(len1>0)CORE_zlarfx2(PlasmaRight, len1 , conj(*V(J3)), conj(*TAU(J3)), A(i   , pt), ELTLDD(vA, i),     A(i,    pt+1), ELTLDD(vA, i) );
                if(len2>0)CORE_zlarfx2(PlasmaRight, len2 , conj(*V(J3)), conj(*TAU(J3)), A(t2st, pt), ELTLDD(vA, t2st),  A(t2st, pt+1), ELTLDD(vA, t2st) );
            }
        }
    } /* end of else for the upper case */

    return PLASMA_SUCCESS;
}
