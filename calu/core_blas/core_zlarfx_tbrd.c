/**
 *
 * @file core_zlarfx_tbrd.c
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

/***************************************************************************//**
 * @ingroup CORE_PLASMA_Complex64_t
 *******************************************************************************
 *
 *  Purpose
 *  =======
 *
 *  CORE_zlarfx2 applies a complex elementary reflector H to a complex m by n
 *  matrix C, from either the left or the right. H is represented in the
 *  form
 *
 *        H = I - tau * v * v'
 *
 *  where tau is a complex scalar and v is a complex vector.
 *
 *  If tau = 0, then H is taken to be the unit matrix
 *
 *  This version uses inline code if H has order < 11.
 *
 *  Arguments
 *  =========
 *
 *  @param[in] side
 *          @arg PlasmaLeft : form  H * C
 *          @arg PlasmaRight: form  C * H
 *
 * @param[in] N
 *          The number of columns of C1 and C2 if side = PlasmaLeft.
 *          The number of rows of C1 and C2 if side = PlasmaRight.
 *
 * @param[in] V
 *          The double complex V in the representation of H.
 *
 * @param[in] TAU
 *          The value tau in the representation of H.
 *
 * @param[in,out] C1
 *          dimension (LDC1,N), if side = PlasmaLeft
 *          dimension (LDC1,1), if side = PlasmaRight
 *          On entry, the m by n matrix C1.
 *          On exit, C1 is overwritten by the matrix H * C1 if SIDE = PlasmaLeft,
 *          or C1 * H if SIDE = PlasmaRight.
 *
 * @param[in] LDC1
 *          The leading dimension of the array C1.
 *          LDC1 >= max(1,N), if side == PlasmaRight.
 *          LDC1 >= 1, otherwise.
 *
 * @param[in,out] C2
 *          dimension (LDC2,N), if side = PlasmaLeft
 *          dimension (LDC2,1), if side = PlasmaRight
 *          On entry, the m by n matrix C2.
 *          On exit, C2 is overwritten by the matrix H * C2 if SIDE = PlasmaLeft,
 *          or C2 * H if SIDE = PlasmaRight.
 *
 * @param[in] LDC2
 *          The leading dimension of the array C2.
 *          LDC2 >= max(1,N), if side == PlasmaRight.
 *          LDC2 >= 1, otherwise.
 *
 *  =====================================================================
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
CORE_zlarfx2(PLASMA_enum side, int N,
             PLASMA_Complex64_t V,
             PLASMA_Complex64_t TAU,
             PLASMA_Complex64_t *C1, int LDC1,
             PLASMA_Complex64_t *C2, int LDC2)
{
    PLASMA_Complex64_t V2, T2, SUM;
    int j;

    if (TAU == (PLASMA_Complex64_t)0.0)
        return PLASMA_SUCCESS;

    /*
     * Special code for 2 x 2 Householder where V1 = I
    */
    if(side==PlasmaLeft){
        V2 = conj(V);
        T2 = TAU*conj(V2);
        for (j = 0; j < N ; j++, C1+=LDC1 ) {
            SUM = *C1 + V2 * (*C2);
            *C1 = *C1 - SUM*TAU;
            *C2 = *C2 - SUM*T2;
            C2 += LDC2;
        }
    }
    else {
        V2 = V;
        T2 = TAU*conj(V2);
        for (j = 0; j < N ; j++, C1++){
            SUM = *C1 + V2 * (*C2);
            *C1 = *C1 - SUM*TAU;
            *C2 = *C2 - SUM*T2;
            C2++;
        }
    }

    return PLASMA_SUCCESS;
}


/***************************************************************************//**
 * @ingroup CORE_PLASMA_Complex64_t
 *******************************************************************************
 *
 *  Purpose
 *  =======
 *
 *  CORE_zlarfx2c applies a complex elementary reflector H to a diagonal corner
 *  C=[C1, C2, C3], from both the left and the right side. C = H * C * H.
 *  It is used in the case of Hermetian.
 *  If PlasmaLower, a left apply is followed by a right apply.
 *  If PlasmaUpper, a right apply is followed by a left apply.
 *  H is represented in the form
 *
 *  This routine is a special code for a corner C diagonal block  C1
 *                                                                C2  C3
 *
 *
 *        H = I - tau * v * v'
 *
 *  where tau is a complex scalar and v is a complex vector.
 *
 *  If tau = 0, then H is taken to be the unit matrix
 *
 *  This version uses inline code if H has order < 11.
 *
 *  Arguments
 *  =========
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] V
 *          The double complex V in the representation of H.
 *
 * @param[in] TAU
 *          The value tau in the representation of H.
 *
 * @param[in,out] C1
 *          On entry, the element C1.
 *          On exit, C1 is overwritten by the result H * C * H.
 *
 * @param[in,out] C2
 *          On entry, the element C2.
 *          On exit, C2 is overwritten by the result H * C * H.
  *
 * @param[in,out] C3
 *          On entry, the element C3.
 *          On exit, C3 is overwritten by the result H * C * H.
 *
 *  =====================================================================
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
CORE_zlarfx2c(PLASMA_enum uplo,
              PLASMA_Complex64_t V,
              PLASMA_Complex64_t TAU,
              PLASMA_Complex64_t *C1,
              PLASMA_Complex64_t *C2,
              PLASMA_Complex64_t *C3)
{
    PLASMA_Complex64_t T2, SUM, TEMP;

    /* Quick return */
    if (TAU == (PLASMA_Complex64_t)0.0)
        return PLASMA_SUCCESS;

    /*
     *        Special code for a diagonal block  C1
     *                                           C2  C3
     */
    if(uplo==PlasmaLower) {
        /*
         *  Do the corner Left then Right (used for the lower case
         *  tridiag) L and R for the 2x2 corner
         *             C(N-1, N-1)  C(N-1,N)        C1  TEMP
         *             C(N  , N-1)  C(N  ,N)        C2  C3
         *  For Left : use conj(TAU) and V.
         *  For Right: nothing, keep TAU and V.
         *  Left 1 ==> C1
         *             C2
         */
        TEMP = conj(*C2); /*  copy C2 here before modifying it. */
        T2   = conj(TAU) * V;
        SUM  = *C1 + conj(V) * (*C2);
        *C1  = *C1 - SUM * conj(TAU);
        *C2  = *C2 - SUM * T2;
        /*  Left 2 ==> TEMP */
        /*             C3 */
        SUM  = TEMP + conj(V) * (*C3);
        TEMP = TEMP - SUM * conj(TAU);
        *C3  = *C3  - SUM * T2;
        /*  Right 1 ==>  C1 TEMP.  NB: no need to compute corner (2,2)=TEMP */
        T2  = TAU * conj(V);
        SUM = *C1 + V*TEMP;
        *C1 = *C1 - SUM*TAU;
        /*  Right 2 ==> C2 C3 */
        SUM = *C2   + V*(*C3);
        *C2 = *C2   - SUM*TAU;
        *C3 = *C3   - SUM*T2;
    }
    else {
        /*
         * Do the corner Right then Left (used for the upper case tridiag)
         *             C(N-1, N-1)  C(N-1,N)        C1    C2
         *             C(N  , N-1)  C(N  ,N)        TEMP  C3
         *  For Left : use TAU       and conj(V).
         *  For Right: use conj(TAU) and conj(V).
         *  Right 1 ==> C1 C2
         */
        V    = conj(V);
        TEMP = conj(*C2); /*  copy C2 here before modifying it. */
        T2   = conj(TAU) * conj(V);
        SUM  = *C1 + V   * (*C2);
        *C1  = *C1 - SUM * conj(TAU);
        *C2  = *C2 - SUM * T2;
        /*  Right 2 ==> TEMP C3 */
        SUM  = TEMP + V   * (*C3);
        TEMP = TEMP - SUM * conj(TAU);
        *C3  = *C3  - SUM * T2;
        /*  Left 1 ==> C1 */
        /*             TEMP. NB: no need to compute corner (2,1)=TEMP */
        T2  = TAU * V;
        SUM = *C1 + conj(V) * TEMP;
        *C1 = *C1 - SUM * TAU;
        /*  Left 2 ==> C2 */
        /*             C3 */
        SUM = *C2 + conj(V) * (*C3);
        *C2 = *C2 - SUM * TAU;
        *C3 = *C3 - SUM * T2;
    }

    return PLASMA_SUCCESS;
}


/***************************************************************************//**
 * @ingroup CORE_PLASMA_Complex64_t
 *******************************************************************************
 *
 *  Purpose
 *  =======
 *
 *  CORE_zlarfx2c applies a complex elementary reflector H to a diagonal corner
 *  C=[C1, C2, C3], from both the left and the right side. C = H * C * H.
 *  It is used in the case of general matrices, where it create a nnz at the
 *  NEW_NNZ position, then it eliminate it and update the reflector V and TAU.
 *  If PlasmaLower, a left apply is followed by a right apply.
 *  If PlasmaUpper, a right apply is followed by a left apply.
 *  H is represented in the form
 *
 *  This routine is a special code for a corner C diagonal block  C1  NEW_NNZ
 *                                                                C2  C3
 *
 *
 *        H = I - tau * v * v'
 *
 *  where tau is a complex scalar and v is a complex vector.
 *
 *  If tau = 0, then H is taken to be the unit matrix
 *
 *  This version uses inline code if H has order < 11.
 *
 *  Arguments
 *  =========
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in, out] V
 *          On entry, the double complex V in the representation of H.
 *          On exit, the double complex V in the representation of H,
 *          updated by the elimination of the NEW_NNZ created by the
 *          left apply in case of PlasmaLower or the right apply in
 *          case of PlasmaUpper.
 *
 * @param[in,out] TAU
 *          On entry, the value tau in the representation of H.
 *          On exit, the value tau in the representation of H,
 *          updated by the elimination of the NEW_NNZ created by the
 *          left apply in case of PlasmaLower or the right apply in
 *          case of PlasmaUpper.
 *
 * @param[in,out] C1
 *          On entry, the element C1.
 *          On exit, C1 is overwritten by the result H * C * H.
 *
 * @param[in,out] C2
 *          On entry, the element C2.
 *          On exit, C2 is overwritten by the result H * C * H.
  *
 * @param[in,out] C3
 *          On entry, the element C3.
 *          On exit, C3 is overwritten by the result H * C * H.
 *
 *  =====================================================================
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
CORE_zlarfx2ce(PLASMA_enum uplo,
               PLASMA_Complex64_t *V,
               PLASMA_Complex64_t *TAU,
               PLASMA_Complex64_t *C1,
               PLASMA_Complex64_t *C2,
               PLASMA_Complex64_t *C3)
{
    PLASMA_Complex64_t T2, SUM, TEMP, VIN, TAUIN;

    /* Quick return */
    if (*TAU == (PLASMA_Complex64_t)0.0)
        return PLASMA_SUCCESS;

    /*
     *        Special code for a diagonal block  C1
     *                                           C2  C3
     */
    if(uplo==PlasmaLower){
        /*
         *  Do the corner for the lower case BIDIAG ==> Left then will
         *  create a new nnz. eliminate it and modify V TAU and then
         *  Right L and R for the 2x2 corner
         *             C(N-1, N-1)  C(N-1,N)        C1  TEMP
         *             C(N  , N-1)  C(N  ,N)        C2  C3
         */
        VIN   = *V;
        TAUIN = conj(*TAU);
        /*  Left 1 ==> C1 */
        /*             C2 */
        VIN  = conj(VIN);
        T2   = TAUIN * conj(VIN);
        SUM  = *C1   + VIN*(*C2);
        *C1  = *C1   - SUM*TAUIN;
        *C2  = *C2   - SUM*T2;
        /*   new nnz at TEMP and update C3 */
        SUM  =        VIN * (*C3);
        TEMP =      - SUM * TAUIN;
        *C3  = *C3  - SUM * T2;
        /*  generate Householder to annihilate the nonzero created at TEMP */
        *V    = TEMP;
        LAPACKE_zlarfg_work( 2, C1, V, 1, TAU);
        VIN   = conj(*V);
        TAUIN = conj(*TAU);
        /*  Right 1 ==> C2 C3 */
        /* VIN     = VIN */
        T2  = TAUIN * conj(VIN);
        SUM = *C2   + VIN*(*C3);
        *C2 = *C2   - SUM*TAUIN;
        *C3 = *C3   - SUM*T2;
    }else if(uplo==PlasmaUpper){
        /*
         * Do the corner for the upper case BIDIAG ==> Right then will
        *  create a new nnz. eliminate it and modify V TAU and then
        *  Left
        *             C(N-1, N-1)  C(N-1,N)        C1    C2
        *             C(N  , N-1)  C(N  ,N)        TEMP  C3
        *  For Left : use conj(TAU) and V.
        *  For Right: use conj(TAU) and conj(V) as input.
        */
        VIN   = conj(*V);
        TAUIN = conj(*TAU);
        /*  Right 1 ==> C1 C2 */
        /* VIN     = VIN */
        T2   = TAUIN*conj(VIN);
        SUM  = *C1  + VIN*(*C2);
        *C1  = *C1  - SUM*TAUIN;
        *C2  = *C2  - SUM*T2;
        /*   new nnz at TEMP and update C3 */
        SUM  =        VIN * (*C3);
        TEMP =      - SUM * TAUIN;
        *C3  = *C3  - SUM * T2;
        /*  generate Householder to annihilate the nonzero created at TEMP */
        *V    = TEMP;
        LAPACKE_zlarfg_work( 2, C1, V, 1, TAU);
        VIN   = *V;
        TAUIN = conj(*TAU);
        /*  apply from the Left using the NEW V TAU to the remaining 2 elements [C2 C3] */
        /*  Left 2 ==> C2 */
        /*             C3 */
        VIN = conj(VIN);
        T2  = TAUIN*conj(VIN);
        SUM = *C2 + VIN*(*C3);
        *C2 = *C2 - SUM*TAUIN;
        *C3 = *C3 - SUM*T2;
    }
    return PLASMA_SUCCESS;
}

