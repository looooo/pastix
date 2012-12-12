/**
 *
 * @file pzhbcpy_t2bl.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n)   BLKADDR(A, PLASMA_Complex64_t, m, n)
#define AB(m,n)   &(AB[(n) + LDAB*((j)*nb) ])
/***************************************************************************//**
 *  Parallel copy of a band matrix from full NxN tile storage to band storage (NxLDAB).
 *  As this function is internal and the space is the same for either Lower or Upper so 
 *  ALWAYS it convert to Lower band and then the bulge chasing will 
 *  always work with a Lower band matrix
 **/
void plasma_pzhbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, 
                               PLASMA_Complex64_t *AB, int LDAB, 
                               PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int j;
    int ldaj, ldx;
    int nb=A.mb;
    int tempjm, tempjn;
    int minmn = min(A.mt, A.nt);

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /*=============================================
     * NOTE :
     * this function transform the Lower/Upper Tile 
     * band matrix to LOWER Band storage matrix.
     * For Lower it copy it directly.
     * For Upper it conjtransposed during the copy.
     *=============================================*/
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    ldx = LDAB-1;
    
    /* copy Lower to Lower */
    if ( uplo == PlasmaLower ) {
       for (j = 0; j < minmn; j++){
           tempjm = j == A.mt-1 ? A.m-j*A.mb : A.mb;
           tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
           ldaj = BLKLDD(A, j);
           QUARK_CORE_zlacpy(
               plasma->quark, &task_flags,
               PlasmaLower, tempjm, tempjn, nb,
               A(j, j), ldaj, AB(j, 0), ldx);

           if( j<minmn-1 ){
               tempjm = (j+1) == A.mt-1 ? A.m-(j+1)*A.mb : A.mb;
               tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
               ldaj = BLKLDD(A, j+1);
               QUARK_CORE_zlacpy(
                   plasma->quark, &task_flags,
                   PlasmaUpper, tempjm, tempjn, nb,
                   A(j+1, j), ldaj, AB(j, nb), ldx);
           }
       }
    }
    /* conjtranspose Upper when copying it to Lower */
    else if ( uplo == PlasmaUpper ) {
       for (j = 0; j < minmn; j++){
           tempjm = j == A.mt-1 ? A.m-j*A.mb : A.mb;
           tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
           ldaj = BLKLDD(A, j);
           QUARK_CORE_zlatro(
               plasma->quark, &task_flags,
               PlasmaUpper, PlasmaConjTrans, tempjm, tempjn, nb,
               A(j, j), ldaj, AB(j, 0), ldx);

           if(j<minmn-1){
               tempjm = j == A.mt-1 ? A.m-j*A.mb : A.mb;
               tempjn = (j+1) == A.nt-1 ? A.n-(j+1)*A.nb : A.nb;
               ldaj = BLKLDD(A, j);
               QUARK_CORE_zlatro(
                   plasma->quark, &task_flags,
                   PlasmaLower, PlasmaConjTrans, tempjm, tempjn, nb,
                   A(j, j+1), ldaj, AB(j, nb), ldx);
           }
       }
    }
}
#undef AB
#undef A
