/**
 *
 * @file coeftab_zdiff.c
 *
 * Precision dependent routines to differentiate two solver matrix structures when debuging.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2017-04-28
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "lapacke.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

/**
 *******************************************************************************
 *
 * @brief Compare two column blocks in full-rank format.
 *
 * The second cblk is overwritten by the difference of the two column blocks.
 * The frobenius norm of the difference is computed and the functions returns 0
 * if the result:
 *      || B - A || / ( || A || * eps )
 *
 * is below 10. Otherwise, an error message is printed and 1 is returned.
 *
 *******************************************************************************
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix that matches the A matrix in
 *          stucture.
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          (B-A).
 *
 *******************************************************************************
 *
 * @return 0 if the test is passed, >= 0 otherwise.
 *
 *******************************************************************************/
int
coeftab_zdiffcblk( const SolverCblk *cblkA,
                   SolverCblk *cblkB )
{
    pastix_complex64_t *lcoefA = cblkA->lcoeftab;
    pastix_complex64_t *ucoefA = cblkA->ucoeftab;
    pastix_complex64_t *lcoefB = cblkB->lcoeftab;
    pastix_complex64_t *ucoefB = cblkB->ucoeftab;
    pastix_int_t        ncols  = cblk_colnbr( cblkA );
    pastix_int_t        stride = cblkA->stride;
    double normL, normU, normfAL, normfAU, normcAL, normcAU, resL, resU, eps;
    int factoLU = (ucoefA == NULL) ? 0 : 1;
    int rc = 0;

    assert( ncols  == cblk_colnbr( cblkB ) );
    assert( stride == cblkB->stride );
    assert( !(cblkA->cblktype & CBLK_LAYOUT_2D) ); /* Not yet implemented */

    eps = LAPACKE_dlamch_work( 'e' );

    normfAL = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                   lcoefA, stride, NULL );
    normcAL = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                   lcoefB, stride, NULL );
    core_zgeadd( PastixNoTrans, stride, ncols,
                 -1., lcoefA, stride,
                  1., lcoefB, stride );

    normL = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', stride, ncols,
                            lcoefB, stride );
    resL = (normfAL == 0.) ? 0. : (normL / (normfAL * eps));

    if ( factoLU ) {
        normfAU = LAPACKE_zlange( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                  ucoefA, stride );
        normcAU = LAPACKE_zlange( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                  ucoefB, stride );
        core_zgeadd( PastixNoTrans, stride, ncols,
                     -1., ucoefA, stride,
                      1., ucoefB, stride );
        normU  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', stride, ncols,
                                 ucoefB, stride );
        resU = (normfAU == 0.) ? 0. : (normU / (normfAU * eps));
    }
    else {
        resU = 0.;
    }

    if ( resL > 10 ) {
        fprintf(stderr, "KO on L: ||full(A)||_f=%e, ||comp(A)||_f=%e, ||comp(A)-full(A)||_0=%e, ||comp(A)-full(A)||_0 / (||full(A)||_2 * eps)=%e\n",
                normfAL, normcAL, normL, resL );
        rc++;
    }
    /* else { */
    /*     fprintf(stderr, "Ok on L: ||full(A)||_f=%e, ||comp(A)||_f=%e, ||comp(A)-full(A)||_0=%e, ||comp(A)-full(A)||_0 / (||full(A)||_2 * eps)=%e\n", */
    /*             normfAL, normcAL, normL, resL ); */
    /* } */

    if ( factoLU && (resU > 10) ) {
        fprintf(stderr, "KO on U: ||full(A)||_f=%e, ||comp(A)||_f=%e, ||comp(A)-full(A)||_0=%e, ||comp(A)-full(A)||_0 / (||full(A)||_2 * eps)=%e\n",
                normfAU, normcAU, normU, resU );
        rc++;
    }
    /* else { */
    /*     fprintf(stderr, "Ok on U: ||full(A)||_f=%e, ||comp(A)||_f=%e, ||comp(A)-full(A)||_0=%e, ||comp(A)-full(A)||_0 / (||full(A)||_2 * eps)=%e\n", */
    /*             normfAU, normcAU, normU, resU ); */
    /* } */

    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Compare two solver matrices full-rank format.
 *
 * The second solver matrix is overwritten by the difference of the two
 * matrices.  The frobenius norm of the difference of each clolumn block is
 * computed and the functions returns 0 if the result for all the column blocks
 * of:
 *      || B_k - A_k || / ( || A_k || * eps )
 *
 * is below 10. Otherwise, an error message is printed and 1 is returned.
 *
 *******************************************************************************
 *
 * @param[in] solvA
 *          The solver matrix A.
 *
 * @param[inout] cblkB
 *          The solver matrix B.
 *          On exit, B coefficient arrays are overwritten by the result of
 *          (B-A).
 *
 *******************************************************************************
 *
 * @return 0 if the test is passed, >= 0 otherwise.
 *
 *******************************************************************************/
int
coeftab_zdiff( const SolverMatrix *solvA, SolverMatrix *solvB )
{
    SolverCblk *cblkA = solvA->cblktab;
    SolverCblk *cblkB = solvB->cblktab;
    pastix_int_t cblknum;
    int rc       = 0;
    int saved_rc = 0;

    for(cblknum=0; cblknum<solvA->cblknbr; cblknum++, cblkA++, cblkB++) {
        rc += coeftab_zdiffcblk( cblkA, cblkB );
        if ( rc != saved_rc ){
            fprintf(stderr, "CBLK %ld was not correctly compressed\n", cblknum);
            saved_rc = rc;
        }
    }

    return rc;
}
