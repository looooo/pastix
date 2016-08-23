/**
 *
 * @file coeftab_zdiff.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "lapacke.h"
#include "coeftab.h"
#include "pastix_zcores.h"

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
    assert( !(cblkA->cblktype & CBLK_SPLIT) ); /* Not yet implemented */

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

    if ( resU > 10 ) {
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

