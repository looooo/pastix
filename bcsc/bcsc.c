/**
 *
 * @file bcsc.c
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include "spm.h"
#include "solver.h"
#include "bcsc.h"

#include "z_bcsc.h"
#include "c_bcsc.h"
#include "d_bcsc.h"
#include "s_bcsc.h"

static inline pastix_int_t
bcsc_init_coltab( const SolverMatrix  *solvmtx,
                  const pastix_int_t  *newcoltab,
                        pastix_int_t   dof,
                        pastix_bcsc_t *bcsc )
{
    bcsc_format_t *blockcol;
    pastix_int_t index, iter, idxcol, nodeidx, colsize;

    bcsc->cscfnbr = solvmtx->cblknbr;
    MALLOC_INTERN( bcsc->cscftab, bcsc->cscfnbr, bcsc_format_t );

    idxcol = 0;
    blockcol = bcsc->cscftab;
    for (index=0; index<bcsc->cscfnbr; index++, blockcol++)
    {
        pastix_int_t fcolnum = solvmtx->cblktab[index].fcolnum;
        pastix_int_t lcolnum = solvmtx->cblktab[index].lcolnum;

        blockcol->colnbr = (lcolnum - fcolnum + 1);
        MALLOC_INTERN( blockcol->coltab, blockcol->colnbr + 1, pastix_int_t );

        /* Works only for DoF constant */
        assert( fcolnum % dof == 0 );

        blockcol->coltab[0] = idxcol;
        for (iter=0; iter < blockcol->colnbr; iter++)
        {
            nodeidx = ( fcolnum + (iter-iter%dof) ) / dof;

            /* if (g2l != NULL && */
            /*     iter != CSC_COLNBR(thecsc,index) && */
            /*     !COL_IS_LOCAL(g2l[ord->peritab[nodeidx]])) */
            /* { */
            /*     errorPrint("Columns in internal CSCD must be in given CSCD"); */
            /* } */

            colsize = (newcoltab[nodeidx+1] - newcoltab[nodeidx]) * dof;
            blockcol->coltab[iter+1] = blockcol->coltab[iter] + colsize;
        }

        idxcol = blockcol->coltab[blockcol->colnbr];
    }

    MALLOC_INTERN( bcsc->rowtab, idxcol, pastix_int_t);
    MALLOC_INTERN( bcsc->Lvalues, idxcol * pastix_size_of( bcsc->flttype ), char );
    bcsc->Uvalues = NULL;

    return idxcol;
}

void
bcsc_restore_coltab( pastix_bcsc_t *bcsc )
{
    bcsc_format_t *blockcol;
    pastix_int_t index, iter, idxcol, idxcoltmp;

    idxcol = 0;
    blockcol = bcsc->cscftab;
    for (index=0; index<bcsc->cscfnbr; index++, blockcol++)
    {
        for (iter=0; iter <= blockcol->colnbr; iter++)
        {
            idxcoltmp = blockcol->coltab[iter];
            blockcol->coltab[iter] = idxcol;
            idxcol = idxcoltmp;
        }
    }
    return;
}

pastix_int_t
bcsc_init_centralized_coltab( const pastix_spm_t   *spm,
                              const pastix_order_t *ord,
                              const SolverMatrix   *solvmtx,
                                    pastix_bcsc_t  *bcsc )
{
    pastix_int_t  valuesize, baseval;
    pastix_int_t *globcol  = NULL;
    pastix_int_t *colptr = spm->colptr;
    pastix_int_t *rowptr = spm->rowptr;
    int dof = spm->dof;
    int sym = (spm->mtxtype == PastixSymmetric) || (spm->mtxtype == PastixHermitian);

    bcsc->mtxtype = spm->mtxtype;
    baseval = spm->colptr[0];

    /**
     * Allocate and initialize globcol that contains the number of elements in
     * each column of the input matrix
     * Globcol is equivalent to the classic colptr for the internal blocked
     * csc. The blocked csc integrate the perumtation computed within order
     * structure.
     */
    MALLOC_INTERN( globcol, spm->gN+1, pastix_int_t );
    memset( globcol, 0, (spm->gN+1) * sizeof(pastix_int_t) );

    assert( spm->loc2glob == NULL );

    {
        pastix_int_t itercol, newcol;

        for (itercol=0; itercol<spm->gN; itercol++)
        {
            pastix_int_t frow = colptr[itercol]   - baseval;
            pastix_int_t lrow = colptr[itercol+1] - baseval;
            newcol = ord->permtab[itercol];
            globcol[newcol] += lrow - frow;

            assert( (lrow - frow) >= 0 );
            if (sym) {
                pastix_int_t iterrow, newrow;

                for (iterrow=frow; iterrow<lrow; iterrow++)
                {
                    pastix_int_t tmprow = rowptr[iterrow] - baseval;
                    if (tmprow != itercol) {
                        newrow = ord->permtab[tmprow];
                        globcol[newrow]++;
                    }
                }
            }
        }

        /* Compute displacements to update the colptr array */
        {
            pastix_int_t tmp, idx;

            idx = 0;
            for (itercol=0; itercol<=spm->gN; itercol++)
            {
                tmp = globcol[itercol];
                globcol[itercol] = idx;
                idx += tmp;
            }
        }
    }

    valuesize = bcsc_init_coltab( solvmtx, globcol, dof, bcsc );
    memFree_null( globcol );

    return valuesize;
}

void
bcscInitCentralized( const pastix_spm_t   *spm,
                     const pastix_order_t *ord,
                     const SolverMatrix   *solvmtx,
                           pastix_int_t    initAt,
                           pastix_bcsc_t  *bcsc )
{
    pastix_int_t  itercol, itercblk;
    pastix_int_t  cblknbr = solvmtx->cblknbr;
    pastix_int_t  eltnbr  = spm->gN * spm->dof + 1;
    pastix_int_t *col2cblk = NULL;

    bcsc->mtxtype = spm->mtxtype;
    bcsc->flttype = spm->flttype;
    bcsc->gN      = spm->gN;
    bcsc->n       = spm->n;

    assert( spm->loc2glob == NULL );

    /**
     * Initialize the col2cblk array. col2cblk[i] contains the cblk index of the
     * i-th column. col2cblk[i] = -1 if not local.
     */
    {
        SolverCblk *cblk = solvmtx->cblktab;

        MALLOC_INTERN( col2cblk, eltnbr, pastix_int_t );
        for (itercol=0; itercol<eltnbr; itercol++)
        {
            col2cblk[itercol] = -1;
        }

        for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
        {
            for (itercol  = cblk->fcolnum;
                 itercol <= cblk->lcolnum;
                 itercol++ )
            {
                col2cblk[itercol] = itercblk;
            }
        }
    }

    /**
     * Fill in the lower triangular part of the blocked csc with values and
     * rows. The upper triangular part is done later if required through LU
     * factorization.
     */
    switch( spm->flttype ) {
    case PastixFloat:
        s_bcscInitCentralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixDouble:
        d_bcscInitCentralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixComplex32:
        c_bcscInitCentralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixComplex64:
        z_bcscInitCentralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixPattern:
    default:
        fprintf(stderr, "bcscInitCentralized: Error unknown floating type for input spm\n");
    }

    memFree_null(col2cblk);
}

double
bcscInit( const pastix_spm_t   *spm,
          const pastix_order_t *ord,
          const SolverMatrix   *solvmtx,
                pastix_int_t    initAt,
                pastix_bcsc_t  *bcsc )
{
    assert( ord->baseval == 0 );
    assert( ord->vertnbr == spm->n );

    double time = 0.;
    clockStart(time);

    if ( spm->loc2glob == NULL ) {
        bcscInitCentralized( spm, ord, solvmtx, initAt, bcsc );
    }
    else {
        fprintf(stderr, "bcscInit: Distributed SPM not yet supported");
    }

    clockStop(time);
    return time;
}

/******************************************************************************
 * Function: z_CscExit                                                        *
 ******************************************************************************
 *                                                                            *
 * Free the internal CSCd structure.                                          *
 *                                                                            *
 * Parameters:                                                                *
 *   thecsc - Internal CSCd to free.                                          *
 *                                                                            *
 ******************************************************************************/
void
bcscExit( pastix_bcsc_t *bcsc )
{
    if ( bcsc->cscftab != NULL )
    {
        pastix_int_t itercscf;
        for (itercscf = 0; itercscf < bcsc->cscfnbr; itercscf++ )
        {
            memFree_null( bcsc->cscftab[itercscf].coltab );
        }

        memFree_null( bcsc->cscftab );
        memFree_null( bcsc->rowtab );

        if ( (bcsc->Uvalues != NULL) &&
             (bcsc->Uvalues != bcsc->Lvalues) ) {
            memFree_null( bcsc->Uvalues );
        }

        memFree_null( bcsc->Lvalues );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * bcscMatVec - Compute the matrix-vector product
 *         y = alpha * op(A) * x + beta * y,
 * where A is given in the bcsc format, x and y are two vectors of size n, and
 * alpha and beta are two scalars.
 * The op function is specified by the trans parameter and performs the
 * operation as follows:
 *              trans = PastixNoTrans   y := alpha*A       *x + beta*y
 *              trans = PastixTrans     y := alpha*A'      *x + beta*y
 *              trans = PastixConjTrans y := alpha*conj(A')*x + beta*y
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A from the bcsc is transposed, not
 *          transposed or conjugate transposed:
 *            = PastixNoTrans:   A is not transposed;
 *            = PastixTrans:     A is transposed;
 *            = PastixConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] bcsc
 *          The bcsc structure describing the matrix A.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
bcscMatVec(       int            trans,
            const void          *alpha,
            const pastix_bcsc_t *bcsc,
            const void          *x,
            const void          *beta,
                  void          *y )
{
    switch (bcsc->flttype) {
    case PastixFloat:
        return s_bcscGemv( trans, *((const float*)alpha), bcsc, (const float*)x, *((const float*)beta), (float*)y );
    case PastixComplex32:
        return c_bcscGemv( trans, *((const pastix_complex32_t*)alpha), bcsc, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
    case PastixComplex64:
        return z_bcscGemv( trans, *((const pastix_complex64_t*)alpha), bcsc, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
    case PastixDouble:
    default:
        return d_bcscGemv( trans, *((const double*)alpha), bcsc, (const double*)x, *((const double*)beta), (double*)y );
    }
}

int
bcscApplyPerm( pastix_bcsc_t *bcsc,
               pastix_int_t   n,
               void          *b,
               pastix_int_t   ldb,
               pastix_int_t  *perm )
{
    switch( bcsc->flttype ) {
    case PastixComplex64:
        if( PASTIX_SUCCESS != z_bcscApplyPerm( bcsc->gN,
                                               n,
                                               b,
                                               ldb,
                                               perm ))
        {
            return PASTIX_ERR_BADPARAMETER;
        }
        break;

    case PastixComplex32:
        if( PASTIX_SUCCESS != c_bcscApplyPerm( bcsc->gN,
                                               n,
                                               b,
                                               ldb,
                                               perm ))
        {
            return PASTIX_ERR_BADPARAMETER;
        }
        break;

    case PastixFloat:
        if( PASTIX_SUCCESS != s_bcscApplyPerm( bcsc->gN,
                                               n,
                                               b,
                                               ldb,
                                               perm ))
        {
            return PASTIX_ERR_BADPARAMETER;
        }
        break;

    case PastixDouble:
    default:
        if( PASTIX_SUCCESS != d_bcscApplyPerm( bcsc->gN,
                                               n,
                                               b,
                                               ldb,
                                               perm ))
        {
            return PASTIX_ERR_BADPARAMETER;
        }
    }
    return PASTIX_SUCCESS;
}

