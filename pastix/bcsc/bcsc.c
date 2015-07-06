/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.2
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "order.h"
#include "csc.h"
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
bcsc_init_centralized_coltab( const pastix_csc_t  *csc,
                              const Order         *ord,
                              const SolverMatrix  *solvmtx,
                                    pastix_bcsc_t *bcsc )
{
    pastix_int_t  valuesize, baseval;
    pastix_int_t *globcol  = NULL;
    pastix_int_t *colptr = csc->colptr;
    pastix_int_t *rowptr = csc->rowptr;
    int dof = csc->dof;
    int sym = (csc->mtxtype == PastixSymmetric) || (csc->mtxtype == PastixHermitian);

    bcsc->mtxtype = csc->mtxtype;
    bcsc->flttype = csc->flttype;

    baseval = csc->colptr[0];

    /**
     * Allocate and initialize globcol that contains the number of elements in
     * each column of the input matrix
     * Globcol is equivalent to the calssic colptr for the internal blocked
     * csc. The blocked csc integrate the perumtation computed within order
     * structure.
     */
    MALLOC_INTERN( globcol, csc->gN+1, pastix_int_t );
    memset( globcol, 0, (csc->gN+1) * sizeof(pastix_int_t) );

    assert( csc->loc2glob == NULL );

    {
        pastix_int_t itercol, newcol;

        for (itercol=0; itercol<csc->gN; itercol++)
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

        /* Compute displacements */
        {
            pastix_int_t tmp, idx;

            idx = 0;
            for (itercol=0; itercol<=csc->gN; itercol++)
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
bcscInitCentralized( const pastix_csc_t  *csc,
                     const Order         *ord,
                     const SolverMatrix  *solvmtx,
                           pastix_int_t   initAt,
                           pastix_bcsc_t *bcsc )
{
    pastix_int_t  itercol, itercblk;
    pastix_int_t  cblknbr = solvmtx->cblknbr;
    pastix_int_t  eltnbr  = csc->gN * csc->dof + 1;
    pastix_int_t *col2cblk = NULL;

    bcsc->mtxtype = csc->mtxtype;
    bcsc->flttype = csc->flttype;
    bcsc->gN      = csc->gN;
    bcsc->n       = csc->n;

    assert( csc->loc2glob == NULL );

    /**
     * Initialize the col2cblk array. col2cblk[i] contains the cblk index of the
     * i-th column. col2cblk[i] = -1 if not local.
     */
    {
        SolverCblk *cblk = solvmtx->cblktab;

        MALLOC_INTERN( col2cblk, eltnbr, pastix_int_t );
        for (itercol=0; itercol<eltnbr; itercol++)
            col2cblk[itercol] = -1;

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
    switch( csc->flttype ) {
    case PastixPattern:
        bcscInitCentralizedFake( csc, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixFloat:
        s_bcscInitCentralized( csc, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixDouble:
        d_bcscInitCentralized( csc, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixComplex32:
        c_bcscInitCentralized( csc, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case PastixComplex64:
        z_bcscInitCentralized( csc, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    default:
        fprintf(stderr, "bcscInitCentralized: Error unknown floating type for input csc\n");
    }

    memFree_null(col2cblk);
}

void
bcscInit( const pastix_csc_t  *csc,
          const Order         *ord,
          const SolverMatrix  *solvmtx,
                pastix_int_t   initAt,
                pastix_bcsc_t *bcsc )
{
    assert( ord->baseval == 0 );
    assert( ord->vertnbr == csc->n );

    double time = 0.;
    clockStart(time);

    if ( csc->loc2glob == NULL )
        bcscInitCentralized( csc, ord, solvmtx, initAt, bcsc );

    clockStop(time);
    fprintf(stdout, "CscdOrdistrib: %.3g s\n", clockVal(time));
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
             (bcsc->Lvalues != bcsc->Lvalues) ) {
            memFree_null( bcsc->Uvalues );
        }

        memFree_null( bcsc->Lvalues );
    }
}

int
bcscMatVec(      int            trans,
           const void          *alpha,
           const pastix_bcsc_t *bcsc,
           const void          *x,
           const void          *beta,
                 void          *y )
{
    switch (bcsc->mtxtype) {
    case PastixHermitian:
#if defined(PRECISION_z) || defined(PRECISION_c)
        switch (bcsc->flttype) {
        case PastixFloat:
            return s_bcscSymv( *((const float*)alpha), bcsc, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_bcscHemv( *((const pastix_complex32_t*)alpha), bcsc, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_bcscHemv( *((const pastix_complex64_t*)alpha), bcsc, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_bcscSymv( *((const double*)alpha), bcsc, (const double*)x, *((const double*)beta), (double*)y );
        }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */
    case PastixSymmetric:
        switch (bcsc->flttype) {
        case PastixFloat:
            return s_bcscSymv( *((const float*)alpha), bcsc, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_bcscSymv( *((const pastix_complex32_t*)alpha), bcsc, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_bcscSymv( *((const pastix_complex64_t*)alpha), bcsc, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_bcscSymv( *((const double*)alpha), bcsc, (const double*)x, *((const double*)beta), (double*)y );
        }
    case PastixGeneral:
    default:
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
}

