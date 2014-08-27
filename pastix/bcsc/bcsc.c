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

/******************************************************************************
 * Function: z_CscExit                                                          *
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
bcscInitCentralized( const pastix_csc_t  *csc,
                     const Order         *ord,
                     const SolverMatrix  *solvmtx,
                           pastix_int_t   forcetr,
                           pastix_bcsc_t *bcsc )
{
    pastix_int_t  itercol, itercblk, valuesize, baseval;
    pastix_int_t  cblknbr = solvmtx->cblknbr;
    pastix_int_t *globcol  = NULL;
    pastix_int_t *col2cblk = NULL;
    pastix_int_t *colptr = csc->colptr;
    pastix_int_t *rows   = csc->rows;
    int dof = csc->dof;
    int sym = (csc->mtxtype == PastixSymmetric) || (csc->mtxtype == PastixHermitian);

    bcsc->mtxtype = csc->mtxtype;
    bcsc->flttype = csc->flttype;

    /**
     * Allocate and initialize globcol that contains the number of element in
     * each column of the symmetrized and permuted matrix
     * This corresponds to the colptr of the internal blocked csc integrating
     * symmetry and ordering.
     */
    MALLOC_INTERN( globcol, csc->gN+1, pastix_int_t );
    memset( globcol, 0, (csc->gN+1) *sizeof(pastix_int_t) );

    assert( csc->loc2glob == NULL );

    baseval = csc->colptr[0];
    {
        pastix_int_t itercol, newcol;

        for (itercol=0; itercol<csc->gN; itercol++)
        {
            pastix_int_t frow = colptr[itercol]   - baseval;
            pastix_int_t lrow = colptr[itercol+1] - baseval;
            newcol = ord->permtab[itercol];
            globcol[newcol] += lrow - frow;

            if (sym) {
                pastix_int_t iterrow, newrow;

                for (iterrow=frow; iterrow<lrow; iterrow++)
                {
                    pastix_int_t tmprow = rows[iterrow] - baseval;
                    if (tmprow != itercol) {
                        newrow = ord->permtab[tmprow];
                        globcol[newrow]++;
                    }
                }
            }
        }

        globcol[0] = 0;
        for (itercol=1; itercol<=csc->gN; itercol++)
        {
            globcol[itercol] += globcol[itercol-1];
        }
    }

    valuesize = bcsc_init_coltab( solvmtx, globcol, dof, bcsc );
    memFree_null( globcol );

    /**
     * Initialize the col2cblk array. col2cblk[i] contains the cblk index of the
     * i-th column. With distributed CSC, col2cblk[i] = -1 if not local.
     */
    {
        SolverCblk *cblk = solvmtx->cblktab;

        MALLOC_INTERN( col2cblk, (csc->gN * dof)+1, pastix_int_t );
        for (itercol=0; itercol< (csc->gN * dof)+1; itercol++)
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
    case PastixFloat:
        bcsc_sInitLvalues( csc, ord, solvmtx, col2cblk, bcsc, csc->avals, bcsc->Lvalues );
    case PastixDouble:
        bcsc_dInitLvalues( csc, ord, solvmtx, col2cblk, bcsc, csc->avals, bcsc->Lvalues );
    case PastixComplex32:
        bcsc_cInitLvalues( csc, ord, solvmtx, col2cblk, bcsc, csc->avals, bcsc->Lvalues );
    case PastixComplex64:
        bcsc_zInitLvalues( csc, ord, solvmtx, col2cblk, bcsc, csc->avals, bcsc->Lvalues );
    default:
        fprintf(stderr, "bcscInitCentralized: Error unknow floating type for input csc\n");
    }

    /* /\** */
    /*  * Fill-in the upper part of the matrix when required for LU factorization */
    /*  *\/ */
    /* if (forcetr) */
    /* { */
    /*     if (sym) */
    /*     { */
    /*         /\** */
    /*          * If PastixHermitian, conjugate is computed later if required to */
    /*          * save memory space. */
    /*          *\/ */
    /*         bcsc->Uvalues = bcsc->Lvalues; */
    /*     } */
    /*     else */
    /*     { */
    /*         MALLOC_INTERN( bcsc->Uvalues, valuesize * pastix_size_of( bcsc->flttype ), char ); */
    /*         MALLOC_INTERN( trowtab, valuesize, pastix_int_t); */
    /*         MALLOC_INTERN( trscltb, solvmtx->cblknbr, pastix_int_t *); */

    /*         for (index=0; index<solvmtx->cblknbr; index++) */
    /*         { */
    /*             MALLOC_INTERN(trscltb[index], */
    /*                           CSC_COLNBR(thecsc,index)+1, pastix_int_t); */
    /*             for (iter=0; iter<(CSC_COLNBR(thecsc,index)+1); iter++) */
    /*             { */
    /*                 trscltb[index][iter] = CSC_COL(thecsc,index,iter); */
    /*             } */
    /*         } */
    /*     } */
    /* } */
}





void bcscInit( const pastix_csc_t  *csc,
               const Order         *ord,
               const SolverMatrix  *solvmtx,
               pastix_int_t   forcetr,
               pastix_bcsc_t *bcsc )
{
    assert( ord->baseval == 0 );
    assert( ord->vertnbr == csc->n );

    if ( csc->loc2glob == NULL )
        bcscInitCentralized( csc, ord, solvmtx, forcetr, bcsc );
}
