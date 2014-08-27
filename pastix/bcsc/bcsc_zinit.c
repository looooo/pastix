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
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "order.h"
#include "csc.h"
#include "bcsc.h"

    /**
     * Fill in the lower triangular part of the blocked csc with values and
     * rows. The upper triangular part is done later if required through LU
     * factorization.
     */
void
bcsc_zInitLvalues( const pastix_csc_t       *csc,
                   const Order              *ord,
                   const SolverMatrix       *solvmtx,
                   const pastix_int_t       *col2cblk,
                         pastix_bcsc_t      *bcsc,
                         pastix_complex64_t *avals,
                         pastix_complex64_t *Lvalues )
{
    pastix_int_t itercblk, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = csc->dof;

    baseval = csc->colptr[0];

    for (itercol=0; itercol<csc->gN; itercol++)
    {
        bcsc_format_t *bcblk;
        pastix_int_t fcolnum, frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;
        itercblk = col2cblk[ itercol2 ];

        /* The block column is not stored locally, we skip it */
        if (itercblk == -1)
            continue;

        bcblk   = &(bcsc->cscftab[itercblk]);
        fcolnum = solvmtx->cblktab[itercblk].fcolnum;

        frow = csc->colptr[itercol]   - baseval;
        lrow = csc->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t iterrow  = csc->rows[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;
            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol - fcolnum;
                pastix_int_t rowidx = iterrow2;
                pastix_int_t pos = bcblk->coltab[ colidx ];

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++, pos++)
                {
                    bcsc->rowtab[ pos ] = rowidx;
                    Lvalues[ pos ] = avals[ ival ];
                }

                bcblk->coltab[ colidx ] += dof;
                assert( bcblk->coltab[ colidx ] <= bcblk->coltab[ colidx+1 ] );
            }
        }
    }
}
