/**
 * @file bcsc_fake.c
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * TODO: Use random generator from PLASMA and check that L and U are correctly
 * setup. 1. or 2. Those functions should be much simpler than this.
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
#include "d_bcsc.h"

/**
 * Fill in the lower triangular part of the blocked csc with values and
 * rows. The upper triangular part is done later if required through LU
 * factorization.
 */
static inline void
bcscInitFakeA( const pastix_csc_t  *csc,
               const Order         *ord,
               const SolverMatrix  *solvmtx,
               const pastix_int_t  *col2cblk,
                     pastix_bcsc_t *bcsc )
{
    double *Lvalues = (double*)(bcsc->Lvalues);
    pastix_int_t itercblk, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = csc->dof;

    baseval = csc->colptr[0];

    /**
     * Initialize the value of the matrix A in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<csc->gN; itercol++)
    {
        pastix_int_t *coltab;
        pastix_int_t  fcolnum, frow, lrow;
        pastix_int_t  itercol2 = ord->permtab[itercol] * dof;
        itercblk = col2cblk[ itercol2 ];

        /* The block column is not stored locally, we skip it */
        if (itercblk == -1)
            continue;

        coltab  = bcsc->cscftab[itercblk].coltab;
        fcolnum = solvmtx->cblktab[itercblk].fcolnum;

        frow = csc->colptr[itercol]   - baseval;
        lrow = csc->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t iterrow  = csc->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;
            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol - fcolnum;
                pastix_int_t rowidx = iterrow2;
                pastix_int_t pos = coltab[ colidx ];

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++, pos++)
                {
                    bcsc->rowtab[ pos ] = rowidx;
                    Lvalues[ pos ] = (rowidx >= colidx) ? 1. : 2.;
                }

                coltab[ colidx ] += dof;
                assert( coltab[ colidx ] <= coltab[ colidx+1 ] );
            }
        }
    }
}

void
bcscInitFakeLt( const pastix_csc_t  *csc,
                const Order         *ord,
                const SolverMatrix  *solvmtx,
                const pastix_int_t  *col2cblk,
                      pastix_bcsc_t *bcsc )
{
    double *Lvalues = (double*)(bcsc->Lvalues);
    pastix_int_t itercblk, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = csc->dof;

    baseval = csc->colptr[0];

    /**
     * Initialize the value of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<csc->gN; itercol++)
    {
        pastix_int_t frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;

        frow = csc->colptr[itercol]   - baseval;
        lrow = csc->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t *coltab;
            pastix_int_t fcolnum;
            pastix_int_t iterrow  = csc->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            itercblk = col2cblk[ iterrow2 ];

            /* The block column is not stored locally, we skip it */
            if ((itercblk == -1) || (iterrow == itercol))
                continue;

            coltab  = bcsc->cscftab[itercblk].coltab;
            fcolnum = solvmtx->cblktab[itercblk].fcolnum;

            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol;
                pastix_int_t rowidx = iterrow2 - fcolnum;
                pastix_int_t pos;

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++, pos++)
                {
                    pos = coltab[ rowidx ];

                    bcsc->rowtab[ pos ] = colidx;
                    Lvalues[ pos ] = (colidx >= rowidx) ? 1. : 2.;

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

void
bcscInitFakeAt( const pastix_csc_t  *csc,
                const Order         *ord,
                const SolverMatrix  *solvmtx,
                const pastix_int_t  *col2cblk,
                      pastix_int_t  *trowtab,
                      pastix_bcsc_t *bcsc )
{
    double *Uvalues = (double*)(bcsc->Uvalues);
    pastix_int_t itercblk, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = csc->dof;

    baseval = csc->colptr[0];

    /**
     * Initialize the value of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<csc->gN; itercol++)
    {
        pastix_int_t frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;

        frow = csc->colptr[itercol]   - baseval;
        lrow = csc->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t *coltab;
            pastix_int_t fcolnum;
            pastix_int_t iterrow  = csc->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            itercblk = col2cblk[ iterrow2 ];

            /* The block column is not stored locally, we skip it */
            if ((itercblk == -1) || (iterrow == itercol))
                continue;

            coltab  = bcsc->cscftab[itercblk].coltab;
            fcolnum = solvmtx->cblktab[itercblk].fcolnum;

            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol;
                pastix_int_t rowidx = iterrow2 - fcolnum;
                pastix_int_t pos;

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++, pos++)
                {
                    pos = coltab[ rowidx ];

                    trowtab[ pos ] = colidx;
                    Uvalues[ pos ] = (colidx >= rowidx) ? 2. : 1.; /*values[ ival ];*/

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

void
bcscInitCentralizedFake( const pastix_csc_t  *csc,
                         const Order         *ord,
                         const SolverMatrix  *solvmtx,
                         const pastix_int_t  *col2cblk,
                               int            initAt,
                               pastix_bcsc_t *bcsc )
{
    pastix_int_t valuesize;

    valuesize = bcsc_init_centralized_coltab( csc, ord, solvmtx, bcsc );

    /**
     * Initialize the blocked structure of the matrix A
     */
    bcscInitFakeA( csc, ord, solvmtx, col2cblk, bcsc );
    if ( csc->mtxtype == PastixSymmetric ) {
        bcscInitFakeLt( csc, ord, solvmtx, col2cblk, bcsc );
    }

    /* Restore the correct coltab arrays */
    bcsc_restore_coltab( bcsc );

    /* Sort the csc */
    d_bcscSort( bcsc, bcsc->rowtab, bcsc->Lvalues );

    if ( csc->mtxtype == PastixGeneral ) {
	/* A^t is not required if only refinment is performed */
        if (initAt) {
            pastix_int_t *trowtab;
            MALLOC_INTERN( bcsc->Uvalues, valuesize * sizeof( double ), char );
            MALLOC_INTERN( trowtab, valuesize, pastix_int_t);

            bcscInitFakeAt( csc, ord, solvmtx, col2cblk, trowtab, bcsc );

            /* Restore the correct coltab arrays */
            bcsc_restore_coltab( bcsc );

	    /* Sort the transposed csc */
	    d_bcscSort( bcsc, trowtab, bcsc->Uvalues );
	    memFree( trowtab );
        }
    }
    else {
        /* In case of PastixHermitian, conj is applied when used to save memory space */
        bcsc->Uvalues = bcsc->Lvalues;
    }
}
