/**
 * @file z_bcsc_init.c
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
#include "spm.h"
#include "solver.h"
#include "bcsc.h"
#include "z_bcsc.h"

/**
 * Fill in the lower triangular part of the blocked csc with values and
 * rows. The upper triangular part is done later if required through LU
 * factorization.
 */
static inline void
z_bcscInitA( const pastix_csc_t  *csc,
             const Order         *ord,
             const SolverMatrix  *solvmtx,
             const pastix_int_t  *col2cblk,
                   pastix_bcsc_t *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(csc->values);
    pastix_complex64_t *Lvalues = (pastix_complex64_t*)(bcsc->Lvalues);
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
                    Lvalues[ pos ] = values[ ival ];
                }

                coltab[ colidx ] += dof;
                assert( coltab[ colidx ] <= coltab[ colidx+1 ] );
            }
        }
    }
}

static inline void
z_bcscInitLt( const pastix_csc_t  *csc,
              const Order         *ord,
              const SolverMatrix  *solvmtx,
              const pastix_int_t  *col2cblk,
                    pastix_bcsc_t *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(csc->values);
    pastix_complex64_t *Lvalues = (pastix_complex64_t*)(bcsc->Lvalues);
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
                    Lvalues[ pos ] = values[ ival ];

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

#if defined(PRECISION_z) || defined(PRECISION_c)
static inline void
z_bcscInitLh( const pastix_csc_t  *csc,
              const Order         *ord,
              const SolverMatrix  *solvmtx,
              const pastix_int_t  *col2cblk,
                    pastix_bcsc_t *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(csc->values);
    pastix_complex64_t *Lvalues = (pastix_complex64_t*)(bcsc->Lvalues);
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
                    Lvalues[ pos ] = conj( values[ ival ] );

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

void
z_bcscInitAt( const pastix_csc_t  *csc,
              const Order         *ord,
              const SolverMatrix  *solvmtx,
              const pastix_int_t  *col2cblk,
                    pastix_int_t  *trowtab,
                    pastix_bcsc_t *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(csc->values);
    pastix_complex64_t *Uvalues = (pastix_complex64_t*)(bcsc->Uvalues);
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
            if (itercblk == -1)
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
                     idofrow++, ival++, rowidx++)
                {
                    pos = coltab[ rowidx ];

                    trowtab[ pos ] = colidx;
                    Uvalues[ pos ] = values[ ival ];

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

void
z_bcscSort( const pastix_bcsc_t *bcsc,
            pastix_int_t        *rowtab,
            pastix_complex64_t  *valtab )
{
    bcsc_format_t *blockcol;
    pastix_int_t itercblk, itercol, size;
    void *sortptr[2];

    blockcol = bcsc->cscftab;
    for (itercblk=0; itercblk<bcsc->cscfnbr; itercblk++, blockcol++)
    {
        for (itercol=0; itercol<blockcol->colnbr; itercol++)
        {
            int i;
            sortptr[0] = (void*)(rowtab + blockcol->coltab[itercol]);
            sortptr[1] = (void*)(valtab + blockcol->coltab[itercol]);

            size = blockcol->coltab[itercol+1] - blockcol->coltab[itercol];
            for (i=0; i<size; i++) {
                assert( rowtab[ blockcol->coltab[itercol] + i ] != -1);
            }

            z_qsortIntFloatAsc( sortptr, size );
        }
    }
}

void
z_bcscInitCentralized( const pastix_csc_t  *csc,
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
    z_bcscInitA( csc, ord, solvmtx, col2cblk, bcsc );
    if ( csc->mtxtype == PastixSymmetric ) {
        z_bcscInitLt( csc, ord, solvmtx, col2cblk, bcsc );
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    else if ( csc->mtxtype == PastixHermitian ) {
        z_bcscInitLh( csc, ord, solvmtx, col2cblk, bcsc );
    }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    /* Restore the correct coltab arrays */
    bcsc_restore_coltab( bcsc );

    /* Sort the csc */
    z_bcscSort( bcsc, bcsc->rowtab, bcsc->Lvalues );

    if ( csc->mtxtype == PastixGeneral ) {
	/* A^t is not required if only refinment is performed */
        if (initAt) {
            pastix_int_t *trowtab, i;
            MALLOC_INTERN( bcsc->Uvalues, valuesize * pastix_size_of( bcsc->flttype ), char );
            MALLOC_INTERN( trowtab, valuesize, pastix_int_t);

            for (i=0; i<valuesize; i++) {
                trowtab[i] = -1;
            }

            z_bcscInitAt( csc, ord, solvmtx, col2cblk, trowtab, bcsc );

            /* Restore the correct coltab arrays */
            bcsc_restore_coltab( bcsc );

	    /* Sort the transposed csc */
	    z_bcscSort( bcsc, trowtab, bcsc->Uvalues );
	    memFree( trowtab );
        }
    }
    else {
        /* In case of PastixHermitian, conj is applied when used to save memory space */
        bcsc->Uvalues = bcsc->Lvalues;
    }
}
