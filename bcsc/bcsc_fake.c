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
#include "spm.h"
#include "solver.h"
#include "bcsc.h"
#include "d_bcsc.h"

static double alpha = 10.;
static double beta  = 1.;

/**
 * Fill in the lower triangular part of the blocked csc with values and
 * rows. The upper triangular part is done later if required through LU
 * factorization.
 */
static inline void
bcscInitFakeA( const pastix_spm_t   *spm,
               const pastix_order_t *ord,
               const SolverMatrix   *solvmtx,
               const pastix_int_t   *col2cblk,
                     pastix_bcsc_t  *bcsc )
{
    double *Lvalues = (double*)(bcsc->Lvalues);
    pastix_int_t itercblk, itercol, baseval;
    pastix_int_t i, idofcol, idofrow;
    int dof = spm->dof;
    double iival;

    baseval = spm->colptr[0];

    /**
     * Initialize the value of the matrix A in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<spm->gN; itercol++)
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

        frow = spm->colptr[itercol]   - baseval;
        lrow = spm->colptr[itercol+1] - baseval;

        if ( spm->mtxtype == PastixGeneral ) {
            iival = pastix_imax( lrow - frow - 1, 1 ) * alpha;
        }
        else {
            iival = pastix_imax( 2 * (lrow - frow - 1), 1 ) * alpha;
        }

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t iterrow  = spm->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol - fcolnum;
                pastix_int_t rowidx = iterrow2;
                pastix_int_t pos = coltab[ colidx ];

                for (idofrow = 0; idofrow < dof;
                     idofrow++, rowidx++, pos++)
                {
                    bcsc->rowtab[ pos ] = rowidx;
                    if ( rowidx == colidx ) {
                        Lvalues[ pos ] = iival;
                    }
                    else {
                        Lvalues[ pos ] = - beta;
                    }
                }

                coltab[ colidx ] += dof;
                assert( coltab[ colidx ] <= coltab[ colidx+1 ] );
            }
        }
    }
}

void
bcscInitFakeLt( const pastix_spm_t   *spm,
                const pastix_order_t *ord,
                const SolverMatrix   *solvmtx,
                const pastix_int_t   *col2cblk,
                      pastix_bcsc_t  *bcsc )
{
    double *Lvalues = (double*)(bcsc->Lvalues);
    pastix_int_t itercblk, itercol, baseval;
    pastix_int_t i, idofcol, idofrow;
    int dof = spm->dof;
    double iival;

    baseval = spm->colptr[0];

    /**
     * Initialize the value of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<spm->gN; itercol++)
    {
        pastix_int_t frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;

        frow = spm->colptr[itercol]   - baseval;
        lrow = spm->colptr[itercol+1] - baseval;

        iival = pastix_imax( 2 * (lrow - frow - 1), 1 ) * alpha;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t *coltab;
            pastix_int_t fcolnum;
            pastix_int_t iterrow  = spm->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            itercblk = col2cblk[ iterrow2 ];

            /* The block column is not stored locally, we skip it */
            if ((itercblk == -1) || (iterrow == itercol))
                continue;

            coltab  = bcsc->cscftab[itercblk].coltab;
            fcolnum = solvmtx->cblktab[itercblk].fcolnum;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol;
                pastix_int_t rowidx = iterrow2 - fcolnum;
                pastix_int_t pos;

                for (idofrow = 0; idofrow < dof;
                     idofrow++, rowidx++, pos++)
                {
                    pos = coltab[ rowidx ];

                    bcsc->rowtab[ pos ] = colidx;
                    if ( rowidx == colidx ) {
                        Lvalues[ pos ] = iival;
                    }
                    else {
                        Lvalues[ pos ] = - beta;
                    }
                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

void
bcscInitFakeAt( const pastix_spm_t   *spm,
                const pastix_order_t *ord,
                const SolverMatrix   *solvmtx,
                const pastix_int_t   *col2cblk,
                      pastix_int_t   *trowtab,
                      pastix_bcsc_t  *bcsc )
{
    double *Uvalues = (double*)(bcsc->Uvalues);
    pastix_int_t itercblk, itercol, baseval;
    pastix_int_t i, idofcol, idofrow;
    int dof = spm->dof;
    double iival;

    baseval = spm->colptr[0];

    /**
     * Initialize the value of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<spm->gN; itercol++)
    {
        pastix_int_t frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;

        frow = spm->colptr[itercol]   - baseval;
        lrow = spm->colptr[itercol+1] - baseval;

        iival = pastix_imax( lrow - frow - 1, 1 ) * alpha;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t *coltab;
            pastix_int_t fcolnum;
            pastix_int_t iterrow  = spm->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            itercblk = col2cblk[ iterrow2 ];

            /* The block column is not stored locally, we skip it */
            if ((itercblk == -1) || (iterrow == itercol))
                continue;

            coltab  = bcsc->cscftab[itercblk].coltab;
            fcolnum = solvmtx->cblktab[itercblk].fcolnum;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol;
                pastix_int_t rowidx = iterrow2 - fcolnum;
                pastix_int_t pos;

                for (idofrow = 0; idofrow < dof;
                     idofrow++, rowidx++, pos++)
                {
                    pos = coltab[ rowidx ];

                    trowtab[ pos ] = colidx;
                    if ( rowidx == colidx ) {
                        Uvalues[ pos ] = iival;
                    }
                    else {
                        Uvalues[ pos ] = - beta;
                    }

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

void
bcscInitCentralizedFake( const pastix_spm_t   *spm,
                         const pastix_order_t *ord,
                         const SolverMatrix   *solvmtx,
                         const pastix_int_t   *col2cblk,
                               int             initAt,
                               pastix_bcsc_t  *bcsc )
{
    pastix_int_t valuesize;

    bcsc->flttype = PastixDouble;
    valuesize = bcsc_init_centralized_coltab( spm, ord, solvmtx, bcsc );

    /*
     * Read environment values for alpha/beta
     */
    {
        char *str = pastix_getenv( "PASTIX_FAKE_ALPHA" );
        double value;

        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                alpha = value;
            }
            pastix_cleanenv( str );
        }

        str = pastix_getenv( "PASTIX_FAKE_BETA" );
        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                beta = value;
            }
            pastix_cleanenv( str );
        }
    }

    /*
     * Initialize the blocked structure of the matrix A
     */
    bcscInitFakeA( spm, ord, solvmtx, col2cblk, bcsc );
    if ( spm->mtxtype == PastixSymmetric ) {
        bcscInitFakeLt( spm, ord, solvmtx, col2cblk, bcsc );
    }

    /* Restore the correct coltab arrays */
    bcsc_restore_coltab( bcsc );

    /* Sort the csc */
    d_bcscSort( bcsc, bcsc->rowtab, bcsc->Lvalues );

    if ( spm->mtxtype == PastixGeneral ) {
	/* A^t is not required if only refinment is performed */
        if (initAt) {
            pastix_int_t *trowtab;
            MALLOC_INTERN( bcsc->Uvalues, valuesize * sizeof( double ), char );
            MALLOC_INTERN( trowtab, valuesize, pastix_int_t);

            bcscInitFakeAt( spm, ord, solvmtx, col2cblk, trowtab, bcsc );

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
