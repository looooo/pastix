/**
 *
 * @file bcsc_zinit.c
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @date 2021-01-03
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include <spm.h>
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "bcsc_z.h"

static inline pastix_complex64_t
__fct_id( pastix_complex64_t val ) {
    return val;
}

static inline pastix_complex64_t
__fct_conj( pastix_complex64_t val ) {
#if defined(PRECISION_c) || defined(PRECISION_z)
    return conj( val );
#else
    /* This function should not be called in this case */
    (void)val;
    assert(0);
#endif
}


/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize the values in the block csc stored in the given spm.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc values field is updated.
 *
 *******************************************************************************/
static inline void
bcsc_zinit_A( const spmatrix_t     *spm,
              const pastix_order_t *ord,
              const SolverMatrix   *solvmtx,
              const pastix_int_t   *col2cblk,
                    pastix_bcsc_t  *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Lvalues = (pastix_complex64_t*)(bcsc->Lvalues);

    pastix_int_t *colptr   = spm->colptr;
    pastix_int_t *loc2glob = spm->loc2glob;
    pastix_int_t  dof      = spm->dof;
    SolverCblk   *cblk;
    pastix_int_t *coltab;
    pastix_int_t  itercblk;
    pastix_int_t  baseval, frow, lrow;
    pastix_int_t  i, j, ig, jg, ip, jp;
    pastix_int_t  ival, idofcol, idofrow;

    baseval = spm->baseval;
    assert(spm->n == spm->gN);
    /**
     * Initialize the values of the matrix A in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for ( j=0; j<spm->n; j++, colptr++, loc2glob++ )
    {
        jg = (spm->loc2glob == NULL) ? j : *loc2glob - baseval;
        jp = ord->permtab[jg] * dof;

        itercblk = col2cblk[jp];
        if ( itercblk == -1 ) {
            continue;
        }
        cblk   = solvmtx->cblktab + itercblk;
        coltab = bcsc->cscftab[cblk->bcscnum].coltab;
        frow   = colptr[0] - baseval;
        lrow   = colptr[1] - baseval;
        for ( i=frow; i<lrow; i++ )
        {
            ig   = spm->rowptr[i] - baseval;
            ip   = ord->permtab[ig] * dof;
            ival = i * dof * dof;

            for ( idofcol = 0; idofcol < dof; idofcol++ )
            {
                pastix_int_t colidx = jp + idofcol - cblk->fcolnum;
                pastix_int_t rowidx = ip;
                pastix_int_t pos = coltab[ colidx ];

                for ( idofrow = 0; idofrow < dof;
                      idofrow++, ival++, rowidx++, pos++ )
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

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize the values in the block csc (upper part) for a matrix since
 * only one side has been initialized by bcsc_zinit_A()
 *
 * This routine will initialize either :
 *      The Symmetric upper part (L^t)
 *      The Hermitian upper part (L^h)
 *      The transpose part of A  (A^t -> U)
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[inout] rowtab
 *          The row tab of the bcsc OR
 *          The row tab associated to the transposition of A.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc values field is updated.
 *
 *******************************************************************************/
static inline void
bcsc_zinit_At( const spmatrix_t     *spm,
               const pastix_order_t *ord,
               const SolverMatrix   *solvmtx,
               const pastix_int_t   *col2cblk,
                     pastix_int_t   *rowtab,
                     pastix_bcsc_t  *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Uvalues;

    pastix_int_t *colptr   = spm->colptr;
    pastix_int_t *loc2glob = spm->loc2glob;
    pastix_int_t  dof      = spm->dof;
    SolverCblk   *cblk;
    pastix_int_t *coltab;
    pastix_int_t  i, j, ig, jg, ip, jp;
    pastix_int_t  itercblk, baseval, frow, lrow;
    pastix_int_t  ival, idofcol, idofrow;

    spm_complex64_t (*_bcsc_conj)(spm_complex64_t) = NULL;

    /* We're working on U */
    if ( spm->mtxtype == SpmGeneral ) {
        _bcsc_conj = __fct_id;
        Uvalues = (pastix_complex64_t*)(bcsc->Uvalues);
    }
    /* L^[t|h] part of the matrix */
    else {
        /*
         * precision_generator/sub.py will change SpmHermitian to SpmSymmetric
         * Don't use else or else if.
         */
        if( spm->mtxtype == SpmHermitian ){
            _bcsc_conj = __fct_conj;
        }
        if( spm->mtxtype == SpmSymmetric ){
            _bcsc_conj = __fct_id;
        }
        Uvalues = (pastix_complex64_t*)(bcsc->Lvalues);
    }

    baseval = spm->baseval;
    assert(spm->n == spm->gN);
    /**
     * Initialize the values of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for ( j=0; j<spm->n; j++, colptr++, loc2glob++ )
    {
        jg = (spm->loc2glob == NULL) ? j : *loc2glob - baseval;
        jp = ord->permtab[jg] * dof;

        frow = colptr[0] - baseval;
        lrow = colptr[1] - baseval;
        for ( i=frow; i<lrow; i++ )
        {
            ig = spm->rowptr[i]-baseval;
            ip = ord->permtab[ig] * dof;

            itercblk = col2cblk[ ip ];
            /*
             *    The block column is not stored locally
             * OR We're on a diagonal block of a symmetric matrix
             * -> we skip it
             */
            if ( ( itercblk == -1 ) ||
                 ( (ig == jg) && (spm->mtxtype != SpmGeneral) ) ) {
                continue;
            }

            cblk   = solvmtx->cblktab + itercblk;
            coltab = bcsc->cscftab[cblk->bcscnum].coltab;
            ival   = i * dof * dof;

            for ( idofcol = 0; idofcol < dof; idofcol++ )
            {
                pastix_int_t colidx = jp + idofcol;
                pastix_int_t rowidx = ip - cblk->fcolnum;
                pastix_int_t pos;

                for ( idofrow = 0; idofrow < dof;
                      idofrow++, ival++, rowidx++ )
                {
                    pos = coltab[ rowidx ];

                    rowtab[ pos ]  = colidx;
                    Uvalues[ pos ] = _bcsc_conj( values[ ival ] );

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Sort the block csc subarray associated to each column block
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *
 * @param[in] rowtab
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] valtab
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 *******************************************************************************/
static inline void
bcsc_zsort( const pastix_bcsc_t *bcsc,
            pastix_int_t        *rowtab,
            pastix_complex64_t  *valtab )
{
    bcsc_cblk_t *blockcol;
    pastix_int_t itercblk, itercol, size;
    void *sortptr[2];

    blockcol = bcsc->cscftab;
    for ( itercblk=0; itercblk<bcsc->cscfnbr; itercblk++, blockcol++ )
    {
        for ( itercol=0; itercol<blockcol->colnbr; itercol++ )
        {
            int i;
            sortptr[0] = (void*)(rowtab + blockcol->coltab[itercol]);
            sortptr[1] = (void*)(valtab + blockcol->coltab[itercol]);

            size = blockcol->coltab[itercol+1] - blockcol->coltab[itercol];
            for ( i=0; i<size; i++ ) {
                assert( rowtab[ blockcol->coltab[itercol] + i ] != -1);
            }

            z_qsortIntFloatAsc( sortptr, size );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize a centralize pastix_complex64_t block csc.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[in] initAt
 *          A flag to enable/disable the initialization of A'
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************/
void
bcsc_zinit( const spmatrix_t     *spm,
            const pastix_order_t *ord,
            const SolverMatrix   *solvmtx,
            const pastix_int_t   *col2cblk,
                  int             initAt,
                  pastix_bcsc_t  *bcsc )
{
    pastix_int_t valuesize;

    bcsc->flttype = spm->flttype;
    valuesize = bcsc_init_global_coltab( spm, ord, solvmtx, bcsc );

    /**
     * Initialize the blocked structure of the matrix A
     */
    bcsc_zinit_A( spm, ord, solvmtx, col2cblk, bcsc );
    if ( spm->mtxtype != SpmGeneral ) {
        bcsc_zinit_At( spm, ord, solvmtx, col2cblk, bcsc->rowtab, bcsc );
    }

    /* Restore the correct coltab arrays */
    bcsc_restore_coltab( bcsc );

    /* Sort the csc */
    bcsc_zsort( bcsc, bcsc->rowtab, bcsc->Lvalues );

    if ( spm->mtxtype == SpmGeneral ) {
	    /* A^t is not required if only refinement is performed */
        if (initAt) {
            pastix_int_t *trowtab, i;
            MALLOC_INTERN( bcsc->Uvalues, valuesize * pastix_size_of( bcsc->flttype ), char );
            MALLOC_INTERN( trowtab, valuesize, pastix_int_t);

            for (i=0; i<valuesize; i++) {
                trowtab[i] = -1;
            }

            bcsc_zinit_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc );

            /* Restore the correct coltab arrays */
            bcsc_restore_coltab( bcsc );

            /* Sort the transposed csc */
            bcsc_zsort( bcsc, trowtab, bcsc->Uvalues );
            memFree( trowtab );
        }
    }
    else {
        /* In case of SpmHermitian, conj is applied when used to save memory space */
        bcsc->Uvalues = bcsc->Lvalues;
    }
}
