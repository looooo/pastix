
/************************************************************/
/**                                                        **/
/**   NAME       : kass.c                                  **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10/02/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "kass.h"

#define print_one(fmt, ...)    if( procnum == 0) fprintf(stdout, fmt, __VA_ARGS__)

extern double nnz(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
extern double recursive_sum(pastix_int_t a, pastix_int_t b,
                            double (*fval)(pastix_int_t, const SymbolMatrix *, const Dof *),
                            const SymbolMatrix * symbmtx, const Dof * dofptr);

int
symbolKass(int             ilu,
           int             levelk,
           int             rat_cblk,
           int             rat_blas,
           SymbolMatrix   *symbmtx,
           pastix_graph_t *csc,
           Order          *orderptr,
           MPI_Comm        pastix_comm)
{
    kass_csr_t graphPA, graphL;
    pastix_int_t snodenbr;
    pastix_int_t *snodetab   = NULL;
    pastix_int_t *streetab    = NULL;
    pastix_int_t *ia         = NULL;
    pastix_int_t  i, j, n;
    pastix_int_t *perm;
    pastix_int_t *invp;
    pastix_int_t *invp2;
    pastix_int_t  newcblknbr;
    pastix_int_t *newrangtab = NULL;
    pastix_int_t  nnzA, nnzL;
    Clock timer;
    double nnzS;
    int procnum;

    MPI_Comm_rank(pastix_comm, &procnum);

    /* Check parameters correctness */
    if ( (orderptr->rangtab != NULL) && (ilu == API_NO ) )
    {
        errorPrintW("Kass cannot be called for Direct factorization and with supernodes already found");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Convert Fortran to C numbering */
    graphBase( csc, 0 );
    orderBase( orderptr, 0 );

    /*
     * If the supernode partition is not provided by the ordering library,
     * we compute it from scratch.
     * If we do incomplete factorization, we drop the supernode factorization
     * given by Scotch and recompute a new one.
     */
    if ( (orderptr->rangtab == NULL) || (ilu == API_YES ) )
    {
        Clock timer;

        MALLOC_INTERN(streetab, csc->n, pastix_int_t);
        clockStart(timer);
        orderFindSupernodes( csc->n, csc->colptr, csc->rows, orderptr, streetab );
        clockStop(timer);
        pastix_print(procnum, 0, "Time to find the supernodes %.3g s \n", clockVal(timer));
    }

    n  = csc->n;
    ia = csc->colptr;
    perm     = orderptr->permtab;
    invp     = orderptr->peritab;
    snodenbr = orderptr->cblknbr;
    snodetab = orderptr->rangtab;

    /* Create the graph of P A */
    kass_csrInit( n, &graphPA );
    kass_csrGenPA( csc, perm, &graphPA );

    pastix_print(procnum, 0,
                 "Level of fill = %ld\n"
                 "Amalgamation ratio: cblk = %d, blas = %d\n",
                 (long)levelk, rat_cblk, rat_blas);

    /*
     * Compute the graph of the factorized matrix L
     */
    /* Direct Factorization */
    if((ilu == API_NO) || (levelk == -1))
    {
        /*
         * (Re)compute the streetab
         */
        if (streetab == NULL) { MALLOC_INTERN( streetab, n, pastix_int_t ); }

        clockStart(timer);
        nnzL = kassFactDirect( &graphPA, snodenbr, snodetab, streetab, &graphL );
        clockStop(timer);
        pastix_print(procnum, 0,
                     "Time to compute scalar symbolic direct factorization  %.3g s\n",
                     clockVal(timer));
    }
    /* ILU(k) Factorization */
    else
    {
        clockStart(timer);
        nnzL = kassFactLevel( &graphPA, levelk, &graphL );
        clockStop(timer);
        pastix_print(procnum, 0,
                     "Time to compute scalar symbolic factorization of ILU(%ld) %.3g s\n",
                     (long)levelk, clockVal(timer));

        /*
         * Compute the treetab array for amalgamation step
         */
        {
            pastix_int_t *treetab;
            assert(streetab != NULL);

            MALLOC_INTERN(treetab, graphL.n, pastix_int_t);
            memcpy( treetab, streetab, graphL.n * sizeof(pastix_int_t) );
            for(j=0;j<snodenbr;j++)
            {
                for(i=snodetab[j]; i<snodetab[j+1]-1; i++)
                    streetab[i] = i+1;

                /* Generic version */
                if( (treetab[j] == -1) ||
                    (treetab[j] == j ) )
                {
                    streetab[i] = -1;
                }
                else
                {
                    streetab[i] = snodetab[treetab[j]];
                }
            }
            memFree(treetab);
        }

        snodenbr = -1;
        snodetab = NULL;
    }

    nnzA = ( kass_csrGetNNZ( &graphPA ) + n ) / 2;
    kass_csrClean( &graphPA );

    pastix_print( procnum, 0,
                  "Scalar nnza = %ld nnzlk = %ld, fillrate0 = %.3g \n",
                  (long)nnzA, (long)nnzL, (double)nnzL / (double)nnzA );

    /*
     * Amalgamate the blocks
     */
    clockStart(timer);
    MALLOC_INTERN(invp2, n, pastix_int_t);

    amalgamate( (double)rat_cblk / 100.,
                (double)rat_blas / 100.,
                &graphL, nnzL,
                snodenbr, snodetab, streetab,
                &newcblknbr, &newrangtab,
                invp2, pastix_comm );

    if( streetab != NULL ) memFree(streetab);

    if( orderptr->rangtab != NULL ) {
        memFree(orderptr->rangtab);
        orderptr->cblknbr = 0;
    }
    orderptr->cblknbr = newcblknbr;
    orderptr->rangtab = newrangtab;

    /*
     * invp2 is the invp vector generated for PA
     * we compute the invp and perm vectors of A
     */
    for(i=0;i<n;i++)
        invp2[i] = invp[invp2[i]];

    memcpy(invp, invp2, sizeof(pastix_int_t)*n);
    for(i=0;i<n;i++)
        perm[invp[i]] = i;
    memFree(invp2);
    clockStop(timer);

    /* Let's build the symbol matrix */
    kassBuildSymbol( &graphL, newcblknbr, newrangtab, symbmtx );
    kass_csrClean( &graphL );

    /********************************************************/
    /** ADD BLOCKS IN ORDER TO GET A REAL ELIMINATION TREE **/
    /********************************************************/
    if (levelk != -1) {
        nnzS = recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, NULL);

        pastix_print(procnum, 0, "Number of block in the non patched symbol matrix = %ld \n",
                     (long)symbmtx->bloknbr);
        pastix_print(procnum, 0, "Number of non zero in the non patched symbol matrix = %g, fillrate1 %.3g \n",
                     nnzS+n, (nnzS+n)/(ia[n]/2.0 +n) );

        if(symbolCheck(symbmtx) != 0) {
            errorPrint("SymbolCheck on symbol matrix before patch failed !!!");
            ASSERT(0, MOD_KASS);
        }

        /* Patch the symbol matrix */
        kassPatchSymbol( symbmtx );
    }

    nnzS = recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, NULL);

    pastix_print(procnum, 0, "Time to compute the amalgamation of supernodes %.3g s\n", clockVal(timer));
    pastix_print(procnum, 0, "Number of cblk in the amalgamated symbol matrix = %ld \n", (long)newcblknbr);
    pastix_print(procnum, 0, "Number of block in final symbol matrix = %ld \n",
                 (long)symbmtx->bloknbr);
    pastix_print(procnum, 0, "Number of non zero in final symbol matrix = %g, fillrate2 %.3g \n",
                 nnzS+n, (nnzS+n)/(ia[n]/2.0 +n));

    if( symbolCheck(symbmtx) != 0 ) {
        errorPrint("SymbolCheck on final symbol matrix failed !!!");
        ASSERT(0, MOD_KASS);
    }

    return PASTIX_SUCCESS;
}
