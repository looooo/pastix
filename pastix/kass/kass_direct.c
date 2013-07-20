
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
#include "dof.h"
#include "symbol.h"
/* #include "symbol_cost.h" */
#include "order.h"
#include "fax.h"
/* #include "ifax.h" */
#include "sparRow.h"
#include "sort_row.h"
#include "SF_level.h"

#include "compact_graph.h"
#include "amalgamate.h"
#include "find_supernodes.h"
#include "kass.h"

#define print_one(fmt, ...)    if( procnum == 0) fprintf(stdout, fmt, __VA_ARGS__)

extern double nnz(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
extern double recursive_sum(pastix_int_t a, pastix_int_t b,
                            double (*fval)(pastix_int_t, const SymbolMatrix *, const Dof *),
                            const SymbolMatrix * symbmtx, const Dof * dofptr);

void Build_SymbolMatrix(csptr P, pastix_int_t cblknbr, pastix_int_t *rangtab, SymbolMatrix *symbmtx);
void Patch_SymbolMatrix(SymbolMatrix *symbmtx);

int kass2(int            ilu,
          int            levelk,
          int            rat,
          SymbolMatrix * symbmtx,
          pastix_graph_t * csc,
          Order        * orderptr,
          MPI_Comm       pastix_comm)
{
    kass_csr_t graphPA, graphL;
    pastix_int_t snodenbr;
    pastix_int_t *snodetab   = NULL;
    pastix_int_t *streetab    = NULL;
    pastix_int_t *ia         = NULL;
    pastix_int_t *ja         = NULL;
    pastix_int_t  i, j, n;
    csptr mat;
    pastix_int_t *perm;
    pastix_int_t *invp;
    pastix_int_t *invp2;
    pastix_int_t  newcblknbr;
    pastix_int_t *newrangtab = NULL;
    Dof dofstr;
    Clock timer;
    double nnzS;
    int procnum;

    MPI_Comm_rank(pastix_comm, &procnum);

#ifdef DEBUG_KASS
    print_one("--- kass begin ---\n");
#endif

    /* Check parameters correctness */
    if ( (orderptr->rangtab != NULL) && (ilu == API_NO ) ) {
        errorPrintW("Kass cannot be called for Direct factorization and with supernodes already found");
        return PASTIX_ERR_BADPARAMETER;
    }

    /*** Convert Fortran to C numbering ***/
    graphBase( csc, 0 );
    orderBase( orderptr, 0 );

    /*
     * If the supernode partition is not provided by the ordering library,
     * we compute it from scratch.
     * If we do incomplete factorization, we drop the supernode factorization
     * given by Scotch and recompute the partition.
     */
    if ( (orderptr->rangtab == NULL) || (ilu == API_YES ) )
    {
        Clock timer;

        MALLOC_INTERN(streetab, csc->n, pastix_int_t);
        clockStart(timer);
        orderFindSupernodes( csc->n, csc->colptr, csc->rows, orderptr, streetab );
        clockStop(timer);
        pastix_print(procnum, 0, "Time to find the supernode (direct) %.3g s \n", clockVal(timer));
    }

    n  = csc->n;
    ia = csc->colptr;
    ja = csc->rows;
    perm     = orderptr->permtab;
    invp     = orderptr->peritab;
    snodenbr = orderptr->cblknbr;
    snodetab = orderptr->rangtab;

    /****************************************/
    /*  Convert the graph                   */
    kass_csrInit( n, &graphPA );
    kass_csrGenPA( csc, perm, &graphPA );

    /*
     Copy the csc in a new format that split each row in a sub array, with an
     extra array that stores the nnz number.  Use the fact that the CSC is
     symmetrized, so CSC = CSR format.
     */
    {
        pastix_int_t *tmpj = NULL;
        pastix_int_t  ind;

        MALLOC_INTERN(mat, 1, struct SparRow);
        initCS(mat, n);
        MALLOC_INTERN(tmpj, n, pastix_int_t);
        /**** Convert and permute the matrix in sparrow form  ****/
        /**** The diagonal is not present in the CSR matrix, we have to put it in the matrix ***/
        bzero(tmpj, sizeof(pastix_int_t)*n);
        for(i=0;i<n;i++) {
            /*** THE GRAPH DOES NOT CONTAIN THE DIAGONAL WE ADD IT ***/
            tmpj[0] = i;
            ind = 1;
            for(j=ia[i];j<ia[i+1];j++)
                tmpj[ind++] = ja[j];

            mat->nnzrow[i] = ind;
            MALLOC_INTERN(mat->ja[i], ind, pastix_int_t);
            memcpy(mat->ja[i], tmpj, sizeof(pastix_int_t)*ind);
            mat->ma[i] = NULL;
        }
        CS_Perm(mat, perm);
        /*** Reorder the matrix ***/
        sort_row(mat);
        memFree(tmpj);
    }

    {
        pastix_int_t nnzL;
        csptr P;

        /*** Compute the ILU(k) pattern of the quotient matrix ***/
        MALLOC_INTERN(P, 1, struct SparRow);
        initCS(P, n);
        pastix_print(procnum, 0,
                     "Level of fill = %ld\n"
                     "Amalgamation ratio = %d \n",
                     (long)levelk, rat);

        /*
         * Direct Factorization
         */
        if((ilu == API_NO) || (levelk == -1))
        {
            /*
             * (Re)compute the streetab
             */
            if (streetab == NULL) { MALLOC_INTERN( streetab, n, pastix_int_t ); }

            clockStart(timer);
            nnzL = SF_Direct(&graphPA, snodenbr, snodetab, streetab, &graphL );
            clockStop(timer);
            pastix_print(procnum, 0,
                         "Time to compute scalar symbolic direct factorization  %.3g s\n",
                         clockVal(timer));
        }
        /*
         * ILU(k) Factorization
         */
        else
        {
            clockStart(timer);
            nnzL = SF_level(2, mat, levelk, P);
            clockStop(timer);
            pastix_print(procnum, 0,
                         "Time to compute scalar symbolic factorization of ILU(%ld) %.3g s\n",
                         (long)levelk, clockVal(timer));
        }
        pastix_print(procnum, 0,
                     "Scalar nnza = %ld nnzlk = %ld, fillrate0 = %.3g \n",
                     (long)( CSnnz(mat) + n)/2, (long)nnzL, (double)nnzL/(double)( (CSnnz(mat)+n)/2.0 ));


        kass_csrClean( &graphPA );
        cleanCS(mat);
        memFree(mat);

        /*
         * Sort the rows of the symbolic matrix
         */
        sort_row(P);

        MALLOC_INTERN(invp2, n, pastix_int_t);

        clockStart(timer);
        if((ilu == API_NO) || (levelk == -1))
        {
            amalgamate2( (double)rat / 100.,
                         &graphL, nnzL,
                         snodenbr, snodetab, streetab,
                         &newcblknbr, &newrangtab,
                         invp2, pastix_comm );
        }
        else
        {
            /*
             * Compute the "k-supernodes"
             */
            pastix_int_t *treetab;

            assert(streetab != NULL);

            MALLOC_INTERN(treetab, P->n, pastix_int_t);
            for(j=0;j<snodenbr;j++)
            {
                for(i=snodetab[j]; i<snodetab[j+1]-1; i++)
                    treetab[i] = i+1;

                /* Generic version */
                if( (streetab[j] == -1) ||
                    (streetab[j] == j ) )
                {
                    treetab[i] = -1;
                }
                else
                {
                    treetab[i] = snodetab[streetab[j]];
                }
            }

            /* NEW ILUK + DIRECT */
            amalgamate( (double)rat / 100.,
                        P, -1, NULL, treetab,
                        &newcblknbr, &newrangtab,
                        invp2, pastix_comm );

            memFree(treetab);
        }

        /* invp2 is the invp vector of P */
        for(i=0;i<n;i++)
            invp2[i] = invp[invp2[i]];

        memcpy(invp, invp2, sizeof(pastix_int_t)*n);
        for(i=0;i<n;i++)
            perm[invp[i]] = i;

        clockStop(timer);

        pastix_print(procnum, 0, "Time to compute the amalgamation of supernodes %.3g s\n", clockVal(timer));
        pastix_print(procnum, 0, "Number of cblk in the amalgamated symbol matrix = %ld \n", (long)newcblknbr);

        //Build_SymbolMatrix(P, newcblknbr, newrangtab, symbmtx);
        kassBuildSymbol( &graphL, newcblknbr, newrangtab, symbmtx );

        pastix_print(procnum, 0, "Number of block in the non patched symbol matrix = %ld \n", (long)symbmtx->bloknbr);

        memFree(invp2);
        cleanCS(P);
        memFree(P);
    }

    if (streetab != NULL ) memFree(streetab);

    dofInit(&dofstr);
    dofConstant(&dofstr, 0, symbmtx->nodenbr, 1);
    nnzS = recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, &dofstr);
    print_one("Number of non zero in the non patched symbol matrix = %g, fillrate1 %.3g \n",
              nnzS+n, (nnzS+n)/(ia[n]/2.0 +n));
    dofExit(&dofstr);

    if(symbolCheck(symbmtx) != 0) {
        errorPrint("SymbolCheck after kass_symbol.");
        ASSERT(0, MOD_KASS);
    }

    /********************************************************/
    /** ADD BLOCKS IN ORDER TO GET A REAL ELIMINATION TREE **/
    /********************************************************/
    //Patch_SymbolMatrix(symbmtx);
    kassPatchSymbol( symbmtx );

    dofInit(&dofstr);
    dofConstant(&dofstr, 0, symbmtx->nodenbr, 1);
    nnzS =  recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, &dofstr);
    dofExit(&dofstr);

    print_one("Number of block in final symbol matrix = %ld \n", (long)symbmtx->bloknbr);
    print_one("Number of non zero in final symbol matrix = %g, fillrate2 %.3g \n",  nnzS+n, (nnzS+n)/(ia[n]/2.0 +n));
    if( symbolCheck(symbmtx) != 0 ) {
        errorPrint("SymbolCheck after Patch_SymbolMatrix.");
        ASSERT(0, MOD_KASS);
    }

#ifdef DEBUG_KASS
    print_one("--- kass end ---\n");
#endif

    memFree(snodetab);
    orderptr->cblknbr = newcblknbr;
    orderptr->rangtab = newrangtab;
}
