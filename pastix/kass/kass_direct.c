
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
#include "SF_Direct.h"

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
        errorPrintW("Kass cannot be called for with Direct factorization and supernodes already found");
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
            nnzL = SF_Direct(mat, snodenbr, snodetab, streetab, P);
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
            amalgamate( (double)rat / 100.,
                        P, snodenbr, snodetab, streetab,
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

        Build_SymbolMatrix(P, newcblknbr, newrangtab, symbmtx);

        pastix_print(procnum, 0, "Number of block in the non patched symbol matrix = %ld \n", (long)symbmtx->bloknbr);

        memFree(invp2);
        cleanCS(P);
        memFree(P);
    }

    if (streetab != NULL ) memFree(streetab);

    dofInit(&dofstr);
    dofConstant(&dofstr, 0, symbmtx->nodenbr, 1);
    nnzS =  recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, &dofstr);
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
    Patch_SymbolMatrix(symbmtx);

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


/* void  Build_SymbolMatrix(csptr P, pastix_int_t cblknbr, pastix_int_t *rangtab, SymbolMatrix *symbmtx) */
/* { */
/*   pastix_int_t i, j, k, l; */
/*   pastix_int_t cblknum; */
/*   pastix_int_t ind; */
/*   pastix_int_t    *tmpj      = NULL; */
/*   double *tmpa      = NULL; */
/*   pastix_int_t    *node2cblk = NULL; */
/*   pastix_int_t    *ja        = NULL; */
/*   pastix_int_t n; */

/*   n = rangtab[cblknbr]; */

/*   /\**** First we transform the P matrix to find the block ****\/ */
/*   MALLOC_INTERN(tmpj, n, pastix_int_t); */
/*   MALLOC_INTERN(tmpa, n, double); */
/*   MALLOC_INTERN(node2cblk, n, pastix_int_t); */


/*   for(k=0;k<cblknbr;k++) */
/*     for(i=rangtab[k];i<rangtab[k+1];i++) */
/*       node2cblk[i] = k; */

/*   for(k=0;k<cblknbr;k++) */
/*     { */
/*       /\*i = rangtab[k];*\/ /\*OLD VERSION QUAND P pas recompacte *\/ */
/*       i = k; */
/* #ifdef DEBUG_KASS */
/*       ASSERT(P->nnzrow[i] >= (rangtab[k+1]-rangtab[k]), MOD_KASS); */

/*       for(l=0;l<rangtab[k+1]-rangtab[k];l++) */
/*         { */
/*           ASSERT(P->ja[i][l] == rangtab[k]+l, MOD_KASS); */
/*           ASSERT(node2cblk[P->ja[i][l]] == i, MOD_KASS); */
/*         } */
/* #endif */
/*       ja = P->ja[i]; */
/*       j = 0; */
/*       ind = 0; */
/*       while(j<P->nnzrow[i]) */
/*         { */
/*           cblknum = node2cblk[ja[j]]; */
/*           l=j+1; */
/*           while(l < P->nnzrow[i] && ja[l] == ja[l-1]+1 && node2cblk[ja[l]] == cblknum) */
/*             l++; */

/*           tmpj[ind] = ja[j]; */
/*           tmpa[ind] = (double)(l-j); */
/*           j = l; */
/*           ind++; */
/*         } */

/*       memFree(P->ja[i]); */
/*       P->nnzrow[i] = ind; */
/*       MALLOC_INTERN(P->ja[i], ind, pastix_int_t); */
/*       MALLOC_INTERN(P->ma[i], ind, double); */
/*       memcpy(P->ja[i], tmpj, sizeof(pastix_int_t)*ind); */
/*       memcpy(P->ma[i], tmpa, sizeof(double)*ind); */

/*     } */


/*   memFree(tmpj); */
/*   memFree(tmpa); */

/* #ifdef DEBUG_KASS */
/*   for(k=0;k<cblknbr;k++) */
/*     { */
/*       /\*i = rangtab[k];*\/ */
/*       i = k; */
/*       assert(P->nnzrow[i] > 0); */

/*       if(P->ma[i][0] != (double)(rangtab[k+1]-rangtab[k])) */
/*         print_one("Cblk %ld ma %ld rg %ld \n", k, (pastix_int_t)P->ma[i][0],rangtab[k+1]-rangtab[k]); */

/*       assert(P->ma[i][0] == (double)(rangtab[k+1]-rangtab[k])); */
/*     } */
/* #endif */

/*   /\**********************************\/ */
/*   /\*** Compute the symbol matrix ****\/ */
/*   /\**********************************\/ */
/*   symbmtx->baseval = 0; */
/*   symbmtx->cblknbr = cblknbr; */

/*   ind = 0; */
/*   symbmtx->bloknbr = CSnnz(P); */
/*   symbmtx->nodenbr = rangtab[cblknbr]; */

/*   MALLOC_INTERN(symbmtx->cblktab, cblknbr+1,        SymbolCblk); */
/*   MALLOC_INTERN(symbmtx->bloktab, symbmtx->bloknbr, SymbolBlok); */

/*   ind = 0; */
/*   for(k=0;k<cblknbr;k++) */
/*     { */
/*       symbmtx->cblktab[k].fcolnum = rangtab[k]; */
/*       symbmtx->cblktab[k].lcolnum = rangtab[k+1]-1; */
/*       symbmtx->cblktab[k].bloknum = ind; */
/*       /\*l = rangtab[k];*\/ /\** OLD VERSION **\/ */
/*       l = k; */
/*       for(i=0;i<P->nnzrow[l];i++) */
/*         { */
/*           j = P->ja[l][i]; */
/*           symbmtx->bloktab[ind].frownum = j; */
/*           symbmtx->bloktab[ind].lrownum = j+(pastix_int_t)(P->ma[l][i])-1; */
/*           symbmtx->bloktab[ind].cblknum = node2cblk[j]; */
/*           symbmtx->bloktab[ind].levfval = 0; */
/*           ind++; */
/*         } */
/* #ifdef DEBUG_KASS */
/*       assert(symbmtx->bloktab[symbmtx->cblktab[k].bloknum].frownum == symbmtx->cblktab[k].fcolnum); */
/*       assert(symbmtx->bloktab[symbmtx->cblktab[k].bloknum].lrownum == symbmtx->cblktab[k].lcolnum); */
/*       assert(symbmtx->bloktab[symbmtx->cblktab[k].bloknum].cblknum == k); */
/* #endif */


/*     } */
/*   /\*  virtual cblk to avoid side effect in the loops on cblk bloks *\/ */
/*   symbmtx->cblktab[cblknbr].fcolnum = symbmtx->cblktab[cblknbr-1].lcolnum+1; */
/*   symbmtx->cblktab[cblknbr].lcolnum = symbmtx->cblktab[cblknbr-1].lcolnum+1; */
/*   symbmtx->cblktab[cblknbr].bloknum = ind; */

/* #ifdef DEBUG_KASS */
/*   if(ind != symbmtx->bloknbr) */
/*     fprintf(stderr, "ind %ld bloknbr %ld \n", ind, symbmtx->bloknbr); */
/*   assert(ind == symbmtx->bloknbr); */
/* #endif */


/*   memFree(node2cblk); */
/* } */


/* void Patch_SymbolMatrix(SymbolMatrix *symbmtx) */
/* { */
/*   pastix_int_t i,j, k; */
/*   pastix_int_t vroot; */
/*   pastix_int_t        *father     = NULL; /\** For the cblk of the symbol matrix **\/ */
/*   SymbolBlok *newbloktab = NULL; */
/*   SymbolCblk *cblktab    = NULL; */
/*   SymbolBlok *bloktab    = NULL; */
/*   csptr Q; */


/*   cblktab = symbmtx->cblktab; */
/*   bloktab = symbmtx->bloktab; */

/*   MALLOC_INTERN(father, symbmtx->cblknbr, pastix_int_t); */
/*   MALLOC_INTERN(newbloktab, symbmtx->bloknbr + symbmtx->cblknbr, SymbolBlok); */

/*   MALLOC_INTERN(Q, 1, struct SparRow); */
/*   initCS(Q, symbmtx->cblknbr); */

/*   /\* Count how many extra-diagonal bloks are facing each diagonal blok */
/*    *\/ */
/*   for(i=0;i<symbmtx->cblknbr;i++) */
/*     for(j=cblktab[i].bloknum+1;j<cblktab[i+1].bloknum;j++) */
/*       Q->nnzrow[bloktab[j].cblknum]++; */

/*   /\* Allocate nFacingBlok integer for each diagonal blok *\/ */
/*   for(i=0;i<symbmtx->cblknbr;i++) */
/*     { */
/*       MALLOC_INTERN(Q->ja[i], Q->nnzrow[i], pastix_int_t); */
/*       Q->ma[i] = NULL; */
/*     } */

/*   for(i=0;i<symbmtx->cblknbr;i++) */
/*     Q->nnzrow[i] = 0; */

/*   /\* Q->ja[k] will contain, for each extra-diagonal facing blok */
/*    * of the column blok k, its column blok. */
/*    *\/ */
/*   for(i=0;i<symbmtx->cblknbr;i++) */
/*     for(j=cblktab[i].bloknum+1;j<cblktab[i+1].bloknum;j++) */
/*       { */
/*         k = bloktab[j].cblknum; */
/*         Q->ja[k][Q->nnzrow[k]++] = i; */
/*       } */

/*   for(i=0;i<Q->n;i++) */
/*     father[i] = -1; */

/*   for(i=0;i<Q->n;i++) */
/*     { */
/*       /\* for each blok facing diagonal blok i, */
/*        * belonging to column blok k. */
/*        * */
/*        *\/ */
/*       for(j=0;j<Q->nnzrow[i];j++) */
/*         { */
/*           k = Q->ja[i][j]; */
/* #ifdef DEBUG_KASS */
/*           assert(k<i); */
/* #endif */
/*           vroot = k; */
/*           while(father[vroot] != -1 && father[vroot] != i) */
/*             vroot = father[vroot]; */
/*           father[vroot] = i; */

/*         } */
/*     } */

/*   for(i=0;i<Q->n;i++) */
/*     if(father[i] == -1) */
/*       father[i]=i+1; */

/*   cleanCS(Q); */
/*   memFree(Q); */



/*   k = 0; */
/*   for(i=0;i<symbmtx->cblknbr-1;i++) */
/*     { */
/*       pastix_int_t odb, fbloknum; */

/*       fbloknum = cblktab[i].bloknum; */
/*       memcpy(newbloktab+k, bloktab + fbloknum, sizeof(SymbolBlok)); */
/*       cblktab[i].bloknum = k; */
/*       k++; */
/*       odb = cblktab[i+1].bloknum-fbloknum; */
/*       if(odb <= 1 || bloktab[fbloknum+1].cblknum != father[i]) */
/*         { */
/*           /\** Add a blok toward the father **\/ */
/*           newbloktab[k].frownum = cblktab[ father[i] ].fcolnum; */
/*           newbloktab[k].lrownum = cblktab[ father[i] ].fcolnum; /\** OIMBE try lcolnum **\/ */
/*           newbloktab[k].cblknum = father[i]; */
/* #ifdef DEBUG_KASS */
/*           if(father[i] != i) */
/*             assert(cblktab[father[i]].fcolnum > cblktab[i].lcolnum); */
/* #endif */

/*           newbloktab[k].levfval = 0; */
/*           k++; */
/*         } */


/*       if( odb > 1) */
/*         { */
/*           memcpy(newbloktab +k, bloktab + fbloknum+1, sizeof(SymbolBlok)*(odb-1)); */
/*           k+=odb-1; */
/*         } */

/*     } */
/*   /\** Copy the last one **\/ */
/*   memcpy(newbloktab+k, bloktab + symbmtx->cblktab[symbmtx->cblknbr-1].bloknum, sizeof(SymbolBlok)); */
/*   cblktab[symbmtx->cblknbr-1].bloknum = k; */
/*   k++; */
/*   /\** Virtual cblk **\/ */
/*   symbmtx->cblktab[symbmtx->cblknbr].bloknum = k; */

/* #ifdef DEBUG_KASS */
/*   assert(k >= symbmtx->bloknbr); */
/*   assert(k < symbmtx->cblknbr+symbmtx->bloknbr); */
/* #endif */
/*   symbmtx->bloknbr = k; */
/*   memFree(symbmtx->bloktab); */
/*   MALLOC_INTERN(symbmtx->bloktab, k, SymbolBlok); */
/*   memcpy( symbmtx->bloktab, newbloktab, sizeof(SymbolBlok)*symbmtx->bloknbr); */
/*   /\*  virtual cblk to avoid side effect in the loops on cblk bloks *\/ */
/*   cblktab[symbmtx->cblknbr].bloknum = k; */

/*   memFree(father); */
/*   memFree(newbloktab); */
/* } */
