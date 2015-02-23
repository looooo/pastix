/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: symbol.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol.c                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the general purpose     **/
/**                routines for the symbolic matrix.       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 0.1  : from : 03 dec 1998     **/
/**                                 to     03 dec 1998     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL

#include "common.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine initializes the given
*** symbolic block matrix structure.
*** It returns:
*** - 0  : in all cases.
+*/

int
symbolInit (
SymbolMatrix * const        symbptr)
{
  memset (symbptr, 0, sizeof (SymbolMatrix));
#ifdef STARPU_GET_TASK_CTX
  symbptr->starpu_subtree_nbr=1;
#endif
  return (0);
}

/*+ This routine frees the contents
*** of the given symbolic block matrix.
*** It returns:
*** - VOID  : in all cases.
+*/

void
symbolExit (
SymbolMatrix * const        symbptr)
{
  if (symbptr->cblktab != NULL)
    memFree_null (symbptr->cblktab);
  if (symbptr->bloktab != NULL)
    memFree_null (symbptr->bloktab);

#ifdef SYMBOL_DEBUG
  symbolInit (symbptr);
#endif /* SYMBOL_DEBUG */
}

/*+ This routine reallocates the arrays
*** of the given symbolic block matrix.
*** It returns:
*** - VOID  : in all cases.
+*/

void
symbolRealloc (
SymbolMatrix * const        symbptr)
{
  SymbolCblk *        cblktab = NULL;
  SymbolBlok *        bloktab = NULL;

  if ((cblktab = (SymbolCblk *) memAlloc ((symbptr->cblknbr + 1) * sizeof (SymbolCblk))) == NULL)
    return;                                       /* Cannot move smallest array */
  memcpy  (cblktab, symbptr->cblktab, (symbptr->cblknbr + 1) * sizeof (SymbolCblk));
  memFree (symbptr->cblktab);                     /* Move column block array */
  symbptr->cblktab = cblktab;

  if ((bloktab = (SymbolBlok *) memAlloc (symbptr->bloknbr * sizeof (SymbolBlok))) == NULL)
    return;                                       /* Cannot move array */
  memcpy  (bloktab, symbptr->bloktab, symbptr->bloknbr * sizeof (SymbolBlok));
  memFree (symbptr->bloktab);                     /* Move column block array */
  symbptr->bloktab = bloktab;
}

/** Get face block for task E2 **/
pastix_int_t
symbolGetFacingBloknum(const SymbolMatrix *symbptr,
                       pastix_int_t bloksrc,
                       pastix_int_t bloknum,
                       pastix_int_t startsearch,
                       int ricar)
{
    SymbolBlok *bsrc;
    SymbolBlok *bdst;
    pastix_int_t i, fcblknum, fbloknum, lbloknum;

    fcblknum = symbptr->bloktab[bloksrc].cblknum;
    fbloknum = symbptr->cblktab[fcblknum].bloknum;
    lbloknum = symbptr->cblktab[fcblknum+1].bloknum;

    if(startsearch < fbloknum )
        startsearch = fbloknum;

    assert( startsearch < lbloknum );

    /* Block in original column block */
    bsrc = (symbptr->bloktab) + bloknum;

    /* Search for the facing block in the facing column block */
    bdst = (symbptr->bloktab) + startsearch;

    if(ricar == 0)
    {
        for(i=startsearch; i<lbloknum; i++, bdst++ )
            if( bdst->lrownum >= bsrc->frownum)
                break;

        /* We should always exit the loop in non ilu(k) mode */
        assert( (bdst->frownum <= bsrc->frownum) &&
                (bdst->lrownum >= bsrc->lrownum) );

        return i;
    }
    else
    {
        for(i=startsearch; i<lbloknum; i++, bdst++)
        {
            if( ((bsrc->frownum >= bdst->frownum) && (bsrc->frownum <= bdst->lrownum)) ||
                ((bsrc->lrownum >= bdst->frownum) && (bsrc->lrownum <= bdst->lrownum)) ||
                ((bsrc->frownum <= bdst->frownum) && (bsrc->lrownum >= bdst->lrownum)) )
                return i;  /** We found the first block that matches **/

            if(bsrc->lrownum < bdst->frownum)
            {
                return -1;
            }
        }
    }
    return -1;
}

void
symbolPrintStats( const SymbolMatrix *symbptr )
{
    SymbolCblk *cblk;
    SymbolBlok *blok;
    pastix_int_t itercblk;
    pastix_int_t cblknbr, bloknbr;
    pastix_int_t cblkmin, cblkmax;
    pastix_int_t blokmin, blokmax;
    double cblkavg, blokavg;

    cblknbr = symbptr->cblknbr;
    bloknbr = symbptr->bloknbr - cblknbr;
    cblkmin = 99999999999;
    cblkmax = 0;
    cblkavg = 0;
    blokmin = 99999999999;
    blokmax = 0;
    blokavg = 0;

    cblk = symbptr->cblktab;
    blok = symbptr->bloktab;

    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;

        pastix_int_t colnbr = cblk->lcolnum - cblk->fcolnum + 1;

        cblkmin = pastix_imin( cblkmin, colnbr );
        cblkmax = pastix_imax( cblkmax, colnbr );
        cblkavg += colnbr;
        blok++;

        /* Only extra diagonal */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            pastix_int_t rownbr = blok->lrownum - blok->frownum + 1;

            blokmin = pastix_imin( blokmin, rownbr );
            blokmax = pastix_imax( blokmax, rownbr );
            blokavg += rownbr;
        }
    }

    cblkavg = cblkavg / (double)cblknbr;
    blokavg = blokavg / (double)bloknbr;

    fprintf(stdout,
            "------ Stats Symbol Matrix ----------\n"
            " Number of cblk  : %ld\n"
            " Number of blok  : %ld\n"
            " Cblk min width  : %ld\n"
            " Cblk max width  : %ld\n"
            " Cblk avg width  : %lf\n"
            " Blok min height : %ld\n"
            " Blok max height : %ld\n"
            " Blok avg height : %lf\n"
            "-------------------------------------\n",
            cblknbr, bloknbr,
            cblkmin, cblkmax, cblkavg,
            blokmin, blokmax, blokavg );

    fprintf(stdout,
            "& %ld & %ld & %ld & %lf & %ld & %ld & %ld & %lf\n",
            cblknbr, cblkmin, cblkmax, cblkavg,
            bloknbr, blokmin, blokmax, blokavg );
}

void
symbolCheckGregoire( const SymbolMatrix *symbptr )
{
    SymbolCblk *cblk;
    SymbolBlok *blok;
    pastix_int_t itercblk, iterblok;
    pastix_int_t cblknbr, bloknbr;

    cblknbr = symbptr->cblknbr;
    bloknbr = symbptr->bloknbr - cblknbr;

    cblk = symbptr->cblktab;
    blok = symbptr->bloktab;

    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;
        pastix_int_t previous_line = -1;
        pastix_int_t previous_blok = -1;

        blok++;

        /* Only extra diagonal */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            if (blok->frownum == (previous_line + 1) && blok->cblknum == previous_blok){
                printf("Symbolic Factorization is NOT OK\n");
                exit(1);
            }
            previous_line = blok->lrownum;
            previous_blok = blok->cblknum;
        }
    }
    return;
}

void
symbolDependencies( const SymbolMatrix *symbptr )
{
    SymbolCblk *cblk;
    SymbolBlok *blok;
    pastix_int_t itercblk, iterblok;
    pastix_int_t cblknbr, bloknbr;

    cblknbr = symbptr->cblknbr;
    bloknbr = symbptr->bloknbr - cblknbr;

    cblk = symbptr->cblktab;
    blok = symbptr->bloktab;

    int i, j, l;
    pastix_int_t k;

    pastix_int_t nb_seps = symbptr->cblknbr;
    int n = symbptr->nodenbr;

    pastix_int_t **seps;
    seps = malloc( nb_seps * sizeof(pastix_int_t*));

    /* Init each separator */
    seps[0] = NULL;
    for (i=1; i < nb_seps; i++){
        SymbolCblk *cblk  = symbptr->cblktab + i;
        pastix_int_t size = cblk->lcolnum - cblk->fcolnum + 1;
        seps[i] = malloc(i * size * sizeof(pastix_int_t));
        memset(seps[i], 0, i * size * sizeof(pastix_int_t));
    }

    pastix_int_t *sort;
    sort = malloc(n * sizeof(pastix_int_t));
    memset(sort, 0, n * sizeof(pastix_int_t));

    /* Fill in lines for each extra-diagonal block */
    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;
        pastix_int_t fill;
        blok++;
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            SymbolCblk *cblk_tmp  = symbptr->cblktab + blok->cblknum;
            pastix_int_t lastdiag = cblk_tmp->fcolnum;
            for (fill = blok->frownum; fill<=blok->lrownum; fill++){
                seps[blok->cblknum][(fill-lastdiag)*blok->cblknum + itercblk] = 1;
            }
        }
    }


    /* Compute new ordering for each separator */
    for (l=1; l < nb_seps; l++){
        SymbolCblk *cblk  = symbptr->cblktab + l;
        pastix_int_t size = cblk->lcolnum - cblk->fcolnum + 1;
        pastix_int_t diag = cblk->fcolnum;

        pastix_int_t current = 0;
        for (i=0; i<size; i++){
            /* Line i is not integrated yet */
            if (sort[diag+i] == 0){
                k=0;
                while(k < l && seps[l][l*i+k] == 0){k++;}

                /* If the line is not empty */
                if (seps[l][l*i+k] != 0){
                    sort[diag+i] = current;
                    if (current == 0){
                        sort[diag+i] = -1;
                    }
                    current++;

                    for (j=i+1; j<size; j++){

                        /* Line j is not integrated yet */
                        if (sort[diag+j] == 0){
                            k = 0;
                            while (k < l && seps[l][l*i+k] == seps[l][l*j+k]){
                                k++;
                            }

                            /* If i and j have the same dependencies */
                            if (k == l){
                                sort[diag+j] = current;
                                current++;
                            }
                        }
                    }
                }
            }
        }

        /* Fill in empty lines */
        for (i=0; i<size; i++){
            if (sort[diag+i] == 0){
                sort[diag+i] = current++;
            }

            /* Special treatment for first line */
            if (sort[diag+i] == -1){
                sort[diag+i] = 0;
            }
            sort[diag+i]+= diag;
        }
    }

    /* Complete with idendity otherwise */
    for (i=0; i<n; i++){
        if (sort[i] == 0){
            sort[i] = i;
        }
    }

    FILE *file = fopen("perm.txt", "w+");
    for(i=0; i<n; i++){
        fprintf(file, "%d\n", (int)sort[i]);
    }
    fclose(file);

    /* Destroy separators */
    for (i=0; i < nb_seps; i++){
        free(seps[i]);
    }
    free(seps);
    free(sort);
    return;
}
