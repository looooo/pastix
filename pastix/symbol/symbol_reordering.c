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
#include "order.h"

void
symbolNewOrdering( const SymbolMatrix *symbptr, Order *order )
{
    SymbolCblk *cblk;
    SymbolBlok *blok;
    pastix_int_t itercblk, itercblk_contribute;
    Clock timer;
    int i, j;

    double time_compute_vectors  = 0.;
    double time_update_perm      = 0.;

    pastix_int_t cblknbr = symbptr->cblknbr;
    pastix_int_t LARGE_NB = 999999;

    /* First lower supernodes contribution for a given supernode */
    pastix_int_t *cblk_begin = malloc(cblknbr * sizeof(pastix_int_t));

    /* First index to consider in a column when dealing with a particular supernode */
    /* Should be replaced by an extended row-major Symbolic Matrix */
    pastix_int_t *col_begin  = calloc(cblknbr, sizeof(pastix_int_t));

    /* When looking over the symbolic structure for a given supernode, the first blok to consider in each previous supernode */
    int *nb_bloks_done = calloc(cblknbr, sizeof(int));


    /* INITIALIZATION STEP */
    {
        clockStart(timer);
        cblk_begin[0] = 0;
        for (i=1; i<cblknbr; i++){
            cblk_begin[i] = LARGE_NB;
        }

        cblk = symbptr->cblktab;
        blok = symbptr->bloktab;

        /* Find the first lower supernode to consider for each supernode */
        for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++){
            pastix_int_t iterblok = cblk[0].bloknum + 1;
            pastix_int_t lbloknum = cblk[1].bloknum;
            blok ++;

            for( ; iterblok < lbloknum; iterblok++, blok++){
                int begin = blok->frownum;
                int end   = blok->lrownum;
                pastix_int_t itercblk_c = dichotomic_search(begin, order);

                /* If the current off-diagonal block is within the considered supernode */
                if (begin >= order->rangtab[itercblk_c] &&
                    end   <  order->rangtab[itercblk_c+1]){
                    if (cblk_begin[itercblk_c] > itercblk){
                        cblk_begin[itercblk_c] = itercblk;
                    }
                }
                else{
                    printf("FATAL ERROR in dichotomic search\n");
                    exit(1);
                }
            }
        }
        clockStop(timer);
        printf("Time for initialization  %lf s\n", clockVal(timer));
    }


    /* INITIALIZING nb_bloks_done */
    {
        cblk = symbptr->cblktab;
        int nb_blok = 0;
        for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++){
            nb_bloks_done[itercblk] = nb_blok;

            pastix_int_t iterblok = cblk[0].bloknum + 1;
            pastix_int_t lbloknum = cblk[1].bloknum;
            nb_blok += lbloknum-iterblok+1;
        }
    }


    /* MAIN LOOP: for each supernode, compute the new ordering */
    /* WARNING: a first pass is to be done to reduce nb_seps for each line and then decrease the memory cost */
    for(itercblk=0; itercblk<cblknbr; itercblk++){
        clockStart(timer);
        int nb_seps  = itercblk-cblk_begin[itercblk]; /* Number of lower supernodes to consider */
        int size     = order->rangtab[itercblk+1] - order->rangtab[itercblk]; /* Supernode size */
        int *vectors = calloc(nb_seps * size, sizeof(int)); /* Hamming vectors */
        int *wmatrix = calloc(size * size, sizeof(int));    /* Hamming distances */
        int first    = cblk_begin[itercblk];

        /* Current index to fill in vectors for each supernode's line */
        int *current_index = calloc(size, sizeof(int));

        cblk = symbptr->cblktab;
        blok = symbptr->bloktab;

        /* Get the current cblk/blok numbers */
        if (first != LARGE_NB && first != 0){
            cblk += first;
            blok += nb_bloks_done[first];
        }

        /* For each related supernode */
        for(itercblk_contribute=first; itercblk_contribute<itercblk; itercblk_contribute++, cblk++){

            pastix_int_t offset   = col_begin[itercblk_contribute];
            pastix_int_t iterblok = cblk[0].bloknum + 1 + offset;
            pastix_int_t lbloknum = cblk[1].bloknum;

            /* Move to the correct offset in blok structure */
            blok += 1 + offset;

            /* Look over each off-diagonal bloks */
            for( ; iterblok < lbloknum; iterblok++, blok++){
                SymbolCblk *cblk_contribute = symbptr->cblktab + itercblk_contribute;
                int begin = blok->frownum;
                int end   = blok->lrownum;

                /* If the blok contributes */
                if (begin >= order->rangtab[itercblk] &&
                    end   <  order->rangtab[itercblk+1]){

                    for (i=blok->frownum; i<=blok->lrownum; i++){
                        int index = i - order->rangtab[itercblk];

                        vectors[index * nb_seps + current_index[index]] = itercblk_contribute;
                        current_index[index]++;
                    }
                    col_begin[itercblk_contribute]++;
                }
                /* Stop if the considered supernode is overlapped */
                else if (end < order->rangtab[itercblk+1]){
                    iterblok = lbloknum;
                }
                /* If the node belongs to a previous supernode, it should have been seen before due to col_begin use */
                else if (begin < order->rangtab[itercblk]){
                    printf("FATAL ERROR: the block should not appear here\n");
                    exit(1);
                }
            }
        }

        clockStop(timer);
        time_compute_vectors += clockVal(timer);

        /* Update the permutation for the considered supernode*/
        clockStart(timer);
        update_perm(size, wmatrix, order, itercblk, nb_seps, vectors, current_index);
        clockStop(timer);
        time_update_perm += clockVal(timer);

        free(vectors);
        free(wmatrix);
        free(current_index);
    }


    printf("Time to compute vectors  %lf s\n", time_compute_vectors);
    printf("Time to update  perm     %lf s\n", time_update_perm);

    /* Update the permutation */
    for (i=0; i<symbptr->nodenbr; i++) {
        order->permtab[ order->peritab[i] ] = i;
    }

    free(nb_bloks_done);
    free(cblk_begin);
    free(col_begin);
}

int hamming_distance_symbol( int n, int *vectors, int xi, int xj, int *current, int stop)
{
    int *end1 = vectors + n * xi + current[xi];
    int *end2 = vectors + n * xj + current[xj];
    int sum = 0;

    int *set1 = vectors + n * xi;
    int *set2 = vectors + n * xj;

    while((set1 < end1) && (set2 < end2))
    {
        if( *set1 == *set2)
        {
            set1++;
            set2++;
        }
        else if( *set1 < *set2 )
        {
            sum ++;
            set1++;
        }
        else if( *set1 > *set2 )
        {
            sum ++;
            set2++;
        }
        else
        {
            assert(0);
        }

        /* The computation is stopped if sum overlapped a given limit */
        if (sum > stop){
            return sum;
        }
    }

    sum += end1-set1;
    sum += end2-set2;
    return sum;
}


void update_perm(int sn_nvertex, int *wmatrix, Order *order, int sn_id,
                 int nb_seps, int *vectors, int *current){

    if ( sn_nvertex < 3 ) {
        return;
    }

    int  i, j, k, l, elected;
    int *tmpinvp = malloc(sn_nvertex*sizeof(int));
    int *tmplen  = malloc(sn_nvertex*sizeof(int));

  int AWK_STOP = 20;
    int stop = AWK_STOP;

    memset(tmplen, 0, sn_nvertex*sizeof(int));

    tmpinvp[0] = 0;
    tmpinvp[1] = 1;

    wmatrix[ 1 * sn_nvertex + 0 ] = hamming_distance_symbol(nb_seps, vectors, 1, 0, current, stop);

    tmplen[0] = wmatrix[ 1 * sn_nvertex + 0 ];
    tmplen[1] = wmatrix[ 1 * sn_nvertex + 0 ];

    for(i=2; i<sn_nvertex; i++) {
        /* Start by adding the row in first position */

        wmatrix[ i * sn_nvertex + tmpinvp[0] ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[0], current, stop);
        wmatrix[ i * sn_nvertex + tmpinvp[1] ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[1], current, stop);

        int minl =
            wmatrix[ i * sn_nvertex + tmpinvp[0] ] +
            wmatrix[ i * sn_nvertex + tmpinvp[1] ] - tmplen[0];
        int mpos = 1;

        int stop_loc = 1000;
        for(j=1; j<i-1; j++ ){

            if (wmatrix[ i * sn_nvertex + tmpinvp[j]] == 0)
                wmatrix[ i * sn_nvertex + tmpinvp[j]   ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[j], current, stop);

            if (wmatrix[ i * sn_nvertex + tmpinvp[j+1]] == 0)
                wmatrix[ i * sn_nvertex + tmpinvp[j+1] ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[j+1], current, stop);

            l = wmatrix[ i * sn_nvertex + tmpinvp[j]   ] +
                wmatrix[ i * sn_nvertex + tmpinvp[j+1] ] - tmplen[j];

            if ( l < minl ) {
                minl = l; mpos = j+1;

                /* WARNING: SUPPLEMENTARY TEST */
                if (l <= 0){
                    j = i;
                }
            }
        }

        /* Test between last and first */
        wmatrix[ i * sn_nvertex + tmpinvp[i-1] ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[i-1], current, stop);

        l = wmatrix[ i * sn_nvertex + tmpinvp[i-1] ] +
            wmatrix[ i * sn_nvertex + tmpinvp[0  ] ] - tmplen[i-1];
        if ( l < minl ) {
            minl = l; mpos = i;
        }

        if (mpos > 0){
            wmatrix[ i * sn_nvertex + tmpinvp[mpos-1] ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[mpos-1], current, stop);
            tmplen[mpos-1] = wmatrix[ i * sn_nvertex + tmpinvp[mpos-1] ];
        }

        if (mpos < i)
        {
            int tmpi, tmpl;
            k = i;

            wmatrix[ i * sn_nvertex + tmpinvp[mpos] ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[mpos], current, stop);
            l = wmatrix[ i * sn_nvertex + tmpinvp[mpos] ];

            /* Insert the line in the tmpinvp/tmplen arrays */
            for(j=mpos; j<i+1; j++ )
            {
                tmpi = tmpinvp[j];
                tmpl = tmplen[j];

                tmpinvp[j] = k;
                tmplen[j]  = l;

                k = tmpi;
                l = tmpl;
            }
        }
        else {
            tmpinvp[i] = i;

            wmatrix[ i * sn_nvertex + tmpinvp[0] ] = hamming_distance_symbol(nb_seps, vectors, i, tmpinvp[0], current, stop);

            tmplen[i]  = wmatrix[ i * sn_nvertex + tmpinvp[0] ];
        }
    }

    /* Search the best split line */
    {
        elected = 1;
        l = tmplen[0];
        for(i=1; i<sn_nvertex; i++) {
            if (tmplen[i] > l) {
                l = tmplen[i];
                elected = i+1;
            }
        }
    }

    pastix_int_t *sn_connected;
    MALLOC_INTERN(sn_connected, sn_nvertex, pastix_int_t);
    {
        pastix_int_t *peritab = order->peritab + order->rangtab[sn_id];
        for (i=0; i<sn_nvertex; i++)
        {
            sn_connected[i] = peritab[ tmpinvp[(i + elected)%sn_nvertex] ];
        }
        memcpy( peritab, sn_connected, sn_nvertex * sizeof(pastix_int_t) );
    }
    free(sn_connected);
    free(tmpinvp);
    free(tmplen);
}

pastix_int_t dichotomic_search( pastix_int_t node, const Order *order )
{
    pastix_int_t cblknum;
    pastix_int_t first, last, middle;

    first   = 0;
    last    = order->cblknbr;
    middle  = (last+first) / 2;
    cblknum = -1;

    while( last - first > 0 ) {
        if ( node >= order->rangtab[middle] ) {
            if ( node < order->rangtab[middle+1] ) {
                cblknum = middle;
                break;
            }
            first = middle;
        }
        else {
            last = middle;
        }
        middle = (last+first) / 2;
    }
    assert( (order->rangtab[cblknum]   <= node) &&
            (order->rangtab[cblknum+1] >  node) );
    return cblknum;
}
