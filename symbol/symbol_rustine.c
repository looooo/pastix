/**
 *
 * @file symbol_rustine.c
 *
 * PaStiX symbol structure functions.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 */
#include "common.h"
#include "symbol.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Patch the symbolic structure to make sure there is a tree, and not a
 * forest.
 *
 *******************************************************************************
 *
 * @param[inout] matrsymb
 *          The symbolic matrix structure with the patch.
 *
 * @param[in]    matrsymb2
 *          The imput symolic matrix that may require a patch
 *
 *******************************************************************************/
void
pastixSymbolRustine( symbol_matrix_t *       matrsymb,
                     symbol_matrix_t * const matrsymb2 )
{
    pastix_int_t i,iter,add,cblknum,bloknum,bloknum2;
    symbol_blok_t *bloktmp = NULL;
    symbol_cblk_t *cblktmp = NULL;

    MALLOC_INTERN(bloktmp, matrsymb->bloknbr+matrsymb->cblknbr, symbol_blok_t);
    MALLOC_INTERN(cblktmp, matrsymb->cblknbr+1,                 symbol_cblk_t);
    for (i=0;i<matrsymb->cblknbr+1;i++)
    {
        cblktmp[i].fcolnum = matrsymb->cblktab[i].fcolnum;
        cblktmp[i].lcolnum = matrsymb->cblktab[i].lcolnum;
        cblktmp[i].bloknum = matrsymb->cblktab[i].bloknum;
        cblktmp[i].brownum = matrsymb->cblktab[i].brownum;
    }

    iter=0,add=0;
    for (cblknum=0;cblknum<matrsymb->cblknbr-1;cblknum++)
    {
        /* recopier le bloc diagonal */
        bloknum = matrsymb->cblktab[cblknum].bloknum;
        bloktmp[iter].frownum = matrsymb->bloktab[bloknum].frownum;
        bloktmp[iter].lrownum = matrsymb->bloktab[bloknum].lrownum;
        bloktmp[iter].lcblknm = matrsymb->bloktab[bloknum].lcblknm;
        bloktmp[iter].fcblknm = matrsymb->bloktab[bloknum].fcblknm;
        iter++;

        bloknum  = matrsymb->cblktab[cblknum].bloknum+1;
        bloknum2 = matrsymb2->cblktab[cblknum].bloknum+1;

        if (bloknum == matrsymb->cblktab[cblknum+1].bloknum)
        {
            /* pas d'extra diag */
            if (matrsymb == matrsymb2)
            {
                add++;
#ifdef RUSTIN_ADD_NEXT_CBLK
                bloktmp[iter].frownum = matrsymb->cblktab[cblknum+1].fcolnum;
                bloktmp[iter].lrownum = matrsymb->cblktab[cblknum+1].fcolnum;
                bloktmp[iter].lcblknm = cblknum;
                bloktmp[iter].fcblknm = cblknum+1;
#else
                bloktmp[iter].frownum = matrsymb->cblktab[matrsymb->cblknbr-1].fcolnum;
                bloktmp[iter].lrownum = matrsymb->cblktab[matrsymb->cblknbr-1].fcolnum;
                bloktmp[iter].lcblknm = cblknum;
                bloktmp[iter].fcblknm = matrsymb->cblknbr-1;
#endif
                iter++;
            }
            else
            {
                assert( bloknum2 != matrsymb2->cblktab[cblknum+1].bloknum );
                add++;
                bloktmp[iter].frownum = matrsymb2->bloktab[bloknum2].frownum;
                bloktmp[iter].lrownum = matrsymb2->bloktab[bloknum2].frownum;
                bloktmp[iter].lcblknm = matrsymb2->bloktab[bloknum2].lcblknm;
                bloktmp[iter].fcblknm = matrsymb2->bloktab[bloknum2].fcblknm;
                iter++;
            }
        }
        else
        {
            if (matrsymb->bloktab[bloknum].fcblknm !=
                matrsymb2->bloktab[bloknum2].fcblknm)
            {
                /* le premier extra diag ne va pas */
                add++;
                bloktmp[iter].frownum = matrsymb2->bloktab[bloknum2].frownum;
                bloktmp[iter].lrownum = matrsymb2->bloktab[bloknum2].frownum;
                bloktmp[iter].lcblknm = matrsymb2->bloktab[bloknum2].lcblknm;
                bloktmp[iter].fcblknm = matrsymb2->bloktab[bloknum2].fcblknm;
                iter++;
            }

            /* on recopie tous les blocs extra du bloc-colonne */
            for (bloknum = matrsymb->cblktab[cblknum].bloknum+1;
                 bloknum < matrsymb->cblktab[cblknum+1].bloknum; bloknum++)
            {
                bloktmp[iter].frownum = matrsymb->bloktab[bloknum].frownum;
                bloktmp[iter].lrownum = matrsymb->bloktab[bloknum].lrownum;
                bloktmp[iter].lcblknm = matrsymb->bloktab[bloknum].lcblknm;
                bloktmp[iter].fcblknm = matrsymb->bloktab[bloknum].fcblknm;
                iter++;
            }
        }

        cblktmp[cblknum+1].bloknum += add;
    }

    bloktmp[iter].frownum = matrsymb->cblktab[matrsymb->cblknbr-1].fcolnum;
    bloktmp[iter].lrownum = matrsymb->cblktab[matrsymb->cblknbr-1].lcolnum;
    bloktmp[iter].lcblknm = matrsymb->cblknbr-1;
    bloktmp[iter].fcblknm = matrsymb->cblknbr-1;
    cblktmp[matrsymb->cblknbr].bloknum+=add;

    memFree_null(matrsymb->bloktab);
    memFree_null(matrsymb->cblktab);
    matrsymb->bloktab = bloktmp;
    matrsymb->cblktab = cblktmp;
    matrsymb->bloknbr += add;
    assert( add < matrsymb->cblknbr );
}
