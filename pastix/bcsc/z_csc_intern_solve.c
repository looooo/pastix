/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/

#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "z_csc.h"
#include "symbol.h"
#include "z_ftgt.h"
#include "z_updown.h"
#include "queue.h"
#include "bulles.h"
#include "z_solver.h"

#include "sopalin_acces.h"
#include "z_csc_intern_solve.h"

#ifdef DEBUG_RAFF
#  define CSC_LOG
#endif

/*
 * Function: z_Csc2solv_cblk
 *
 * Copy the part of the internal CSCd corresponding to
 * the column bloc itercblk into the z_SolverMatrix structure
 * coeftab which will be used to compute the decomposition.
 *
 * Used in NUMA mode.
 *
 * Parameters:
 *   cscmtx   - The internal CSCd matrix.
 *   datacode - The z_SolverMatrix structure used during decomposition.
 *   trandcsc - The internal CSCd transpose used in LU decomposition.
 *   itercblk - Column bloc number in which we had the internal CSCd.
 */
void z_Csc2solv_cblk(const z_CscMatrix *cscmtx,
                     z_SolverMatrix    *datacode,
                     pastix_complex64_t           *trandcsc,
                     pastix_int_t              itercblk)
{
  pastix_int_t itercoltab;
  pastix_int_t iterbloc;
  pastix_int_t coefindx;
  pastix_int_t iterval;

#ifdef CSC_LOG
  fprintf(stdout, "-> Csc2solv \n");
#endif

  if (itercblk < CSC_FNBR(cscmtx)){
    for (itercoltab=0;
         itercoltab < CSC_COLNBR(cscmtx,itercblk);
         itercoltab++)
      {
        for (iterval = CSC_COL(cscmtx,itercblk,itercoltab);
             iterval < CSC_COL(cscmtx,itercblk,itercoltab+1);
             iterval++)
          {
            if (CSC_ROW(cscmtx,iterval) >=
                SYMB_FCOLNUM(itercblk))
              {
                iterbloc = SYMB_BLOKNUM(itercblk);

                ASSERTDBG(iterbloc < SYMB_BLOKNBR, MOD_SOPALIN);
                while (( iterbloc < SYMB_BLOKNUM(itercblk+1)) &&
                       (( SYMB_LROWNUM(iterbloc) < CSC_ROW(cscmtx,iterval)) ||
                        ( SYMB_FROWNUM(iterbloc) > CSC_ROW(cscmtx,iterval))))
                  {
                    iterbloc++;
                  }

                if ( iterbloc < SYMB_BLOKNUM(itercblk+1) )
                  {
                    coefindx = SOLV_COEFIND(iterbloc);

                    coefindx += CSC_ROW(cscmtx,iterval) - SYMB_FROWNUM(iterbloc);

                    coefindx += SOLV_STRIDE(itercblk)*itercoltab;
                    SOLV_COEFTAB(itercblk)[coefindx] = CSC_VAL(cscmtx,iterval);

                    if (trandcsc != NULL && iterbloc != SYMB_BLOKNUM(itercblk))
                      {
                        if (cscmtx->type == 'H')
                          SOLV_UCOEFTAB(itercblk)[coefindx] = CONJ_FLOAT(trandcsc[iterval]);
                        else
                          SOLV_UCOEFTAB(itercblk)[coefindx] = trandcsc[iterval];
                      }
                  }
                else printf("ILU: csc2solv drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                            (long)datacode->cblktab[itercblk].fcolnum+
                            (long)itercoltab,(long)itercoltab,
                            (long)CSC_ROW(cscmtx,iterval),(long)iterval,
                            (long)itercblk,
                            (long)datacode->cblktab[itercblk].fcolnum,
                            (long)datacode->cblktab[itercblk].lcolnum);
              }
          }
      }
  }
#ifdef CSC_LOG
  fprintf(stdout, "<- Csc2solv \n");
#endif
}
