#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "csc.h"
#include "symbol.h"
#include "ftgt.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"

#include "sopalin_acces.h"
#include "csc_intern_solve.h"

#ifdef DEBUG_RAFF
#  define CSC_LOG
#endif

/*
 * Function: Csc2solv_cblk
 *
 * Copy the part of the internal CSCd corresponding to
 * the column bloc itercblk into the SolverMatrix structure
 * coeftab which will be used to compute the decomposition.
 *
 * Used in NUMA mode.
 *
 * Parameters:
 *   cscmtx   - The internal CSCd matrix.
 *   datacode - The SolverMatrix structure used during decomposition.
 *   trandcsc - The internal CSCd transpose used in LU decomposition.
 *   itercblk - Column bloc number in which we had the internal CSCd.
 */
void Csc2solv_cblk(const CscMatrix *cscmtx,
                   SolverMatrix    *datacode,
                   pastix_float_t           *trandcsc,
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
                            (long)datacode->symbmtx.cblktab[itercblk].fcolnum+
                            (long)itercoltab,(long)itercoltab,
                            (long)CSC_ROW(cscmtx,iterval),(long)iterval,
                            (long)itercblk,
                            (long)datacode->symbmtx.cblktab[itercblk].fcolnum,
                            (long)datacode->symbmtx.cblktab[itercblk].lcolnum);
              }
          }
      }
  }
#ifdef CSC_LOG
  fprintf(stdout, "<- Csc2solv \n");
#endif
}
