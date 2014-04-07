#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include "common.h"
#include "csc.h"
#include "queue.h"
#include "bulles.h"
#include "stack.h"

#include "sopalin_compute.h"

/*#include "symbol.h"
#include "ftgt.h"
#include "updown.h"
#include "solver.h" */

#define LOW 0
#define UPPER 1


void CscMatrix_RowMult(CscMatrix *csc, pastix_float_t *scaletab) {
  int k, col, i, index;
  index = 0;
  for(k=0; k<CSC_FNBR(csc); k++)
    {
      for(col=0; col<CSC_COLNBR(csc,k); col++)
        {
          for(i=0; i<CSC_COL(csc,k,col+1) - CSC_COL(csc,k,col); i++)
            {
              CSC_VAL(csc, index) *= scaletab[CSC_ROW(csc, index)];
              index++;
            }
        }
    }
}

void CscMatrix_ColMult(CscMatrix *csc, pastix_float_t *scaletab) {
  int k, col, colsize, colbegin, colnum;
  colnum = 0;
  for(k=0; k<CSC_FNBR(csc); k++)
    {
      for(col=0; col<CSC_COLNBR(csc,k); col++)
        {
          colsize = CSC_COL(csc,k,col+1) - CSC_COL(csc,k,col);
          colbegin = CSC_COL(csc,k,col);
          SOPALIN_SCAL(colsize, scaletab[colnum], &CSC_VAL(csc,colbegin), iun);
          colnum++;
        }
    }
}

void SolverMatrix_RowMult(SolverMatrix *datacode, pastix_float_t *scaletab, int lu) {
  int k, bloknum, i, colspan, stride, first_row;
  pastix_float_t *ptr;


  for(k=0; k<SYMB_CBLKNBR; k++)
    {
      colspan = CBLK_COLNBR(k);
      stride = SOLV_STRIDE(k);

      for(bloknum = SYMB_BLOKNUM(k); bloknum < SYMB_BLOKNUM(k+1); bloknum++)
        {
          ptr = ((lu==LOW)?SOLV_COEFTAB(k):SOLV_UCOEFTAB(k)) + SOLV_COEFIND(bloknum);

          first_row = SYMB_FROWNUM(bloknum);

          for(i = first_row; i <= SYMB_LROWNUM(bloknum); i++)
            {
              SOPALIN_SCAL(colspan, scaletab[i], &ptr[i-first_row], stride);
            }
        }
    }
}

void SolverMatrix_ColMult(SolverMatrix *datacode, pastix_float_t *scaletab, int lu) {
  int k, colspan, stride, first_col, m;
  pastix_float_t *ptr;


  for(k=0; k<SYMB_CBLKNBR; k++)
    {
      colspan = CBLK_COLNBR(k);
      stride  = SOLV_STRIDE(k);

      ptr = ((lu==LOW)?SOLV_COEFTAB(k):SOLV_UCOEFTAB(k)) + SOLV_COEFIND(SYMB_BLOKNUM(k));

      first_col = SYMB_FCOLNUM(k);

      for(m = 0; m < colspan; m++)
        {
          SOPALIN_SCAL(stride, scaletab[first_col+m], &ptr[m*stride], iun);
        }
    }
}

void SolverMatrix_DiagMult(SolverMatrix *datacode, pastix_float_t *scaletab) {
  int bloknum, stride, colspan, k, i;
  pastix_float_t *ptr, *scalptr;

  for(k=0; k<SYMB_CBLKNBR; k++)
    {
      bloknum = SYMB_BLOKNUM(k);
      stride  = SOLV_STRIDE(k);
      colspan = CBLK_COLNBR(k);

      scalptr = scaletab + SYMB_FCOLNUM(k);
      ptr = SOLV_COEFTAB(k) + SOLV_COEFIND(bloknum);

      for(i=0; i<colspan; i++)
        {
          ptr[i*(stride+1)] *= scalptr[i];
        }
    }
}

/* Unscale the CscMatrix (which is not factorized) */
/* symmetric case */
void CscMatrix_Unscale_Sym(CscMatrix *cscmtx, pastix_float_t *scaletab, pastix_float_t *iscaletab) {
  (void)scaletab;
  CscMatrix_RowMult(cscmtx, iscaletab);
  CscMatrix_ColMult(cscmtx, iscaletab);
}

/* Unscale the CscMatrix (which is not factorized) */
/* unsymmetric case */
void CscMatrix_Unscale_Unsym(CscMatrix *cscmtx, pastix_float_t *scalerowtab, pastix_float_t *iscalerowtab, pastix_float_t *scalecoltab, pastix_float_t *iscalecoltab) {
  (void)scalerowtab; (void)scalecoltab;
  CscMatrix_RowMult(cscmtx, iscalerowtab);
  CscMatrix_ColMult(cscmtx, iscalecoltab);
}

/* Unscale the factorized SolverMatrix */
/* symmetric case */
void SolverMatrix_Unscale_Sym(SolverMatrix *solvmtx, pastix_float_t *scaletab, pastix_float_t *iscaletab) {
  SolverMatrix_RowMult(solvmtx, iscaletab, LOW);
  SolverMatrix_ColMult(solvmtx, scaletab, LOW);

  SolverMatrix_DiagMult(solvmtx, iscaletab);
  SolverMatrix_DiagMult(solvmtx, iscaletab);
}

/* Unscale the factorized SolverMatrix */
/* unsymmetric case */
void SolverMatrix_Unscale_Unsym(SolverMatrix *solvmtx, pastix_float_t *scalerowtab, pastix_float_t *iscalerowtab, pastix_float_t *scalecoltab, pastix_float_t *iscalecoltab) {
  (void)scalecoltab;
  SolverMatrix_RowMult(solvmtx, iscalerowtab, LOW);
  SolverMatrix_ColMult(solvmtx, scalerowtab,  LOW);

  SolverMatrix_RowMult(solvmtx, iscalecoltab, UPPER);
  SolverMatrix_ColMult(solvmtx, iscalerowtab, UPPER);
}

/* Unscale a matrix after it has been factorized */
/* (the SolverMatrix is factorized but the CscMatrix not) */
/* symmetric case */
void Matrix_Unscale_Sym(pastix_data_t * pastix_data,
                        SolverMatrix *solvmtx, pastix_float_t *scaletab, pastix_float_t *iscaletab) {
  SolverMatrix_Unscale_Sym(solvmtx, scaletab, iscaletab);
  CscMatrix_Unscale_Sym(&(pastix_data->cscmtx), scaletab, iscaletab);
}

/* Unscale a matrix after it has been factorized */
/* (the SolverMatrix is factorized but the CscMatrix not) */
/* unsymmetric case */
void Matrix_Unscale_Unsym(pastix_data_t * pastix_data,
                          SolverMatrix *solvmtx, pastix_float_t *scalerowtab, pastix_float_t *iscalerowtab, pastix_float_t *scalecoltab, pastix_float_t *iscalecoltab) {
  SolverMatrix_Unscale_Unsym(solvmtx, scalerowtab, iscalerowtab, scalecoltab, iscalecoltab);
  CscMatrix_Unscale_Unsym(&(pastix_data->cscmtx), scalerowtab, iscalerowtab, scalecoltab, iscalecoltab);
}
