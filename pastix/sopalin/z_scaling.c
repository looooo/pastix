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
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include "common.h"
#include "z_csc.h"
#include "queue.h"
#include "bulles.h"
#include "stack.h"

#include "z_sopalin_compute.h"

/*#include "symbol.h"
#include "z_ftgt.h"
#include "z_updown.h"
#include "z_solver.h" */

#define LOW 0
#define UPPER 1


void z_CscMatrix_RowMult(z_CscMatrix *csc, pastix_complex64_t *scaletab) {
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

void z_CscMatrix_ColMult(z_CscMatrix *csc, pastix_complex64_t *scaletab) {
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

void z_SolverMatrix_RowMult(z_SolverMatrix *datacode, pastix_complex64_t *scaletab, int lu) {
  int k, bloknum, i, colspan, stride, first_row;
  pastix_complex64_t *ptr;


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

void z_SolverMatrix_ColMult(z_SolverMatrix *datacode, pastix_complex64_t *scaletab, int lu) {
  int k, colspan, stride, first_col, m;
  pastix_complex64_t *ptr;


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

void z_SolverMatrix_DiagMult(z_SolverMatrix *datacode, pastix_complex64_t *scaletab) {
  int bloknum, stride, colspan, k, i;
  pastix_complex64_t *ptr, *scalptr;

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

/* Unscale the z_CscMatrix (which is not factorized) */
/* symmetric case */
void z_CscMatrix_Unscale_Sym(z_CscMatrix *cscmtx, pastix_complex64_t *scaletab, pastix_complex64_t *iscaletab) {
  (void)scaletab;
  z_CscMatrix_RowMult(cscmtx, iscaletab);
  z_CscMatrix_ColMult(cscmtx, iscaletab);
}

/* Unscale the z_CscMatrix (which is not factorized) */
/* unsymmetric case */
void z_CscMatrix_Unscale_Unsym(z_CscMatrix *cscmtx, pastix_complex64_t *scalerowtab, pastix_complex64_t *iscalerowtab, pastix_complex64_t *scalecoltab, pastix_complex64_t *iscalecoltab) {
  (void)scalerowtab; (void)scalecoltab;
  z_CscMatrix_RowMult(cscmtx, iscalerowtab);
  z_CscMatrix_ColMult(cscmtx, iscalecoltab);
}

/* Unscale the factorized z_SolverMatrix */
/* symmetric case */
void z_SolverMatrix_Unscale_Sym(z_SolverMatrix *solvmtx, pastix_complex64_t *scaletab, pastix_complex64_t *iscaletab) {
  z_SolverMatrix_RowMult(solvmtx, iscaletab, LOW);
  z_SolverMatrix_ColMult(solvmtx, scaletab, LOW);

  z_SolverMatrix_DiagMult(solvmtx, iscaletab);
  z_SolverMatrix_DiagMult(solvmtx, iscaletab);
}

/* Unscale the factorized z_SolverMatrix */
/* unsymmetric case */
void z_SolverMatrix_Unscale_Unsym(z_SolverMatrix *solvmtx, pastix_complex64_t *scalerowtab, pastix_complex64_t *iscalerowtab, pastix_complex64_t *scalecoltab, pastix_complex64_t *iscalecoltab) {
  (void)scalecoltab;
  z_SolverMatrix_RowMult(solvmtx, iscalerowtab, LOW);
  z_SolverMatrix_ColMult(solvmtx, scalerowtab,  LOW);

  z_SolverMatrix_RowMult(solvmtx, iscalecoltab, UPPER);
  z_SolverMatrix_ColMult(solvmtx, iscalerowtab, UPPER);
}

/* Unscale a matrix after it has been factorized */
/* (the z_SolverMatrix is factorized but the z_CscMatrix not) */
/* symmetric case */
void z_Matrix_Unscale_Sym(z_pastix_data_t * pastix_data,
                        z_SolverMatrix *solvmtx, pastix_complex64_t *scaletab, pastix_complex64_t *iscaletab) {
  z_SolverMatrix_Unscale_Sym(solvmtx, scaletab, iscaletab);
  z_CscMatrix_Unscale_Sym(&(pastix_data->cscmtx), scaletab, iscaletab);
}

/* Unscale a matrix after it has been factorized */
/* (the z_SolverMatrix is factorized but the z_CscMatrix not) */
/* unsymmetric case */
void z_Matrix_Unscale_Unsym(z_pastix_data_t * pastix_data,
                          z_SolverMatrix *solvmtx, pastix_complex64_t *scalerowtab, pastix_complex64_t *iscalerowtab, pastix_complex64_t *scalecoltab, pastix_complex64_t *iscalecoltab) {
  z_SolverMatrix_Unscale_Unsym(solvmtx, scalerowtab, iscalerowtab, scalecoltab, iscalecoltab);
  z_CscMatrix_Unscale_Unsym(&(pastix_data->cscmtx), scalerowtab, iscalerowtab, scalecoltab, iscalecoltab);
}
