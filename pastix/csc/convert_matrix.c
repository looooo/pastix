/**
 * @file convert_matrix.c
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include "drivers.h" 
// ajouter les prototypes dans driver.h

/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_csc_to_csr - Convert a matrix from csc format to csr format and inversly.
 *
 *******************************************************************************
 *
 * @param[out] csc
 *          At enter, contains the matrix in csr format.
 *          At exit, contains the matrix in csc format.
 *          Does not work with cscd.
 *
 *******************************************************************************/
void 
z_csc_to_csr( pastix_csc_t *csc )
{
  int j,k,col,row,nnz;
  int *row_csr;
  int *col_csr;
  void *val_csr;
  int *count;
  pastix_complex64_t val;
  pastix_complex64_t *valptr;
  
  if(csc->fmttype==PastixCSC){
    csc->fmttype=PastixCSR
  }else if(csc->fmttype==PastixCSR){
    csc->fmttype=PastixCSC
  }else if(csc->fmttype==PastixIJV){
    printf("csc_to_csr: Matrix should be in csc or csr format\n");
    return;
  }

  nnz=csc->colptr[csc->gN]-1
  if (NULL == (row_csr = calloc((csc->gN+1)*sizeof(int))))
    printf("csc_to_csr: Not enough memory (row_csr)\n");
  if (NULL == (col_csr = malloc(nnz*sizeof(int))))
    printf("csc_to_csr: Not enough memory (col_csr)\n");
  if (NULL == (val_csr = malloc(nnz*sizeof(pastix_complex64_t))))
    printf("csc_to_csr: Not enough memory (val_csr)\n");
  if (NULL == (count = calloc(csc->gN*sizeof(int))))
    printf("csc_to_csr: Not enough memory (count)\n");
  
  for (j=0;j<=nnz;j++){
    row_csr[csc->rows[j]]+=1;
  }
  row_csr[0]=0;
  for (j=1;j<=csc->gN;j++){
    row_csr[j]+=row_csr[j-1];
  }
  row_csr[0]=1;
  
  assert( row_csr[csc->gN] == nnz+1 );
  
  for (col=1;col<=csc->gN;col++){
    for (k=1;k<=csc->colptr[col]-csc->colptr[col-1];k++){
      row=csc->rows[csc->colptr[col-1]-2+k];
      valptr=csc->avals+csc->colptr[col-1]-2+k;
      val=*valptr;
      col_csr[row_csr[row-1]-1+count[row-1]]=col;
      val_csr[row_csr[row-1]-1+count[row-1]]=val;
      count[row-1]+=1;
    }
  }
  memFree_null(count);
  memFree_null(csc->colptr);
  memFree_null(csc->rows);
  memFree_null(csc->avals);
  csc->colptr=col_csr;
  csc->rows  =row_csr;
  csc->avals =val_csr;
}

/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * z_ijv_to_csr - Convert a matrix from ijv format to csc format.
 *
 *******************************************************************************
 *
 * @param[out] csc
 *          At enter, contains the matrix in ijv format.
 *          At exit, contains the matrix in csc format.
 *
 *******************************************************************************/
void 
z_ijv_to_csc( pastix_csc_t *csc )
{
  int j,k,nnz,limit;
  int *row_csc;
  int *col_csc;
  pastix_complex64_t *val_csc;
  pastix_complex64_t *valptr;
  
  if(csc->fmttype!=PastixIJV){
    printf("ijv_to_csc: Matrix should be in ijv format\n");
    return;
  }
  csc->fmttype=PastixCSC

  nnz=csc->colptr[csc->gN]-1
  if (NULL == (row_csc = calloc((csc->gN+1)*sizeof(int))))
    printf("ijv_to_csc: Not enough memory (row_csc)\n");
  if (NULL == (col_csc = malloc(nnz*sizeof(int))))
    printf("ijv_to_csc: Not enough memory (col_csc)\n");
  if (NULL == (val_csc = malloc(nnz*sizeof(pastix_complex64_t))))
    printf("ijv_to_csc: Not enough memory (val_csc)\n");
  
  for (j=0;j<=nnz;j++){
    col_csc[csc->rows[j]]+=1;
  }
  col_csc[0]=0;
  for (j=1;j<=csc->gN;j++){
    col_csc[j]+=col_csc[j-1];
  }
  col_csc[0]=1;
  
  assert( col_csc[csc->gN] == nnz+1 );
  
  for (i=0;i<=csc->gN;i++){
    j = (int)col_csc[csc->colptr[i]-1]-1;
    limit = (int)col_csc[csc->colptr[i]]-1;
    while(row_csc[j] != 0 && j < limit)
    {
      j++;
    }
    valptr=csc->avals+j;
    row_csc[j] = csc->rows[i];
    val_csc[j] = *valptr;
  }
  memFree_null(csc->colptr);
  memFree_null(csc->rows);
  memFree_null(csc->avals);
  csc->colptr=col_csc;
  csc->rows  =row_csc;
  csc->avals =val_csc;
}