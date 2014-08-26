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
/*
  File: hbread.c

  Interface to the Harwell-Boeing driver in C (iohb.c)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>
#include "z_pastix.h"
#include "common_drivers.h"
#include "z_hbread.h"
#include "iohb.h"

/*
  Function: z_HBRead

  Interface to the Harwell-Boeing driver in C (iohb.c)

  Parameters
    filename - Path to the file to read from
    Nrow     - Number of rows
    Ncol     - Number of columns
    Nnzero   - Number of non zeros
    col      - Index of first element of each column in *row* and *val*
    row      - Row of eah element
    val      - Value of each element
    Type     - Type of the matrix
    RhsType  - Type of the right hand side.
 */
void z_HBRead(char const      *filename,
            pastix_int_t    *Nrow,
            pastix_int_t    *Ncol,
            pastix_int_t    *Nnzero,
            pastix_int_t   **col,
            pastix_int_t   **row,
            pastix_complex64_t **val,
            char           **Type,
            char           **RhsType)
{
  int      i;
  int      nrhs;
  int      tmpNrow;
  int      tmpNcol;
  int      tmpNnzero;
  int     *tmpcol;
  int     *tmprow;
  int      Nrow2;
  int      Ncol2;
  int      Nnzero2;
#define dbl dou ## ble
  dbl  *tmpval; /* hack to avoid redefinition of double... */
  int      ierr;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  readHB_info(filename, &Nrow2, &Ncol2, &Nnzero2, Type, &nrhs);

  *Nrow = Nrow2;
  *Ncol = Ncol2;
  *Nnzero = Nnzero2;

/*   fprintf(stderr,"Matrix in file %s is %ld x %ld, with %ld nonzeros with type %s;\n", */
/*        filename, (long)*Nrow, (long)*Ncol, (long)*Nnzero, *Type); */
/*   fprintf(stderr,"%d right-hand-side(s) available.\n",nrhs); */

/*   printf("RSA: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero); */
#ifdef TYPE_COMPLEX
  fprintf(stderr,"Warning: z_HBRead is a real matrix driver, imaginary part will be 0\n");
  exit(EXIT_FAILURE);
#endif

  tmpNrow=(int)*Nrow;
  tmpNcol=(int)*Ncol;
  tmpNnzero=(int)*Nnzero;


  *col=(pastix_int_t*)malloc((*Nrow+1)*sizeof(pastix_int_t));
  ASSERT(*col!=NULL,MOD_SI);
  tmpcol=(int*)malloc((tmpNrow+1)*sizeof(int));
  ASSERT(tmpcol!=NULL,MOD_SI);
  *row=(pastix_int_t*)malloc(*Nnzero*sizeof(pastix_int_t));
  ASSERT(*row!=NULL,MOD_SI);
  tmprow=(int*)malloc(tmpNnzero*sizeof(int));
  ASSERT(tmprow!=NULL,MOD_SI);
  *val=(pastix_complex64_t*)malloc(*Nnzero*sizeof(pastix_complex64_t));
  ASSERT(*val!=NULL,MOD_SI);

  nrhs=0;
#if (defined PREC_DOUBLE && !defined TYPE_COMPLEX)
  tmpval = *val;
#else
  tmpval = (dbl*)malloc(*Nnzero*sizeof(dbl));
#endif

  ierr = readHB_mat_double(filename, tmpcol, tmprow, tmpval);
  if(ierr == 0) {
    fprintf(stderr, "cannot read matrix (job=2)\n");
  }

#if (!defined PREC_DOUBLE || defined TYPE_COMPLEX)
  for (i = 0; i < *Nnzero; i++)
    (*val)[i] = (pastix_complex64_t)tmpval[i];
#endif
  (*RhsType)[0]='\0';
  myupcase(*Type);
  for (i=0;i<tmpNrow+1;i++) (*col)[i]=(pastix_int_t)(tmpcol[i]);
  for (i=0;i<tmpNnzero;i++) (*row)[i]=(pastix_int_t)(tmprow[i]);
  memFree_null(tmpcol);
  memFree_null(tmprow);
  *Nrow=(pastix_int_t)tmpNrow;
  *Ncol=(pastix_int_t)tmpNcol;
  *Nnzero=(pastix_int_t)tmpNnzero;
}
