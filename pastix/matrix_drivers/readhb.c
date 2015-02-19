/**
 * @file readhb.c
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
#include <stdio.h>
#include <libgen.h>
#include "common.h"
#include "drivers.h"
#include "iohb.h"
/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readHB - Interface to the Harwell-Boeing driver in C (iohb.c)
 * This driver can read only real matrices.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The file containing the matrix.
 *
 * @param[in] csc
 *          At exit, contains the matrix in csc format.
 *
 *******************************************************************************/

void
readHB( const char   *filename,
        pastix_csc_t *csc )
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
  char    *Type;
#define dbl dou ## ble
  dbl  *tmpval; /* hack to avoid redefinition of double... */
  int      ierr;

  Type = (char *) malloc(4*sizeof(char));

  readHB_info(filename, &Nrow2, &Ncol2, &Nnzero2, &Type, &nrhs);

  csc->n = Nrow2;
  csc->gN = csc->n;

/*   fprintf(stderr,"Matrix in file %s is %ld x %ld, with %ld nonzeros with type %s;\n", */
/*        filename, (long)*Nrow, (long)*Ncol, (long)*Nnzero, *Type); */
/*   fprintf(stderr,"%d right-hand-side(s) available.\n",nrhs); */

/*   printf("RSA: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero); */
#ifdef TYPE_COMPLEX
  fprintf(stderr,"Warning: z_HBRead is a real matrix driver, imaginary part will be 0\n");
  exit(EXIT_FAILURE);
#endif

  tmpNrow=(int)csc->n;
  tmpNcol=(int)csc->n;
  tmpNnzero=(int)Nnzero2;


  csc->colptr=(pastix_int_t*)malloc((csc->n +1)*sizeof(pastix_int_t));
  ASSERT(csc->colptr!=NULL,MOD_SI);
  tmpcol=(int*)malloc((tmpNrow+1)*sizeof(int));
  ASSERT(tmpcol!=NULL,MOD_SI);
  csc->rows=(pastix_int_t*)malloc(Nnzero2*sizeof(pastix_int_t));
  ASSERT(csc->rows!=NULL,MOD_SI);
  tmprow=(int*)malloc(tmpNnzero*sizeof(int));
  ASSERT(tmprow!=NULL,MOD_SI);
  csc->avals=(double*)malloc(Nnzero2*sizeof(double));
  ASSERT(csc->avals!=NULL,MOD_SI);

  nrhs=0;
  tmpval = (double*)malloc(tmpNnzero*sizeof(double));

  ierr = readHB_mat_double(filename, tmpcol, tmprow, tmpval);
  if(ierr == 0) {
    fprintf(stderr, "cannot read matrix (job=2)\n");
  }

  memcpy(csc->avals, tmpval, tmpNnzero*sizeof(double));
//   (*RhsType)[0]='\0';

  if (Type[0] == 'C' || Type[0] == 'c')
  {
    csc->flttype=PastixComplex64;
  }else{
    csc->flttype=PastixDouble;
  }

  if (Type[1] == 'S' || Type[1] == 's')
  {
    csc->mtxtype=PastixSymmetric;
  }else if(Type[1] == 'H' || Type[1] == 'h')
  {
    csc->mtxtype=PastixHermitian;
  }else if(Type[1] == 'U' || Type[1] == 'u')
  {
    csc->mtxtype=PastixGeneral;
  }
  for (i=0;i<tmpNrow+1;i++) (csc->colptr)[i]=(pastix_int_t)(tmpcol[i]);
  for (i=0;i<tmpNnzero;i++) (csc->rows)[i]=(pastix_int_t)(tmprow[i]);
  memFree_null(tmpcol);
  memFree_null(tmprow);
  csc->n=(pastix_int_t)tmpNrow;
  csc->gN=(pastix_int_t)tmpNcol;
  csc->fmttype = PastixCSC;
}
