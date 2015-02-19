/**
 * @file readpetsc.c
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
#include "common.h"
#include "drivers.h"
#include "mmio.h"
#include "pastix.h"

static inline
int swap_indians_int(int num)
{
  int byte0, byte1, byte2, byte3;

  byte0 = (num & 0x000000FF) >> 0 ;
  byte1 = (num & 0x0000FF00) >> 8 ;
  byte2 = (num & 0x00FF0000) >> 16 ;
  byte3 = (num & 0xFF000000) >> 24 ;

  return((byte0 << 24) | (byte1 << 16) | (byte2 << 8) | (byte3 << 0));
}

static inline
pastix_complex64_t swap_indians_pastix_float_t(pastix_complex64_t d)
{
  size_t i;
  union
  {
    pastix_complex64_t value;
    char bytes[sizeof(pastix_complex64_t)];
  } in, out;
  in.value = d;
  for (i = 0; i < sizeof(pastix_complex64_t); i++)
    out.bytes[i] = in.bytes[sizeof(pastix_complex64_t)-1-i];

  return out.value;
}

static inline
int swap_indians_int2(int d)
{
  size_t i;
  union
  {
    int value;
    char bytes[sizeof(int)];
  } in, out;
  in.value = d;
  for (i = 0; i < sizeof(int); i++)
    out.bytes[i] = in.bytes[sizeof(int)-1-i];

  return out.value;
}

#define SWAP_INDIANS(i)                                                 \
  ( sizeof(int) == sizeof(i) )?                                         \
  ( (need_convert)?swap_indians_int(i):(i)):                            \
  ( (need_convert)?swap_indians_pastix_float_t(i):(i))

/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readPETSC - Read Matrix written by PETSc in binary format.
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file which contains the matrix stored in PETSc binary format.
 *
 * @param[out] csc
 *          At exit, contains the the matrix in csc format.
 *
 *******************************************************************************/
void
readPETSC( const char   *filename,
           pastix_csc_t *csc )
{
  char           *  buffer;
  int            *  intbuff;
  pastix_complex64_t *  floatbuff;
  unsigned long     fileLen;
  int               i, j, iter, rowsize, baseval;
  pastix_int_t      total;
  pastix_int_t      tmp;
  pastix_int_t      pos;
  pastix_int_t      limit;
	pastix_int_t      Ncol;
	pastix_int_t      Nrow;
	pastix_int_t      Nnzero;
  int               need_convert;
  int             * tempcol;
  int             * temprow;
  pastix_complex64_t  * tempval;
  pastix_complex64_t  * valptr;
  FILE            * file = fopen(filename, "rb");
  int rc;

  if (file==NULL)
  {
    fprintf(stderr,"cannot load %s\n", filename);
    exit(-1);
  }
  fseek(file, 0, SEEK_END);
  fileLen = ftell(file);
  fseek(file, 0, SEEK_SET);

  buffer = (char *)malloc(fileLen+1);
  if (!buffer)
  {
    fprintf(stderr, "Error in z_PETScRead : Not enough memory\n");
    exit(-1);
  }

  rc = fread(buffer, fileLen, 1, file);
  if (rc != (int)fileLen ) {
      perror("Error while reading file");
  }
  fclose(file);
  intbuff = (int*)buffer;
  i = (int)(*(intbuff++));
  if (i == 1211216) need_convert = 0;
  else need_convert = 1;
  i = SWAP_INDIANS(i);
  if (i != 1211216)
  {
    fprintf(stderr, "Error in z_PETScRead : Incorrect file header\n");
    exit(-1);
  }


  Nrow = (pastix_int_t)(SWAP_INDIANS(*(intbuff++)));
  Ncol = (pastix_int_t)(SWAP_INDIANS(*(intbuff++)));
  Nnzero =  (pastix_int_t)(SWAP_INDIANS(*(intbuff++)));
  fprintf(stdout, "%ld %ld %ld %ld\n", (long)i, (long)(Ncol), (long)(Nrow), (long)(Nnzero));
  if (fileLen*sizeof(char) != (4+(Nrow)+(Nnzero))*sizeof(int)+(Nnzero)*sizeof(pastix_complex64_t))
  {
    fprintf(stderr, "Error in z_PETScRead : Incorrect size of file (%ld != %ld) \n",
      (long)(fileLen*sizeof(char)),
      (long)(1*sizeof(char)+(3+(Nrow)+(Nnzero))*sizeof(int)+(Nnzero)*sizeof(pastix_complex64_t)));
    exit(-1);
  }
  if (Nnzero == -1)
  {
    fprintf(stderr, "Error in z_PETScRead : Full matrices not supported\n");
    exit(-1);
  }
  if (Nrow != Ncol)
  {
    fprintf(stderr, "Error in z_PETScRead : Non square matrices not supported\n");
    exit(-1);
  }

  tempcol = (int *) malloc((Nnzero)*sizeof(int));
  temprow = (int *) malloc((Nnzero)*sizeof(int));
  tempval = (pastix_complex64_t  *) malloc((Nnzero)*sizeof(pastix_complex64_t));

  {
    int * ttrow = temprow;

    for (i = 0; i < Nrow; i++)
    {
      rowsize = SWAP_INDIANS(*(intbuff++));
      for (j = 0; j < rowsize; j++)
        *(ttrow++) = i;
    }
  }
  for (i = 0; i < Nnzero; i++)
  {
    tempcol[i] = SWAP_INDIANS(*(intbuff++));
  }
  floatbuff = (pastix_complex64_t*)(intbuff);
  for (i = 0; i < Nnzero; i++)
  {
    tempval[i] = SWAP_INDIANS(*(floatbuff++));
  }
  free(buffer);

  csc->colptr = (pastix_int_t *) malloc((Nrow+1)*sizeof(pastix_int_t));
  memset(csc->colptr,0,(Nrow+1)*sizeof(pastix_int_t));
  csc->rows = (pastix_int_t *) malloc((Nnzero)*sizeof(pastix_int_t));
  memset(csc->rows,0,(Nnzero)*sizeof(pastix_int_t));
  csc->avals = (pastix_complex64_t *) malloc((Nnzero)*sizeof(pastix_complex64_t));
  if ((csc->colptr==NULL) || (csc->rows == NULL) || (csc->avals == NULL))
  {
    fprintf(stderr, "petscread : Not enough memory (Nnzero %ld)\n",(long)Nnzero);
    exit(-1);
  }

  baseval=1; /* Attention on base a 1 */
  for(iter=0; iter<(Nnzero); iter++)
  {
    tempcol[iter]+=baseval;
    temprow[iter]+=baseval;
  }

  for (iter = 0; iter < (Nnzero); iter ++)
  {
    csc->colptr[tempcol[iter]-1]++;
  }


  total = baseval;
  for (iter = 0; iter < (Ncol)+1; iter ++)
  {
    tmp = csc->colptr[iter];
    csc->colptr[iter]=total;
    total+=tmp;
  }

  for (iter = 0; iter < (Nnzero); iter ++)
  {

    pos = csc->colptr[tempcol[iter]-1]-1;
    limit = csc->colptr[tempcol[iter]]-1;
    while(csc->rows[pos] != 0 && pos < limit)
    {
      pos++;
    }
    if (pos == limit)
      fprintf(stderr, "Erreur de lecture\n");

    valptr=csc->avals+pos;
    csc->rows[pos] = temprow[iter];
    *valptr = tempval[iter];
  }

  memFree_null(tempval);
  memFree_null(temprow);
  memFree_null(tempcol);
  csc->n       = Ncol;
  csc->gN      = Ncol;
  csc->flttype = PastixComplex64;
  csc->mtxtype = PastixGeneral;
  csc->fmttype = PastixCSC;

  return;
}
