/**
 * @file readmm.c
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

#define MIN(x,y) (((x)<(y))?(x):(y))
/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * d_readMM - Read a real matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix Market format.
 *
 * @param[out] Values
 *          At exit, contains the values of the matrix.
 *
 *******************************************************************************/
void 
d_readMM( FILE *file,
					pastix_csc_t *csc )
{
	
	int * tempcol;
	int iter,baseval;
	int * temprow;
	double * tempval;
	double * values;
	int total;
	int tmp;
	int pos;
	int Nnzero;
	int limit;
	int tmpncol,tmpnrow;
	
	/* find out size of sparse matrix .... */
	
	if (mm_read_mtx_crd_size(file, &tmpnrow, &tmpncol, &Nnzero) !=0)
		exit(1);
	
	csc->gN = tmpncol;
	csc->n = tmpncol;
	
	/* Allocation memoire */
	tempcol = (int *) malloc((Nnzero)*sizeof(int));
	temprow = (int *) malloc((Nnzero)*sizeof(int));
	tempval = (double *) malloc((Nnzero)*sizeof(double));
	values = (double *) malloc((Nnzero)*sizeof(double));
	
	if ((tempcol==NULL) || (temprow == NULL) || (tempval == NULL))
	{
		fprintf(stderr, "d_MatrixMarketRead : Not enough memory (Nnzero %ld)\n",(long)Nnzero);
		exit(-1);
	}
	
	/* Remplissage */
	{
		long temp1,temp2;
		double re;
		
		for (iter=0; iter<(Nnzero); iter++)
		{
			if (3 != fscanf(file,"%ld %ld %lg\n", &temp1, &temp2, &re))
			{
				fprintf(stderr, "ERROR: reading matrix (line %ld)\n",
								(long int)iter);
				exit(1);
			}
			temprow[iter]=(int)temp1;
			tempcol[iter]=(int)temp2;
			tempval[iter]=(double)(re);
		}
	}
	
	csc->colptr = (int*) malloc((csc->n+1)*sizeof(int));
	memset(csc->colptr,0,(csc->n+1)*sizeof(int));
	csc->rows = (int*) malloc((Nnzero)*sizeof(int));
	memset(csc->rows,0,(Nnzero)*sizeof(int));
	csc->avals = (double*) malloc((Nnzero)*sizeof(double));
	if (((csc->colptr)==NULL) || ((csc->rows) == NULL) || ((csc->avals) == NULL))
	{
		fprintf(stderr, "z_MatrixMarketRead : Not enough memory (Nnzero %ld)\n",(long)Nnzero);
		exit(-1);
	}
	
	/* Detection de la base */
	baseval = 1;
	for(iter=0; iter<(Nnzero); iter++)
		baseval = MIN(baseval, tempcol[iter]);
	if (baseval == 0)
	{
		for(iter=0; iter<(Nnzero); iter++)
		{
			tempcol[iter]++;
			temprow[iter]++;
		}
	}
	
	for (iter = 0; iter < (Nnzero); iter ++)
	{
		csc->colptr[tempcol[iter]-1]++;
	}
	
	baseval=1; /* Attention on base a 1 */
	total = baseval;
	
	for (iter = 0; iter < (csc->gN)+1; iter ++)
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
		
		csc->rows[pos] = temprow[iter];
		values[pos] = tempval[iter];
	}
	memcpy(csc->avals, values, Nnzero*sizeof(double));
	memFree_null(tempval);
	memFree_null(values);
	memFree_null(temprow);
	memFree_null(tempcol);
}

/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * z_readMM - Read a complex matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix Market format.
 *
 * @param[out] Values
 *          At exit, contains the values of the matrix.
 *
 *******************************************************************************/
void
z_readMM( FILE *file,
					pastix_csc_t *csc )
{
	
	int * tempcol;
	int iter,baseval;
	int * temprow;
	double complex * tempval;
	double complex * values;
	int total;
	int tmp;
	int pos;
	int Nnzero;
	int limit;
	int tmpncol,tmpnrow;
	
	/* find out size of sparse matrix .... */
	
	if (mm_read_mtx_crd_size(file, &tmpnrow, &tmpncol, &Nnzero) !=0)
		exit(1);
	
	csc->gN = tmpncol;
	csc->n = tmpncol;
	
	/* Allocation memoire */
	tempcol = (int *) malloc((Nnzero)*sizeof(int));
	temprow = (int *) malloc((Nnzero)*sizeof(int));
	tempval = (double complex *) malloc((Nnzero)*sizeof(double complex));
	values = (double complex *) malloc((Nnzero)*sizeof(double complex));
	
	if ((tempcol==NULL) || (temprow == NULL) || (tempval == NULL))
	{
		fprintf(stderr, "z_MatrixMarketRead : Not enough memory (Nnzero %ld)\n",(long)Nnzero);
		exit(-1);
	}
	
	/* Remplissage */
	{
		long temp1,temp2;
		double re,im;
		im = 0.0;
		
		for (iter=0; iter<(Nnzero); iter++)
		{
			if (4 != fscanf(file,"%ld %ld %lg %lg\n", &temp1, &temp2, &re, &im))
			{
				fprintf(stderr, "ERROR: reading matrix (line %ld)\n",
								(long int)iter);
				exit(1);
			}
			
			temprow[iter]=(int)temp1;
			tempcol[iter]=(int)temp2;
			tempval[iter]=(double complex)(re+im*I);
		}
	}
	
	csc->colptr = (int*) malloc((csc->n+1)*sizeof(int));
	memset(csc->colptr,0,(csc->n +1)*sizeof(int));
	csc->rows = (int*) malloc((Nnzero)*sizeof(int));
	memset(csc->rows,0,(Nnzero)*sizeof(int));
	csc->avals = (double complex *) malloc((Nnzero)*sizeof(double complex));
	if (((csc->colptr)==NULL) || ((csc->rows) == NULL) || ((csc->avals) == NULL))
	{
		fprintf(stderr, "z_MatrixMarketRead : Not enough memory (Nnzero %ld)\n",(long)Nnzero);
		exit(-1);
	}
	
	/* Detection de la base */
	baseval = 1;
	for(iter=0; iter<(Nnzero); iter++)
		baseval = MIN(baseval, tempcol[iter]);
	if (baseval == 0)
	{
		for(iter=0; iter<(Nnzero); iter++)
		{
			tempcol[iter]++;
			temprow[iter]++;
		}
	}
	
	for (iter = 0; iter < (Nnzero); iter ++)
	{
		csc->colptr[tempcol[iter]-1]++;
	}
	
	baseval=1; /* Attention on base a 1 */
	total = baseval;
	
	for (iter = 0; iter < (csc->gN)+1; iter ++)
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
		
		csc->rows[pos] = temprow[iter];
		values[pos] = tempval[iter];
	}
	memcpy(csc->avals, values, Nnzero*sizeof(double complex));
	memFree_null(tempval);
	memFree_null(values);
	memFree_null(temprow);
	memFree_null(tempcol);
}
/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readMM - Read a matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename that contains the matrix stored in Matrix Market format.
 *
 * @param[out] Nrow
 *          At exit, contains the number of rows of the matrix.
 *
 * @param[out] Ncol
 *          At exit, contains the number of columns of the matrix.
 *
 * @param[out] Nnzero
 *          At exit, contains the number of non zero entries of the matrix.
 *
 * @param[out] Type
 *          At exit, contains the type of the matrix.
 *
 * @param[out] RhsType
 *          At exit, contains the type of the right hand side.
 *
 *******************************************************************************/
void 
readMM( const char   *filename,
				pastix_csc_t *csc )
{
	
	char    Type[4];
	char    RhsType[4];
	FILE * file;
	MM_typecode matcode;
	
	RhsType[0] = '\0';
	
	file = fopen (filename,"r");
	if (file==NULL)
	{
		fprintf(stderr,"cannot load %s\n", filename);
		exit(-1);
	}
	
	if (mm_read_banner(file, &matcode) != 0)
	{
		fprintf(stderr,"Could not process Matrix Market banner.\n");
		exit(1);
	}
	
	if (mm_is_complex(matcode))
	{
		Type[0] = 'C';
		csc->flttype=PastixComplex64;
	}else{
		Type[0] = 'R';
		csc->flttype=PastixDouble;
	}
	
	Type[1] = 'U';
	csc->mtxtype=PastixGeneral;
	if (mm_is_symmetric(matcode))
	{
		Type[1] = 'S';
		csc->mtxtype=PastixSymmetric;
	}else{
		if (mm_is_hermitian(matcode))
		{
			Type[1] = 'H';
			csc->mtxtype=PastixHermitian;
		}
	}
	Type[2] = 'A';
	Type[3] = '\0';
	
	if(mm_is_complex(matcode)){
		z_readMM(file,csc);
	}else{
		d_readMM(file,csc);
	}
	
}
