/**
 * @file readmmd.c
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
#include "mmio.h"
/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readMMD - Read a matrix in Distributed Matrix Market format.
 * For more information about matrix market format see mmio.c/mmio.h
 * This driver can read complex and real matrices.
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
readMMD( const char   *filename,
         pastix_csc_t *csc )
{

    char    Type[4];
    FILE * file;
    int * tempcol, *g2l, *templ2g;
    int iter, iter2, baseval, mincol, maxcol;
    int * temprow;
    double * tempval_double;
    double complex * tempval_complex;
    double * valptr_double;
    double complex * valptr_complex;
    int Nnzero;
    int total;
    int tmp;
    int pos;
    int limit;
    MM_typecode matcode;
    int tmpncol,tmpnrow,tmpnnzero;
    int rank, size;
    int tmpint;
    int i;
    char line[BUFSIZ], my_filename_s[BUFSIZ];
    char * my_filename;
    int column;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    file = fopen(filename,"r");
    if (file == NULL)
    {
        fprintf(stderr,"readmmd[%d]: Cannot open the file (%s)\n", rank, filename);
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Read number of files storing the information */
    {
        int nbfiles;
        if (fscanf(line, "%d", &nbfiles) != 1) {
            fprintf(stderr,"readmmd[%d]: Cannot open the file (%s)\n", rank, filename);
            return PASTIX_ERR_IO;
        }

        if (size != nbfiles) {
            if (rank == 0) {
                fprintf(stderr,
                        "readmmd: Does'n support for incompatible number of files and processess\n"
                        "         Please rerun with %d processors\n", nbfiles);
            }
            return PASTIX_ERR_BADPARAMETER;
        }
    }

    for (i=0; i<rank+1; i++)
    {
        if ( fscanf(line, "%s", my_filename_s) != 1 ) {
            fprintf(stderr,"readmmd[%d]: Pb reading the filename for process %d\n", rank, i);
            // TODO: Clean exit while being in distributed
            return PASTIX_ERR_IO;
        }
    }
    fclose(file);

    /* Generate the full name for the local file */
    my_filename = (char*)malloc(sizeof(char)*(strlen(filename)+BUFSIZ));
    strcpy(my_filename, filename);
    sprintf(my_filename,"%s/%s", dirname(my_filename), my_filename_s);

    /* Call the centralized driver */
    rc = readMM( my_filename, csc );
    if ( rc != PASTIX_SUCCESS ) {
        // TODO: Clean exit while being in distributed
        return rc;
    }

    /* Allocation memoire */
    templ2g = (int *) malloc(Nnzero*sizeof(int));

    /* Remplissage */
    {
        long temp1,temp2;
        double re,im;

        if (mm_is_complex(matcode))
        {
            for (iter=0; iter<Nnzero; iter++)
            {
                if (4 != fscanf(file,"%ld %ld %lg %lg\n", &temp1, &temp2, &re, &im))
                {
                    fprintf(stderr, "ERROR: reading matrix (line %ld)\n",
                            (long int)iter);
                    exit(1);
                }

                iter2 = 0;
                column = temp2;
                while ( iter2 < csc->n && templ2g[iter2] < column)
                    iter2++;
                if ( !(iter2 < csc->n && (templ2g)[iter2] == column) )
                {
                    while ( iter2 < csc->n )
                    {
                        tmp = column;
                        column = (templ2g)[iter2];
                        (templ2g)[iter2] = tmp;
                        iter2++;
                    }
                    (templ2g)[iter2] = column;
                    (csc->n)++;
                }

                temprow[iter]=(int)temp1;
                tempcol[iter]=(int)temp2;
                tempval_complex[iter]=(double complex)(re+im*I);
            }
        }
        else
        {
            for (iter=0; iter<Nnzero; iter++)
            {
                if (3 != fscanf(file,"%ld %ld %lg\n", &temp1, &temp2, &re))
                {
                    fprintf(stderr, "ERROR: reading matrix (line %ld)\n",
                            (long int)iter);
                    exit(1);
                }

                iter2 = 0;
                column = temp2;
                while ( iter2 < csc->n && (templ2g)[iter2] < column)
                    iter2++;
                if ( iter2 == csc->n || (templ2g)[iter2] != column )
                {
                    while ( iter2 < csc->n )
                    {
                        tmp = (templ2g)[iter2];
                        (templ2g)[iter2] = column;
                        column = tmp;
                        iter2++;
                    }
                    (templ2g)[iter2] = column;
                    (csc->n)++;
                }

                temprow[iter]=(int)temp1;
                tempcol[iter]=(int)temp2;
                tempval_double[iter]=(double complex)(re);
            }
        }
    }

    csc->loc2glob = (int *) malloc(csc->n*sizeof(int));
    memcpy(csc->loc2glob, templ2g, csc->n*sizeof(int));
    free(templ2g);
    csc->colptr = (int *) calloc(csc->n+1,sizeof(int));
    csc->rowptr = (int *) calloc(Nnzero,sizeof(int));
    if (mm_is_complex(matcode))
    {
        csc->values = (double complex *) malloc(Nnzero*sizeof(double complex));
    }else{
        csc->values = (double *) malloc(Nnzero*sizeof(double));
    }

    if ((csc->colptr==NULL) || (csc->rowptr == NULL) || (csc->values == NULL))
    {
        fprintf(stderr, "z_MatrixMarketRead : Not enough memory (Nnzero %ld)\n",
                (long)Nnzero);
        exit(-1);
    }

    /* Detection de la base */
    mincol = csc->loc2glob[0];
    maxcol = csc->loc2glob[csc->n-1];
    if (sizeof(int) == sizeof(int))
    {
        MPI_Allreduce(&mincol, &baseval, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxcol, &csc->gN, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Allreduce(&mincol, &baseval, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxcol, &csc->gN, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    }
    if (baseval == 0)
    {
        for(iter=0; iter<Nnzero; iter++)
        {
            tempcol[iter]++;
            temprow[iter]++;
        }
        for (iter = 0; iter < csc->n; iter++)
        {
            csc->loc2glob[iter]++;
        }
    }

    if (baseval > 1 || baseval < 0)
    {
        fprintf(stderr, "Baseval > 1 || baseval < 0\n");
        exit(1);
    }
    baseval = 1;
    {
        /* Build loc2gloab */
        int inserted_column = 0;
        int iter2;
        int column;
        for (iter = 0; iter < Nnzero; iter ++)
        {
            iter2 = 0;
            column = tempcol[iter];
            while ( iter2 < inserted_column &&
                    csc->loc2glob[iter2] < column)
                iter2++;
            if (iter2 < inserted_column && csc->loc2glob[iter2] == column)
                continue;

            while ( iter2 < inserted_column )
            {
                tmp = column;
                column = csc->loc2glob[iter2];
                csc->loc2glob[iter2] = tmp;
                iter2++;
            }
            (csc->loc2glob)[iter2] = column;
            inserted_column++;
            if (inserted_column == csc->n)
                break;
        }
        assert(inserted_column == csc->n);
    }

    {
        /* Build loc2glob */
        g2l = malloc(csc->gN*sizeof(int));
        for (iter = 0; iter < csc->gN; iter++)
            g2l[iter] = -1;
        for (iter = 0; iter < (csc->n); iter ++)
            g2l[csc->loc2glob[iter]-1] = iter+1;
    }

    for (iter = 0; iter < Nnzero; iter ++)
    {
        csc->colptr[g2l[tempcol[iter]-1]-1]++;
    }

    total = baseval;
    for (iter = 0; iter < (csc->n)+1; iter ++)
    {
        tmp = csc->colptr[iter];
        csc->colptr[iter]=total;
        total+=tmp;
    }

    for (iter = 0; iter < Nnzero; iter ++)
    {

        pos = csc->colptr[g2l[tempcol[iter]-1]-1]-1;
        limit = csc->colptr[g2l[tempcol[iter]-1]]-1;
        while(csc->rowptr[pos] != 0 && pos < limit)
        {
            pos++;
        }
        if (pos == limit)
            fprintf(stderr, "Read error %ld %ld %ld\n",
                    (long int)csc->colptr[g2l[tempcol[iter]-1]-1]-1,
                    (long int)pos, (long int)limit);

        if (mm_is_complex(matcode))
        {
            valptr_complex=csc->values+pos;
            csc->rowptr[pos] = temprow[iter];
            *valptr_complex = tempval_complex[iter];
        }else{
            valptr_double=csc->values+pos;
            csc->rowptr[pos] = temprow[iter];
            *valptr_double = tempval_double[iter];
        }
    }

    if (mm_is_complex(matcode))
    {
        memFree_null(tempval_complex);
    }else{
        memFree_null(tempval_double);
    }
    memFree_null(temprow);
    memFree_null(tempcol);

    csc->fmttype = PastixCSC;
    free(my_filename);
}
