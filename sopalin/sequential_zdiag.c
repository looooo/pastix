/**
 *
 * @file sequential_zdiag.c
 *
 *  PaStiX routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#include "cblas.h"
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

void
sequential_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data,
                  int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t  i, j, k;
    (void)pastix_data;

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        pastix_complex64_t *coeftab = cblk->lcoeftab;
        pastix_complex64_t *tmp, *lb;
        pastix_int_t size = cblk->lcolnum - cblk->fcolnum + 1;

        lb = b + cblk->lcolidx;

        if( nrhs == 1 ) {
            MALLOC_INTERN( tmp, size, pastix_complex64_t );
            cblas_zcopy( size, coeftab, cblk->stride+1, tmp, 1 );

            /* Compute */
            for (k=0; k<nrhs; k++, lb+=ldb)
            {
                for (j=0; j<size; j++) {
                    lb[j] /= tmp[j];
                }
            }
            memFree_null(tmp);
        }
        else {
            for (j=0; j<size; j++, lb++, coeftab+=(cblk->stride+1)) {
                *lb = (*lb) / (*coeftab);
            }
        }
    }
}
