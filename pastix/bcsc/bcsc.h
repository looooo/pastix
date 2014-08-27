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
#ifndef _BCSC_H
#define _BCSC_H

#include "common.h"

struct bcsc_format_s {
    pastix_int_t   colnbr; /*> Number of columns in the block column.                                    */
    pastix_int_t * coltab; /*> Array of indexes of the start of each column in the row and value arrays. */
};

typedef struct bcsc_format_s bcsc_format_t;

/**
 * Internal column block distributed CSC matrix.
 */
struct pastix_bcsc_s {
    int            mtxtype; /*> Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.           */
    int            flttype; /*> valtab datatype: PastixFloat, PastixDouble, PastixComplex32 or PastixComplex64 */
    pastix_int_t   cscfnbr; /*> Number of fronts/column blocks.                                                */
    bcsc_format_t *cscftab; /*> Array of Block column structures of size cscfnbr. (<pastix_bcscFormat_t>)      */
    pastix_int_t  *rowtab;  /*> Array of rows in the matrix.                                                   */
    void          *Lvalues; /*> Array of values of the lower triangular part w/  diagonal of the matrix.       */
    void          *Uvalues; /*> Array of values of the upper triangular part w/o diagonal of the matrix.       */
};

typedef struct pastix_bcsc_s pastix_bcsc_t;


void bcsc_zInitLvalues( const pastix_csc_t *csc,
                        const Order        *ord,
                        const SolverMatrix *solvmtx,
                        const pastix_int_t *col2cblk,
                        pastix_bcsc_t      *bcsc,
                        pastix_complex64_t *avals,
                        pastix_complex64_t *Lvalues );

void bcsc_cInitLvalues( const pastix_csc_t *csc,
                        const Order        *ord,
                        const SolverMatrix *solvmtx,
                        const pastix_int_t *col2cblk,
                        pastix_bcsc_t      *bcsc,
                        pastix_complex32_t *avals,
                        pastix_complex32_t *Lvalues );

void bcsc_dInitLvalues( const pastix_csc_t *csc,
                        const Order        *ord,
                        const SolverMatrix *solvmtx,
                        const pastix_int_t *col2cblk,
                        pastix_bcsc_t      *bcsc,
                        double *avals,
                        double *Lvalues );

void bcsc_sInitLvalues( const pastix_csc_t *csc,
                        const Order        *ord,
                        const SolverMatrix *solvmtx,
                        const pastix_int_t *col2cblk,
                        pastix_bcsc_t      *bcsc,
                        float *avals,
                        float *Lvalues );


#endif /* CSC_H */
