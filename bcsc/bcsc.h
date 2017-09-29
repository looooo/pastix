/**
 * @file bcsc.h
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _bcsc_h_
#define _bcsc_h_

struct bcsc_format_s {
    pastix_int_t   colnbr; /*> Number of columns in the block column.                                    */
    pastix_int_t * coltab; /*> Array of indexes of the start of each column in the row and value arrays. */
};

typedef struct bcsc_format_s bcsc_format_t;

/**
 * Internal column block distributed CSC matrix.
 */
struct pastix_bcsc_s {
    int            gN;      /*> Global number of vertices                                                      */
    int            n;       /*> Local number of vertices                                                       */
    int            mtxtype; /*> Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.           */
    int            flttype; /*> valtab datatype: PastixFloat, PastixDouble, PastixComplex32 or PastixComplex64 */
    pastix_int_t   cscfnbr; /*> Number of fronts/column blocks.                                                */
    bcsc_format_t *cscftab; /*> Array of Block column structures of size cscfnbr. (<pastix_bcscFormat_t>)      */
    pastix_int_t  *rowtab;  /*> Array of rows in the matrix.                                                   */
    void          *Lvalues; /*> Array of values of the matrix A                                                */
    void          *Uvalues; /*> Array of values of the matrix A^t                                              */
};

pastix_int_t
bcsc_init_centralized_coltab( const pastix_spm_t   *spm,
                              const pastix_order_t *ord,
                              const SolverMatrix   *solvmtx,
                                    pastix_bcsc_t  *bcsc );

void
bcsc_restore_coltab( pastix_bcsc_t *bcsc );

double bcscInit( const pastix_spm_t   *spm,
                 const pastix_order_t *ord,
                 const SolverMatrix   *solvmtx,
                 pastix_int_t          initAt,
                 pastix_bcsc_t        *bcsc );

void bcscExit( pastix_bcsc_t *bcsc );

int bcscMatVec(      int            trans,
               const void          *alpha,
               const pastix_bcsc_t *bcsc,
               const void          *x,
               const void          *beta,
                     void          *y );

int bcscApplyPerm( pastix_bcsc_t *bcsc,
                   pastix_int_t   n,
                   void          *b,
                   pastix_int_t   ldb,
                   pastix_int_t  *perm );

#endif /* _bcsc_h_ */
