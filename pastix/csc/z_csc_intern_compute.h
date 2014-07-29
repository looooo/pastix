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
 * File: z_csc_intern_compute.h
 *
 * Functions computing operations on the CSC.
 *
 */

#ifndef Z_CSC_INTERN_COMPUTE_H
#define Z_CSC_INTERN_COMPUTE_H

/*
 * Function: z_CscNorm1
 *
 * Computes the norm 1 of the internal CSCd
 *
 * Norm 1 is equal to the maximum value of the sum of the
 * Absolute values of the elements of each columns.
 *
 * Parameters:
 *   cscmtx - Internal CSCd to compute the norm of.
 *   comm   - MPI communicator.
 *
 * Returns:
 *   The norm 1 of the internal CSCd.
 */
#define z_CscNorm1 API_CALL(z_CscNorm1)
double z_CscNorm1(const z_CscMatrix *cscmtx,
                  MPI_Comm         comm);

/*
 * Function: z_CscbMAx
 *
 * Computes $$r = b-Ax$$.
 *
 * Parameters:
 *   sopalin_data - Internal structure containing global datas used by PaStiX.
 *   me           - thread ID.
 *   r            - will contains $$b-Ax$$
 *   b            - Vector $$b$$.
 *   cscmtx       - Internal CSCd matrix containing $$A$$.
 *   updovct      - z_UpDownVector structure containing $$x$$.
 *   solvmtx      - Solver matrix.
 *   comm         - MPI communicator.
 *   transpose    - Indicate if we want to transpose A.
 */
#define z_CscbMAx API_CALL(z_CscbMAx)
void z_CscbMAx(z_Sopalin_Data_t       *sopalin_data,
             int                   me,
             volatile pastix_complex64_t       *r,
             const volatile pastix_complex64_t *b,
             const z_CscMatrix      *cscmtx,
             const z_UpDownVector   *updovct,
             const z_SolverMatrix   *solvmtx,
             MPI_Comm              comm,
             pastix_int_t                   transpose);


/*
 * function: z_CscAxPb
 *
 *
 *  Compute the operation $$r=|A||x|+|b|$$
 *
 *  Parameters:
 *     sopalin_data - z_Sopalin_Data_t structure, common to all threads.
 *     me           - thread number
 *     r            - solution (vector commont to all threads)
 *     b            - Added vector (vector commont to all threads)
 *     cscmtx       - Compress Sparse Column matrix *A*
 *     updovct      - x, multiplied vector
 *     solvmtx      - solver matrix to know tho local structure of the matrix
 *     comm         - MPI communicator
 *     transpose    - Indicate if we want to transpose A.
 */
#define z_CscAxPb API_CALL(z_CscAxPb)
void z_CscAxPb(z_Sopalin_Data_t     *sopalin_data,
             int                 me,
             pastix_complex64_t              *r,
             const pastix_complex64_t        *b,
             const z_CscMatrix    *cscmtx,
             const z_UpDownVector *updovct,
             const z_SolverMatrix *solvmtx,
             MPI_Comm            comm,
             pastix_int_t                 transpose);


/*
 *  Function: z_CscBerr
 *
 *  Compute the operation $$berr= max_{i} |r|_{i}/|s|_{i}$$.
 *
 *  Parameters :
 *
 *  sopalin_data - z_Sopalin_Data_t structure, common to all threads.
 *  me           - thread number
 *  r            - vector(s) (common to all threads)
 *  s            - vector(s) (common to all threads)
 *  colnbr       - size of the vectors
 *  smvnbr       - number of vectors in r and s
 *  berr         - berr (local variable)
 *  comm         - MPI communicator
 */
#define z_CscBerr API_CALL(z_CscBerr)
void z_CscBerr(z_Sopalin_Data_t *sopalin_data,
             int            me,
             const pastix_complex64_t   *r,
             const pastix_complex64_t   *s,
             const pastix_int_t      colnbr,
             const pastix_int_t      smxnbr,
             double        *berr,
             MPI_Comm       comm);


/*
 * Function: z_CscNormErr
 *
 * Computes the norm 2 of r and the norm 2 of b and return
 * the quotient of these two vectors.
 *
 * This Function is multithreaded, each thread will compute a part of the norm,
 * it will be gathered between threads, then between MPI processors.
 *
 * Parameters:
 *   sopalin_data - global PaStix informations.
 *   me           - Thread ID.
 *   r            - first vector from which the norm 2 is computed.
 *   b            - second vector from which the norm 2 is computed.
 *   colnbr       - Size of the vectors.
 *   smxnbr       - Number of vectors (multi-right-hand-side method)
 *   comm         - PaStiX MPI communicator.
 */
#define z_CscNormErr API_CALL(z_CscNormErr)
double z_CscNormErr(z_Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_complex64_t *r,
                  const volatile pastix_complex64_t *b,
                  const pastix_int_t             colnbr,
                  const pastix_int_t             smxnbr,
                  MPI_Comm              comm);


/*
 * Function: z_CscNormFro
 *
 * Computes the norm 2 of x
 *
 * This Function is multithreaded, each thread will compute a part of the norm,
 * it will be gathered between threads, then between MPI processors.
 *
 * Parameters:
 *   sopalin_data - global PaStix informations.
 *   me           - Thread ID.
 *   x            - vector from which the norm 2 is computed.
 *   colnbr       - Size of the vectors.
 *   smxnbr       - Number of vectors (multi-right-hand-side method)
 *   comm         - PaStiX MPI communicator.
 */
#define z_CscNormFro API_CALL(z_CscNormFro)
double z_CscNormFro(z_Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_complex64_t *x,
                  const PASTIX_INT             colnbr,
                  const PASTIX_INT             smxnbr,
                  MPI_Comm              comm);



/*
 * Function: z_CscAx
 *
 * Computes *A* times *p* and store the result in *x*.
 *
 * When compiled with SMP_RAFF, this operation is multi-threaded.
 *
 * Parameters:
 *   sopalin_data - Gloabl PaStiX data.
 *   me           - thread ID.
 *   cscmtx       - Internal CSCd matrix, A.
 *   p            - Vector which will be multiplied by the CSCd.
 *   x            - vector which will contain the computation result.
 *   solvmtx      - Solver matrix.
 *   updovct      - Structure used for updown step, it contains information
 *                  about the vectors.
 *   comm         - MPI Communicator.
 *     transpose    - Indicate if we want to transpose A.
 */
#define z_CscAx API_CALL(z_CscAx)
void z_CscAx(z_Sopalin_Data_t       *sopalin_data,
           int                   me,
           const z_CscMatrix      *cscmtx,
           const volatile pastix_complex64_t *p,
           volatile pastix_complex64_t       *x,
           const z_SolverMatrix   *solvmtx,
           const z_UpDownVector   *updovct,
           MPI_Comm              comm,
           PASTIX_INT                   transpose);

/*
 * Function: z_CscGradBeta
 *
 * Computes the scalar product between *r* and *z*
 * and store the result in *beta*.
 *
 * At the end, beta is only on thread 0.
 *
 * Parameters:
 *   sopalin_data - PaStiX data structure.
 *   me           - Thread ID.
 *   r            - first vector of size *colnbr* times *smxnbr*.
 *   z            - second vector of size *colnbr* times *smxnbr*.a
 *   colnbr       - Number of unknowns.
 *   smxnbr       - Number of right-hand-side members.
 *   beta         - Float which will store the solution.
 *   comm         - MPI communicator.
 */
#define z_CscGradBeta API_CALL(z_CscGradBeta)
void z_CscGradBeta(z_Sopalin_Data_t       *sopalin_data,
                 int                   me,
                 const volatile pastix_complex64_t *r,
                 const volatile pastix_complex64_t *z,
                 pastix_int_t                   colnbr,
                 pastix_int_t                   smxnbr,
                 pastix_complex64_t               *beta,
                 MPI_Comm              comm);

/*
 * Function: z_CscGmresBeta
 *
 * Computes the scalar product between *r* and *z*
 * and store the result in *beta*.
 *
 * Parameters:
 *   sopalin_data - PaStiX data structure.
 *   me           - Thread ID.
 *   r            - first vector of size *colnbr* times *smxnbr*.
 *   z            - second vector of size *colnbr* times *smxnbr*.a
 *   colnbr       - Number of unknowns.
 *   smxnbr       - Number of right-hand-side members.
 *   beta         - Float which will store the solution.
 *   comm         - MPI communicator.
 */
#define z_CscGmresBeta API_CALL(z_CscGmresBeta)
void z_CscGmresBeta(z_Sopalin_Data_t                    *sopalin_data,
                    int                                me,
                    const volatile pastix_complex64_t *r,
                    const volatile pastix_complex64_t *z,
                    pastix_int_t                       colnbr,
                    pastix_int_t                       smxnbr,
                    pastix_complex64_t                *beta,
                    MPI_Comm                           comm);


/*
 * Function: z_CscCopy
 *
 * Copy a vector into another vector
 *
 * This Function is multithreaded, each thread will compute a part of the copy,
 * it will be gathered between threads, then between MPI processors.
 *
 * Parameters:
 *   sopalin_data - global PaStix informations.
 *   me           - Thread ID.
 *   x            - vector from which the copy is done.
 *   y            - vector where the copy is done
 *   colnbr       - Size of the vectors.
 *   smxnbr       - Number of vectors (multi-right-hand-side method)
 *   comm         - PaStiX MPI communicator.
 */
#define z_CscCopy API_CALL(z_CscCopy)
void z_CscCopy(z_Sopalin_Data_t              *sopalin_data,
             int                          me,
             const volatile pastix_complex64_t *x,
             volatile pastix_complex64_t       *y,
             const PASTIX_INT             colnbr,
             const PASTIX_INT             smxnbr,
             MPI_Comm                     comm);

/*
 * Function: z_CscScal
 *
 * Multiply a vector by a scalaire
 *
 * This Function is multithreaded, each thread will compute a part of the copy,
 * it will be gathered between threads, then between MPI processors.
 *
 * Parameters:
 *   sopalin_data - global PaStix informations.
 *   me           - Thread ID.
 *   x            - vector from which the copy is done.
 *   y            - vector where the copy is done
 *   colnbr       - Size of the vectors.
 *   smxnbr       - Number of vectors (multi-right-hand-side method)
 *   comm         - PaStiX MPI communicator.
 */
#define z_CscScal API_CALL(z_CscScal)
void z_CscScal(z_Sopalin_Data_t        *sopalin_data,
             int                    me,
             volatile pastix_complex64_t  alpha,
             volatile pastix_complex64_t *x,
             const PASTIX_INT       colnbr,
             const PASTIX_INT       smxnbr,
             MPI_Comm               comm);


/*
 * Function: z_CscAXPY
 *
 * Y<-aX+Y
 *
 * This Function is multithreaded, each thread will compute a part of the operation,
 * it will be gathered between threads, then between MPI processors.
 *
 * Parameters:
 *   sopalin_data - global PaStix informations.
 *   me           - Thread ID.
 *   alpha
 *   x
 *   y
 *   colnbr       - Size of the vectors.
 *   smxnbr       - Number of vectors (multi-right-hand-side method)
 *   comm         - PaStiX MPI communicator.
 */
#define z_CscAXPY API_CALL(z_CscAXPY)
void z_CscAXPY(z_Sopalin_Data_t              *sopalin_data,
             int                          me,
             pastix_complex64_t                 alpha,
             const volatile pastix_complex64_t *x,
             volatile pastix_complex64_t       *y,
             const PASTIX_INT             colnbr,
             const PASTIX_INT             smxnbr,
             MPI_Comm                     comm);

#endif /* Z_CSC_INTERN_COMPUTE_H */
