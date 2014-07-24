/*
 * File: csc_intern_compute.h
 *
 * Functions computing operations on the CSC.
 *
 */

#ifndef CSC_INTERN_COMPUTE_H
#define CSC_INTERN_COMPUTE_H

/*
 * Function: CscNorm1
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
#define CscNorm1 API_CALL(CscNorm1)
double CscNorm1(const CscMatrix *cscmtx,
                MPI_Comm         comm);

/*
 * Function: CscbMAx
 *
 * Computes $$r = b-Ax$$.
 *
 * Parameters:
 *   sopalin_data - Internal structure containing global datas used by PaStiX.
 *   me           - thread ID.
 *   r            - will contains $$b-Ax$$
 *   b            - Vector $$b$$.
 *   cscmtx       - Internal CSCd matrix containing $$A$$.
 *   updovct      - d_UpDownVector structure containing $$x$$.
 *   solvmtx      - Solver matrix.
 *   comm         - MPI communicator.
 *   transpose    - Indicate if we want to transpose A.
 */
#define CscbMAx API_CALL(CscbMAx)
void CscbMAx(Sopalin_Data_t       *sopalin_data,
             int                   me,
             volatile pastix_float_t       *r,
             const volatile pastix_float_t *b,
             const CscMatrix      *cscmtx,
             const d_UpDownVector   *updovct,
             const d_SolverMatrix   *solvmtx,
             MPI_Comm              comm,
             pastix_int_t                   transpose);


/*
 * function: CscAxPb
 *
 *
 *  Compute the operation $$r=|A||x|+|b|$$
 *
 *  Parameters:
 *     sopalin_data - Sopalin_Data_t structure, common to all threads.
 *     me           - thread number
 *     r            - solution (vector commont to all threads)
 *     b            - Added vector (vector commont to all threads)
 *     cscmtx       - Compress Sparse Column matrix *A*
 *     updovct      - x, multiplied vector
 *     solvmtx      - solver matrix to know tho local structure of the matrix
 *     comm         - MPI communicator
 *     transpose    - Indicate if we want to transpose A.
 */
#define CscAxPb API_CALL(CscAxPb)
void CscAxPb(Sopalin_Data_t     *sopalin_data,
             int                 me,
             pastix_float_t              *r,
             const pastix_float_t        *b,
             const CscMatrix    *cscmtx,
             const d_UpDownVector *updovct,
             const d_SolverMatrix *solvmtx,
             MPI_Comm            comm,
             pastix_int_t                 transpose);


/*
 *  Function: CscBerr
 *
 *  Compute the operation $$berr= max_{i} |r|_{i}/|s|_{i}$$.
 *
 *  Parameters :
 *
 *  sopalin_data - Sopalin_Data_t structure, common to all threads.
 *  me           - thread number
 *  r            - vector(s) (common to all threads)
 *  s            - vector(s) (common to all threads)
 *  colnbr       - size of the vectors
 *  smvnbr       - number of vectors in r and s
 *  berr         - berr (local variable)
 *  comm         - MPI communicator
 */
#define CscBerr API_CALL(CscBerr)
void CscBerr(Sopalin_Data_t *sopalin_data,
             int            me,
             const pastix_float_t   *r,
             const pastix_float_t   *s,
             const pastix_int_t      colnbr,
             const pastix_int_t      smxnbr,
             double        *berr,
             MPI_Comm       comm);


/*
 * Function: CscNormErr
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
#define CscNormErr API_CALL(CscNormErr)
double CscNormErr(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_float_t *r,
                  const volatile pastix_float_t *b,
                  const pastix_int_t             colnbr,
                  const pastix_int_t             smxnbr,
                  MPI_Comm              comm);


/*
 * Function: CscNormFro
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
#define CscNormFro API_CALL(CscNormFro)
double CscNormFro(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_float_t *x,
                  const PASTIX_INT             colnbr,
                  const PASTIX_INT             smxnbr,
                  MPI_Comm              comm);



/*
 * Function: CscAx
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
#define CscAx API_CALL(CscAx)
void CscAx(Sopalin_Data_t       *sopalin_data,
           int                   me,
           const CscMatrix      *cscmtx,
           const volatile pastix_float_t *p,
           volatile pastix_float_t       *x,
           const d_SolverMatrix   *solvmtx,
           const d_UpDownVector   *updovct,
           MPI_Comm              comm,
           PASTIX_INT                   transpose);

/*
 * Function: CscGradBeta
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
#define CscGradBeta API_CALL(CscGradBeta)
void CscGradBeta(Sopalin_Data_t       *sopalin_data,
                 int                   me,
                 const volatile pastix_float_t *r,
                 const volatile pastix_float_t *z,
                 pastix_int_t                   colnbr,
                 pastix_int_t                   smxnbr,
                 pastix_float_t               *beta,
                 MPI_Comm              comm);

/*
 * Function: CscGmresBeta
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
#define CscGmresBeta API_CALL(CscGmresBeta)
void CscGmresBeta(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_float_t *r,
                  const volatile pastix_float_t *z,
                  pastix_int_t                   colnbr,
                  pastix_int_t                   smxnbr,
                  double               *beta,
                  MPI_Comm              comm);


/*
 * Function: CscCopy
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
#define CscCopy API_CALL(CscCopy)
void CscCopy(Sopalin_Data_t              *sopalin_data,
             int                          me,
             const volatile pastix_float_t *x,
             volatile pastix_float_t       *y,
             const PASTIX_INT             colnbr,
             const PASTIX_INT             smxnbr,
             MPI_Comm                     comm);

/*
 * Function: CscScal
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
#define CscScal API_CALL(CscScal)
void CscScal(Sopalin_Data_t        *sopalin_data,
             int                    me,
             volatile pastix_float_t  alpha,
             volatile pastix_float_t *x,
             const PASTIX_INT       colnbr,
             const PASTIX_INT       smxnbr,
             MPI_Comm               comm);


/*
 * Function: CscAXPY
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
#define CscAXPY API_CALL(CscAXPY)
void CscAXPY(Sopalin_Data_t              *sopalin_data,
             int                          me,
             pastix_float_t                 alpha,
             const volatile pastix_float_t *x,
             volatile pastix_float_t       *y,
             const PASTIX_INT             colnbr,
             const PASTIX_INT             smxnbr,
             MPI_Comm                     comm);

#endif /* CSC_INTERN_COMPUTE_H */
