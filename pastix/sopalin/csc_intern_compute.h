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
 *   updovct      - UpDownVector structure containing $$x$$.
 *   solvmtx      - Solver matrix.
 *   comm         - MPI communicator.
 *   transpose    - Indicate if we want to transpose A.
 */
#define CscbMAx API_CALL(CscbMAx)
void CscbMAx(Sopalin_Data_t       *sopalin_data,
             int                   me,
             volatile PASTIX_FLOAT       *r,
             const volatile PASTIX_FLOAT *b,
             const CscMatrix      *cscmtx,
             const UpDownVector   *updovct,
             const SolverMatrix   *solvmtx,
             MPI_Comm              comm,
             PASTIX_INT                   transpose);


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
             PASTIX_FLOAT              *r,
             const PASTIX_FLOAT        *b,
             const CscMatrix    *cscmtx,
             const UpDownVector *updovct,
             const SolverMatrix *solvmtx,
             MPI_Comm            comm,
             PASTIX_INT                 transpose);


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
             const PASTIX_FLOAT   *r,
             const PASTIX_FLOAT   *s,
             const PASTIX_INT      colnbr,
             const PASTIX_INT      smxnbr,
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
                  const volatile PASTIX_FLOAT *r,
                  const volatile PASTIX_FLOAT *b,
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
           const volatile PASTIX_FLOAT *p,
           volatile PASTIX_FLOAT       *x,
           const SolverMatrix   *solvmtx,
           const UpDownVector   *updovct,
           MPI_Comm              comm,
           PASTIX_INT                   transpose);

/*
 * Function: CscGradAlpha
 *
 * Computes the scalar product of *r* with *z*,
 * then the scalar product of *x* with *p*
 * and finaly store the quotient in *alpha*.
 *
 * Multi-threaded in SMP_RAFF mode.
 *
 * Parameters:
 *   sopalin_data - Gloabal PaStiX data structure.
 *   me           - Thread ID.
 *   r            - A vector of size *colnbr* times *smxnbr*.
 *   z            - A vector of size *colnbr* times *smxnbr*.
 *   x            - A vector of size *colnbr* times *smxnbr*.
 *   p            - A vector of size *colnbr* times *smxnbr*.
 *   colnbr       - Number of unkowns.
 *   smxnbr       - Number of right-hand-side members.
 *   alpha        - Float which will store the computation result.
 *   comm         - MPI communicator.
 */
#define CscGradAlpha API_CALL(CscGradAlpha)
void CscGradAlpha(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile PASTIX_FLOAT *r,
                  const volatile PASTIX_FLOAT *z,
                  const volatile PASTIX_FLOAT *x,
                  const volatile PASTIX_FLOAT *p,
                  PASTIX_INT                   colnbr,
                  PASTIX_INT                   smxnbr,
                  double               *alpha,
                  MPI_Comm              comm);

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
                 const volatile PASTIX_FLOAT *r,
                 const volatile PASTIX_FLOAT *z,
                 PASTIX_INT                   colnbr,
                 PASTIX_INT                   smxnbr,
                 PASTIX_FLOAT               *beta,
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
                  const volatile PASTIX_FLOAT *r,
                  const volatile PASTIX_FLOAT *z,
                  PASTIX_INT                   colnbr,
                  PASTIX_INT                   smxnbr,
                  double               *beta,
                  MPI_Comm              comm);
#endif /* CSC_INTERN_COMPUTE_H */
