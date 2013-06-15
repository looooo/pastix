/*
 * File: read_matrix.h
 *
 * Definition of a global function to read all type of matrices.
 *
 */
/*
 * typedef: driver_type_enum
 *
 * RSA        - rsa file format.
 * CHB        - chb file format.
 * CCC        - ccc file format.
 * RCC        - rcc file format.
 * OLAF       - olaf file format.
 * PEER       - peer file format.
 * HB         - Harwell-boeing file format.
 * THREEFILES - IJV file format with one file for each element of the triple.
 * MM         - Matrix Market file format.
 * PETSCS     - PETSc binary file format symmetric.
 * PETSCH     - PETSc binary file format hermitian.
 * PETSCU     - PETSc binary file format unsymmetric.
 * CSCD       - CSCd file format.
 * LAPLACIAN  - Generating laplacian.
 * FDUP       - Binary matrix driver from Fabrice dupros.
 */
#ifndef READ_MATRIX_H
#define READ_MATRIX_H
enum driver_type_enum {
  RSA,
  CHB,
  CCC,
  RCC,
  OLAF,
  PEER,
  HB,
  THREEFILES,
  MM,
  MMD,
  PETSCS,
  PETSCU,
  PETSCH,
  CSCD,
  LAPLACIAN,
  FDUP,
  FDUP_DIST
};
/*
  typedef: driver_type_t
 */
typedef enum driver_type_enum driver_type_t;

#ifndef MTX_ISSYM
#  define MTX_ISSYM(a) ((a)[1]=='S')
#  define MTX_ISHER(a) ((a)[1]=='H')
#  define MTX_ISCOM(a) ((a)[0]=='C')
#  define MTX_ISRHX(a) ((a)[2]=='X')
#  define MTX_ISRHS(a) ((a)[0]!='\0')
#endif

int api_str_to_int(char * string, int * value);

/*
 * Function: read_matrix
 *
 * Reads a matrix from a file in the format given by driver_type
 *  (see <driver_type_enum>).
 *
 * Parameters:
 *   filename    - Name of the file to read from.
 *   ncol        - Number of column in the matrix (output).
 *   colptr      - Indexes in row and avals of first element of each column
 *                 of the matrix.(output)
 *   row         - Row of each element of the matrix.
 *   values      - Values of each element of the matrix.
 *   type        - type of the matrix.
 *   rhstype     - type of the right and side.
 *   driver_type - driver to use to read the matrix.
 *   pastix_comm - MPI communicator containing all processes wich call
 *                 <read_matrix>.
 */
int read_matrix(char            *filename,
                pastix_int_t    *ncol,
                pastix_int_t   **colptr,
                pastix_int_t   **row,
                pastix_float_t **values,
                pastix_float_t **rhs,
                char           **type,
                char           **rhstype,
                driver_type_t    driver_type,
                MPI_Comm         pastix_comm);

/*
 * Function: dread_matrix
 *
 * Reads a matrix from a file in the format given by driver_type
 * (see <driver_type_enum>) and distribute the matrix on the MPI processors.
 *
 * Parameters:
 *   filename    - Name of the file to read from.
 *   ncol        - Number of column in the matrix (output).
 *   colptr      - Indexes in rows and avals of first element of each column
 *                 of the matrix.(output)
 *   rows        - Row of each element of the matrix.
 *   values      - Values of each element of the matrix.
 *   type        - type of the matrix.
 *   rhstype     - type of the right and side.
 *   driver_type - driver to use to read the matrix.
 *   pastix_comm - MPI communicator containing all processes wich call
 *                 <read_matrix>.
 */
int dread_matrix(char            *filename,
                 pastix_int_t    *ncol,
                 pastix_int_t   **colptr,
                 pastix_int_t   **rows,
                 pastix_int_t   **loc2glob,
                 pastix_float_t **values,
                 pastix_float_t **rhs,
                 char           **type,
                 char           **rhstype,
                 driver_type_t    driver_type,
                 MPI_Comm         pastix_comm);

/*
 *  Function: comparcouple
 *
 *  Compares 2 couples (i,j).
 *
 *  Parameters:
 *    a - one couple
 *    b - one other couple
 *
 *  Returns:
 *    -1 - if a_i < b_i
 *    1  - if a_i > b_i
 *    -1 - if a_i = b_i and a_j < b_j
 *    1  - else
 */
int comparcouple(const void *a, const void *b);

/*
 *  Function: checkStrucSym
 *
 *  Check if unsymmetric matrix has a symmetric pattern, if not correct it.
 *
 *  Parameters:
 *    n      - size of the matrix
 *    nz     - number of elements
 *    colptr - Index of first element of each column in *row* and *avals*
 *    row    - row of each element
 *    avals  - value of each element
 */
void checkStrucSym(pastix_int_t     n,
                   pastix_int_t    *nz,
                   pastix_int_t   **colptr,
                   pastix_int_t   **row,
                   pastix_float_t **avals);

#endif /* not READ_MATRIX_H */
