/*
  File: cscd_utils.h

  Several operations on CSCD.

 */

#ifndef CSCD_UTILS_H
#define CSCD_UTILS_H

#ifndef _GLIBCXX_HAVE_COMPLEX_H
#  define _GLIBCXX_HAVE_COMPLEX_H 0
#endif

#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX || defined _STLP_template_complex)
#  define PASTIX_HAS_COMPLEX
#endif

#ifdef   __cplusplus
#  if (_GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
#    define  COMPLEX  std::complex<float>
#    define  DCOMPLEX std::complex<double>
#  endif
#else /* not __cplusplus */
#  if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__)
#    define  COMPLEX float complex
#    define  DCOMPLEX double complex
#  endif
#endif /* not __cplusplus */

/*
 * MULTIPLE_TYPE_DEFINE
 *
 * Automaticaly generate function for each pastix_float_t type.
 *
 * This macro is fitted for function not taking floating points arguments.
 */
#ifdef PASTIX_HAS_COMPLEX
#define MULTIPLE_TYPE_DEFINE(functype, funcname, funcargs) \
  functype s_ ## funcname funcargs;                        \
  functype d_ ## funcname funcargs;                        \
  functype c_ ## funcname funcargs;                        \
  functype z_ ## funcname funcargs;
#else
#define MULTIPLE_TYPE_DEFINE(functype, funcname, funcargs)  \
  functype s_ ## funcname funcargs;                         \
  functype d_ ## funcname funcargs;
#endif

/*
 * MULTIPLE_TYPE_DEFINE_F
 *
 * Automaticaly generate function for each pastix_float_t type.
 *
 * This macro is fitted for function taking floating points arguments.
 */
#ifdef PASTIX_HAS_COMPLEX
#define MULTIPLE_TYPE_DEFINE_F(functype,        \
                               funcname,        \
                               funcargs_s,			\
                               funcargs_d,			\
                               funcargs_c,			\
                               funcargs_z)			\
  functype s_ ## funcname funcargs_s;           \
  functype d_ ## funcname funcargs_d;           \
  functype c_ ## funcname funcargs_c;           \
  functype z_ ## funcname funcargs_z;
#else
#define MULTIPLE_TYPE_DEFINE_F(functype,        \
                               funcname,        \
                               funcargs_s,			\
                               funcargs_d,			\
                               funcargs_c,			\
                               funcargs_z)			\
  functype s_ ## funcname funcargs_s;           \
  functype d_ ## funcname funcargs_d;
#endif

/*
 * Enum: CSCD_OPERATIONS
 *
 * Operation when adding CSCD
 *
 * CSCD_ADD  - Add coefficient values.
 * CSCD_KEEP - Keep value from first CSCD.
 * CSCD_MAX  - Keep maximum of first and second CSCD.
 * CSCD_MIN  - Keep minimum of first and second CSCD.
 * CSCD_OVW  - Overwrite with second CSCD value.
 */
enum CSCD_OPERATIONS {
  CSCD_ADD,
  CSCD_KEEP,
  CSCD_MAX,
  CSCD_MIN,
  CSCD_OVW
};
typedef enum CSCD_OPERATIONS CSCD_OPERATIONS_t;

/*
 * Enum: CSC_DISPATCH_OP
 *
 * Operation when dispatching the csc into a cscd
 *
 * CSC_DISP_SIMPLE - Reparts linearly the columns over the proc.
 * CSC_DISP_CYCLIC - Reparts cyclicly the columns over the proc.
 *
 */
enum CSC_DISPATCH_OP {
  CSC_DISP_SIMPLE,
  CSC_DISP_CYCLIC
};
typedef enum CSC_DISPATCH_OP CSCD_DISPATCH_OP_t;

/* Section: Functions */
#if (defined pastix_float_t)
/*
 *  Function: csc_dispatch
 *
 *  Distribute a CSC to a CSCD
 *
 *  Parameters:
 *     gN                - global number of columns
 *     gcolptr           - global starting index of each column in grows ans gavals.
 *     grows             - global rows of each element.
 *     gavals            - global values of each element.
 *     gperm             - global permutation tabular.
 *     ginvp             - global reverse permutation tabular.
 *     lN                - local number of columns (output).
 *     lcolptr           - starting index of each local column (output).
 *     lrowptr           - row number of each local element (output).
 *     lavals            - values of each local element (output).
 *     lrhs              - local part of the right hand side (output).
 *     lperm             - local part of the permutation tabular (output).
 *     loc2glob          - global numbers of local columns (before permutation).
 *     dispatch          - choose how to dispatch the csc
 *     pastix_comm       - PaStiX MPI communicator.
 */
void csc_dispatch(pastix_int_t  gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, pastix_float_t *  gavals,
                  pastix_float_t *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                  pastix_int_t *lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, pastix_float_t ** lavals,
                  pastix_float_t ** lrhs, pastix_int_t ** lperm,
                  pastix_int_t **loc2glob, int dispatch, MPI_Comm pastix_comm);
#endif
MULTIPLE_TYPE_DEFINE_F(void,
                       csc_dispatch,
                       (pastix_int_t  gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, float *  gavals,
                        float *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t *lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, float ** lavals,
                        float ** lrhs, pastix_int_t ** lperm,
                        pastix_int_t **loc2glob, int dispatch, MPI_Comm pastix_comm),
                       (pastix_int_t  gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, double *  gavals,
                        double *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t *lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, double ** lavals,
                        double ** lrhs, pastix_int_t ** lperm,
                        pastix_int_t **loc2glob, int dispatch, MPI_Comm pastix_comm),
                       (pastix_int_t  gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, COMPLEX *  gavals,
                        COMPLEX *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t *lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, COMPLEX ** lavals,
                        COMPLEX ** lrhs, pastix_int_t ** lperm,
                        pastix_int_t **loc2glob, int dispatch, MPI_Comm pastix_comm),
                       (pastix_int_t  gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, DCOMPLEX *  gavals,
                        DCOMPLEX *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t *lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, DCOMPLEX ** lavals,
                        DCOMPLEX ** lrhs, pastix_int_t ** lperm,
                        pastix_int_t **loc2glob, int dispatch, MPI_Comm pastix_comm))

#if (defined pastix_float_t)
/*
 * Function: csc_cyclic_distribution
 *
 * Distribute the CSC cyclicaly.
 *
 * Parameters:
 *   column      - column number to distribute
 *   columnnbr   - Number of colmuns.
 *   pastix_comm - PaStiX MPI communicator
 *
 * Return:
 *   owner of the column (column%commSize)
 */
pastix_int_t csc_cyclic_distribution(pastix_int_t column, pastix_int_t columnnbr, MPI_Comm pastix_comm);
#endif
MULTIPLE_TYPE_DEFINE(pastix_int_t, csc_cyclic_distribution,
                     (pastix_int_t column, pastix_int_t columnnbr, MPI_Comm pastix_comm))

#if (defined pastix_float_t)
/*
 * Function: csc_simple_distribution
 *
 * Distribute the CSC.
 * First columns are for first proc and so on.
 *
 * Parameters:
 *   column      - column number to distribute
 *   columnnbr   - Number of colmuns.
 *   pastix_comm - PaStiX MPI communicator
 *
 * Return:
 *   owner of the column (column/commSize)
 */
pastix_int_t csc_simple_distribution(pastix_int_t column, pastix_int_t columnnbr, MPI_Comm pastix_comm);

#endif
MULTIPLE_TYPE_DEFINE(pastix_int_t, csc_simple_distribution,
                     (pastix_int_t column, pastix_int_t columnnbr, MPI_Comm pastix_comm))

#if (defined pastix_float_t)
/*
 * Function: cscd_symgraph
 *
 * Check if the CSCD graph is symetric.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - Starting index of each columns in *ja* and *a*
 *   ja          - Row of each element.
 *   a           - Values of each element.
 *   newn        - New number of local columns
 *   newia       - Starting index of each columns in *newja* and *newa*
 *   newja       - Row of each element.
 *   newa        - Values of each element.
 *   l2g         - global number of each local column.
 *   malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int cscd_symgraph(pastix_int_t      n, pastix_int_t *     ia, pastix_int_t *     ja, pastix_float_t *     a,
                  pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja, pastix_float_t ** newa,
                  pastix_int_t *     l2g,  MPI_Comm comm);
#endif
MULTIPLE_TYPE_DEFINE(int, cscd_symgraph,
                     (pastix_int_t      n, pastix_int_t *     ia, pastix_int_t *     ja, float *     a,
                      pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja, float ** newa,
                      pastix_int_t *     l2g,  MPI_Comm comm))

#if (defined pastix_float_t)
/*
 * Function: cscd_addlocal
 *
 * Add second cscd to first cscd into third cscd (unallocated)
 *
 * Parameters:
 *   n           - First cscd size
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   l2g         - local 2 global column numbers for first cscd
 *   addn        - CSCD to add size
 *   addia       - CSCD to add starting index of each column in *addja* and *adda*
 *   addja       - Row of each element in second CSCD
 *   adda        - value of each cscd in second CSCD (can be NULL -> add 0)
 *   addl2g      - local 2 global column numbers for second cscd
 *   newn        - new cscd size (same as first)
 *   newia       - CSCD to add starting index of each column in *newja* and *newwa*
 *   newja       - Row of each element in third CSCD
 *   newa        - value of each cscd in third CSCD
 *   OP          - Operation to manage common CSCD coefficients.
 *   dof         - Number of degrees of freedom.
 */

int cscd_addlocal(pastix_int_t   n   , pastix_int_t *  ia   , pastix_int_t *  ja   , pastix_float_t *  a   , pastix_int_t * l2g,
      pastix_int_t   addn, pastix_int_t *  addia, pastix_int_t *  addja, pastix_float_t *  adda, pastix_int_t * addl2g,
      pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja, pastix_float_t ** newa, CSCD_OPERATIONS_t OP, int dof);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_addlocal,
                       (pastix_int_t   n   , pastix_int_t *  ia   , pastix_int_t *  ja   , float *  a   , pastix_int_t * l2g,
                        pastix_int_t   addn, pastix_int_t *  addia, pastix_int_t *  addja, float *  adda, pastix_int_t * addl2g,
                        pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja, float ** newa, CSCD_OPERATIONS_t OP, int dof),
                       (pastix_int_t   n   , pastix_int_t *  ia   , pastix_int_t *  ja   , double *  a   , pastix_int_t * l2g,
                        pastix_int_t   addn, pastix_int_t *  addia, pastix_int_t *  addja, double *  adda, pastix_int_t * addl2g,
                        pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja, double ** newa, CSCD_OPERATIONS_t OP, int dof),
                       (pastix_int_t   n   , pastix_int_t *  ia   , pastix_int_t *  ja,
                        COMPLEX *  a   , pastix_int_t * l2g,
                        pastix_int_t   addn, pastix_int_t *  addia, pastix_int_t *  addja,
                        COMPLEX *  adda, pastix_int_t * addl2g,
                        pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja,
                        COMPLEX ** newa, CSCD_OPERATIONS_t OP, int dof),
                       (pastix_int_t   n   , pastix_int_t *  ia   , pastix_int_t *  ja,
                        DCOMPLEX *  a   , pastix_int_t * l2g,
                        pastix_int_t   addn, pastix_int_t *  addia, pastix_int_t *  addja,
                        DCOMPLEX *  adda, pastix_int_t * addl2g,
                        pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja,
                        DCOMPLEX ** newa, CSCD_OPERATIONS_t OP, int dof))


#if (defined pastix_float_t)
/**
 *   Function: csc2cscd
 *
 *   Transform a csc to a cscd.
 *   Allocate the CSCD.
 *   If grhs == NULL forget right hand side part.
 *   If gperm == NULL forget permutation and reverse permutation part.
 *
 *   Parameters:
 *     gN       - global number of columns
 *     gcolptr  - global starting index of each column in grows ans gavals.
 *     grows    - global rows of each element.
 *     gavals   - global values of each element.
 *     gperm    - global permutation tabular.
 *     ginvp    - global reverse permutation tabular.
 *     lN       - local number of columns.
 *     lcolptr  - starting index of each local column.
 *     lrowptr  - row number of each local element.
 *     lavals   - values of each local element.
 *     lrhs     - local part of the right hand side (output).
 *     lperm    - local part of the permutation tabular (output).
 *     linvp    - local part of the reverse permutation tabular (output).
 *     loc2glob - global numbers of local columns (before permutation).
 */
void  csc2cscd(pastix_int_t gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, pastix_float_t *  gavals,
               pastix_float_t *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
               pastix_int_t lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, pastix_float_t ** lavals,
               pastix_float_t ** lrhs, pastix_int_t ** lperm, pastix_int_t ** linvp,
               pastix_int_t *loc2glob);
#endif
MULTIPLE_TYPE_DEFINE_F(void,  csc2cscd,
                       (pastix_int_t gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, float *  gavals,
                        float *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, float ** lavals,
                        float ** lrhs, pastix_int_t ** lperm, pastix_int_t ** linvp,
                        pastix_int_t *loc2glob),
                       (pastix_int_t gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, double *  gavals,
                        double *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, double ** lavals,
                        double ** lrhs, pastix_int_t ** lperm, pastix_int_t ** linvp,
                        pastix_int_t *loc2glob),
                       (pastix_int_t gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, COMPLEX *  gavals,
                        COMPLEX *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, COMPLEX ** lavals,
                        COMPLEX ** lrhs, pastix_int_t ** lperm, pastix_int_t ** linvp,
                        pastix_int_t *loc2glob),
                       (pastix_int_t gN, pastix_int_t *  gcolptr, pastix_int_t *  grow, DCOMPLEX *  gavals,
                        DCOMPLEX *  grhs, pastix_int_t *  gperm, pastix_int_t *  ginvp,
                        pastix_int_t lN, pastix_int_t ** lcolptr, pastix_int_t ** lrow, DCOMPLEX ** lavals,
                        DCOMPLEX ** lrhs, pastix_int_t ** lperm, pastix_int_t ** linvp,
                        pastix_int_t *loc2glob))

#if (defined pastix_float_t)
/**
 *   Function: cscd2csc
 *
 *   Transform a cscd to a csc.
 *   colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.
 *
 *   Parameters:
 *      lN          - number of local column.
 *      lcolptr     - starting index of each local column in row and avals.
 *      lrow        _ row number of each local element.
 *      lavals      - values of each local element.
 *      lrhs        - local part of the right hand side.
 *      lperm       - local part of the permutation tabular.
 *      linvp       - local part of the reverse permutation tabular.
 *      gN          - global number of columns (output).
 *      gcolptr     - starting index of each column in row2 and avals2 (output).
 *      grow        - row number of each element (output).
 *      gavals      - values of each element (output).
 *      grhs        - global right hand side (output).
 *      gperm       - global permutation tabular (output).
 *      ginvp       - global reverse permutation tabular (output).
 *      loc2glob    - global number of each local column.
 *      pastix_comm - PaStiX MPI communicator.
 *      ndof        - Number of degree f freedom by node.
 *
 */

void  cscd2csc(pastix_int_t  lN, pastix_int_t *  lcolptr, pastix_int_t * lrow, pastix_float_t * lavals,
               pastix_float_t * lrhs, pastix_int_t * lperm, pastix_int_t * linvp,
               pastix_int_t *gN, pastix_int_t ** gcolptr, pastix_int_t **grow, pastix_float_t **gavals,
               pastix_float_t **grhs, pastix_int_t **gperm, pastix_int_t **ginvp,
               pastix_int_t *loc2glob, MPI_Comm pastix_comm, pastix_int_t ndof);
#endif
MULTIPLE_TYPE_DEFINE_F(void,  cscd2csc,
                       (pastix_int_t  lN, pastix_int_t *  lcolptr, pastix_int_t * lrow, float * lavals,
                        float * lrhs, pastix_int_t * lperm, pastix_int_t * linvp,
                        pastix_int_t *gN, pastix_int_t ** gcolptr, pastix_int_t **grow, float **gavals,
                        float **grhs, pastix_int_t **gperm, pastix_int_t **ginvp,
                        pastix_int_t *loc2glob, MPI_Comm pastix_comm, pastix_int_t ndof),
                       (pastix_int_t  lN, pastix_int_t *  lcolptr, pastix_int_t * lrow, double * lavals,
                        double * lrhs, pastix_int_t * lperm, pastix_int_t * linvp,
                        pastix_int_t *gN, pastix_int_t ** gcolptr, pastix_int_t **grow, double **gavals,
                        double **grhs, pastix_int_t **gperm, pastix_int_t **ginvp,
                        pastix_int_t *loc2glob, MPI_Comm pastix_comm, pastix_int_t ndof),
                       (pastix_int_t  lN, pastix_int_t *  lcolptr, pastix_int_t * lrow, COMPLEX * lavals,
                        COMPLEX * lrhs, pastix_int_t * lperm, pastix_int_t * linvp,
                        pastix_int_t *gN, pastix_int_t ** gcolptr, pastix_int_t **grow, COMPLEX **gavals,
                        COMPLEX **grhs, pastix_int_t **gperm, pastix_int_t **ginvp,
                        pastix_int_t *loc2glob, MPI_Comm pastix_comm, pastix_int_t ndof),
                       (pastix_int_t  lN, pastix_int_t *  lcolptr, pastix_int_t * lrow, DCOMPLEX * lavals,
                        DCOMPLEX * lrhs, pastix_int_t * lperm, pastix_int_t * linvp,
                        pastix_int_t *gN, pastix_int_t ** gcolptr, pastix_int_t **grow, DCOMPLEX **gavals,
                        DCOMPLEX **grhs, pastix_int_t **gperm, pastix_int_t **ginvp,
                        pastix_int_t *loc2glob, MPI_Comm pastix_comm, pastix_int_t ndof))

#if (defined pastix_float_t)
/*
 * Function: cscd_redispatch
 *
 * Redistribute the first cscd into a new one using *dl2g*.
 *
 * - gather all new loc2globs on all processors.
 * - allocate *dia*, *dja* and *da*.
 * - Create new CSC for each processor and send it.
 * - Merge all new CSC to the new local CSC with <cscd_addlocal_int>.
 *
 * If communicator size is one, check that n = dn and
 * l2g = dl2g and simply create a copy of the first cscd.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   rhs         - right-hand-side member corresponding to the first CSCD
 *                 (can be NULL)
 *   nrhs        - number of right-hand-side.
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - New cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in new CSCD
 *   da          - value of each cscd in new CSCD
 *   rhs         - right-hand-side member corresponding to the new CSCD
 *   dl2g        - local 2 global column numbers for new cscd
 *   comm        - MPI communicator
 *
 * Returns:
 *   EXIT_SUCCESS - If all goes well
 *   EXIT_FAILURE - If commsize = 1 and *n* != *dn* or *l2g* != *dl2g*.
 */
int cscd_redispatch(pastix_int_t   n, pastix_int_t *   ia, pastix_int_t *   ja, pastix_float_t *   a,
                    pastix_float_t *  rhs,  pastix_int_t nrhs, pastix_int_t *   l2g,
                    pastix_int_t  dn, pastix_int_t ** dia, pastix_int_t ** dja, pastix_float_t ** da,
                    pastix_float_t ** drhs,  pastix_int_t *  dl2g,
                    MPI_Comm comm, pastix_int_t dof);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_redispatch,
                       (pastix_int_t   n, pastix_int_t *   ia, pastix_int_t *   ja, float *   a,
                        float *  rhs,  pastix_int_t nrhs,   pastix_int_t *   l2g,
                        pastix_int_t  dn, pastix_int_t ** dia, pastix_int_t ** dja, float ** da,
                        float ** drhs,  pastix_int_t *  dl2g,
                        MPI_Comm comm, pastix_int_t dof),
                       (pastix_int_t   n, pastix_int_t *   ia, pastix_int_t *   ja, double *   a,
                        double *  rhs,  pastix_int_t nrhs,   pastix_int_t *   l2g,
                        pastix_int_t  dn, pastix_int_t ** dia, pastix_int_t ** dja, double ** da,
                        double ** drhs,  pastix_int_t *  dl2g,
                        MPI_Comm comm, pastix_int_t dof),
                       (pastix_int_t   n, pastix_int_t *   ia, pastix_int_t *   ja, COMPLEX *   a,
                        COMPLEX *  rhs,  pastix_int_t nrhs,   pastix_int_t *   l2g,
                        pastix_int_t  dn, pastix_int_t ** dia, pastix_int_t ** dja, COMPLEX ** da,
                        COMPLEX ** drhs,  pastix_int_t *  dl2g,
                        MPI_Comm comm, pastix_int_t dof),
                       (pastix_int_t   n, pastix_int_t *   ia, pastix_int_t *   ja, DCOMPLEX *   a,
                        DCOMPLEX *  rhs,  pastix_int_t nrhs,   pastix_int_t *   l2g,
                        pastix_int_t  dn, pastix_int_t ** dia, pastix_int_t ** dja, DCOMPLEX ** da,
                        DCOMPLEX ** drhs,  pastix_int_t *  dl2g,
                        MPI_Comm comm, pastix_int_t dof))

#if (defined pastix_float_t)
/*
 *  Function: cscd_save
 *
 *  save a distributed csc to disk.
 *  files are called $(filename) and $(filename)$(RANK)
 *  if filename is NULL then filename = cscd_matrix.
 *
 *  file filename contains the number of processors/files
 *  on first line. Then each line contain the name of each file
 *  (here $(filename)$(RANK)).
 *
 *
 *
 *  Parameters:
 *    n           - Number of local columns
 *    ia          - First cscd starting index of each column in *ja* and *a*
 *    ja          - Row of each element in first CSCD
 *    a           - value of each cscd in first CSCD (can be NULL)
 *    rhs         - Right hand side.
 *    l2g         - local 2 global column numbers for first cscd
 *    dof         - Number of degrees of freedom
 *    filename    - name of the files.
 *    comm        - MPI communicator
 */
int cscd_save(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t * a, pastix_float_t * rhs, pastix_int_t* l2g,
              int dof, const char * filename, MPI_Comm comm);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_save,
                       (pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, float * a, float * rhs, pastix_int_t* l2g,
                        int dof, const char * filename, MPI_Comm comm),
                       (pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, double * a, double * rhs, pastix_int_t* l2g,
                        int dof, const char * filename, MPI_Comm comm),
                       (pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, COMPLEX * a, COMPLEX * rhs,
                        pastix_int_t* l2g,  int dof,  const char * filename, MPI_Comm comm),
                       (pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, DCOMPLEX * a, DCOMPLEX * rhs,
                        pastix_int_t* l2g,  int dof,  const char * filename, MPI_Comm comm))

#if (defined pastix_float_t)
/*
 *  Function: cscd_load
 *
 *  Loads a distributed csc from disk.
 *  if filename is NULL then filename = cscd_matrix.
 *
 *  Parameters:
 *    n           - Number of local columns
 *    ia          - First cscd starting index of each column in *ja* and *a*
 *    ja          - Row of each element in first CSCD
 *    a           - value of each cscd in first CSCD (can be NULL)
 *    rhs         - Right hand side.
 *    l2g         - local 2 global column numbers for first cscd
 *    filename    - name of the files.
 *    comm        - MPI communicator
 */
int cscd_load(pastix_int_t *n, pastix_int_t ** ia, pastix_int_t ** ja, pastix_float_t ** a, pastix_float_t ** rhs, pastix_int_t ** l2g,
              const char * filename, MPI_Comm mpi_comm);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_load,
                       (pastix_int_t *n, pastix_int_t ** ia, pastix_int_t ** ja, float ** a, float ** rhs, pastix_int_t ** l2g,
                        const char * filename, MPI_Comm mpi_comm),
                       (pastix_int_t *n, pastix_int_t ** ia, pastix_int_t ** ja, double ** a, double ** rhs,
                        pastix_int_t ** l2g, const char * filename, MPI_Comm mpi_comm),
                       (pastix_int_t *n, pastix_int_t ** ia, pastix_int_t ** ja, COMPLEX ** a,
                        COMPLEX ** rhs, pastix_int_t ** l2g, const char * filename,
                        MPI_Comm mpi_comm),
                       (pastix_int_t *n, pastix_int_t ** ia, pastix_int_t ** ja, DCOMPLEX ** a,
                        DCOMPLEX ** rhs, pastix_int_t ** l2g, const char * filename,
                        MPI_Comm mpi_comm))

#undef MULTIPLE_TYPE_DEFINE_F
#undef MULTIPLE_TYPE_DEFINE
#undef COMPLEX
#undef DCOMPLEX
#endif /* CSCD_UTILS_H */
