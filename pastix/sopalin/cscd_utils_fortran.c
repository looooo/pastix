/*
  File: pastix_fortran.c

  Interface to the PaStiX API functions.

 */
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif
#include "common.h"
#include "cscd_utils.h"


#if (defined X_ARCHpower_ibm_aix)
#define FORTRAN_CALL(nom) PASTIX_EXTERN_F(nom)
#else
#define FORTRAN_CALL(nom) PASTIX_EXTERN_F(nom ## _)
#endif
/*
  Struct: csc_data_

  Contains the new CSCD

*/
struct csc_data_ {
  pastix_int_t     n;
  pastix_int_t   * colptr;
  pastix_int_t   * rows;
  pastix_float_t * values;
  pastix_float_t * rhs;
  pastix_int_t     nrhs;
  pastix_int_t   * perm;
  pastix_int_t   * l2g;
  pastix_int_t     dof;
};

/*
  Typedef: csc_data_t

  Type coresponding to the struct <csc_data_>
*/
typedef struct csc_data_ csc_data_t;

void FORTRAN_CALL(csc_dispatch_fortran)(csc_data_t ** csc_data,
                                        pastix_int_t         *gN,
                                        pastix_int_t         *gcolptr,
                                        pastix_int_t         *grow,
                                        pastix_float_t       *gavals,
                                        pastix_float_t       *grhs,
                                        pastix_int_t         *gperm,
                                        pastix_int_t         *ginvp,
                                        int         *dispatch,
                                        pastix_int_t         *newn,
                                        pastix_int_t         *newnnz,
                                        MPI_Fint    *fortran_comm)
{
  MPI_Comm        pastix_comm;
  (void)fortran_comm;

  pastix_comm = MPI_Comm_f2c(*fortran_comm);
  MALLOC_INTERN(*csc_data, 1, struct csc_data_);
  (*csc_data)->n      = 0;
  (*csc_data)->colptr = NULL;
  (*csc_data)->rows   = NULL;
  (*csc_data)->values = NULL;
  (*csc_data)->rhs    = NULL;
  (*csc_data)->perm   = NULL;
  (*csc_data)->l2g    = NULL;

  csc_dispatch(*gN, gcolptr, grow, gavals, grhs, gperm, ginvp,
               &((*csc_data)->n), &((*csc_data)->colptr), &((*csc_data)->rows), &((*csc_data)->values),
               &((*csc_data)->rhs), &((*csc_data)->perm),
               &((*csc_data)->l2g), *dispatch, pastix_comm);

  *newn   = (*csc_data)->n;
  *newnnz = (*csc_data)->colptr[(*csc_data)->n]-1;

}

void FORTRAN_CALL(csc_dispatch_fortran_end)(csc_data_t ** csc_data,
                                            pastix_int_t         *lcolptr,
                                            pastix_int_t         *lrow,
                                            pastix_float_t       *lavals,
                                            pastix_float_t       *lrhs,
                                            pastix_float_t       *lperm,
                                            pastix_int_t         *l2g)
{
  pastix_int_t nnz = 0;
  if ((*csc_data)->colptr != NULL) {
    nnz = (*csc_data)->colptr[(*csc_data)->n]-1;
    memcpy(lcolptr, (*csc_data)->colptr, (1+(*csc_data)->n)*sizeof(pastix_int_t));
    free((*csc_data)->colptr);
  }
  if ((*csc_data)->rows != NULL) {
    memcpy(lrow,    (*csc_data)->rows,   nnz*sizeof(pastix_int_t));
    free((*csc_data)->rows);
  }
  if ((*csc_data)->values != NULL) {
    memcpy(lavals,  (*csc_data)->values,   nnz*sizeof(pastix_float_t));
    free((*csc_data)->values);
  }
  if ((*csc_data)->rhs != NULL) {
    memcpy(lrhs, (*csc_data)->rhs, (*csc_data)->n*sizeof(pastix_float_t));
    free((*csc_data)->rhs);
  }
  if ((*csc_data)->perm != NULL) {
    memcpy(lperm, (*csc_data)->perm, (*csc_data)->n*sizeof(pastix_int_t));
    free((*csc_data)->perm);
  }
  if ((*csc_data)->l2g != NULL) {
    memcpy(l2g, (*csc_data)->l2g, (*csc_data)->n*sizeof(pastix_int_t));
    free((*csc_data)->l2g);
  }
}

void FORTRAN_CALL(cscd_redispatch_fortran)(csc_data_t ** csc_data,
                                           pastix_int_t         * n,
                                           pastix_int_t         * ia,
                                           pastix_int_t         * ja,
                                           pastix_float_t       * a,
                                           pastix_float_t       * rhs,
                                           pastix_int_t         * nrhs,
                                           pastix_int_t         * l2g,
                                           pastix_int_t         * newn,
                                           pastix_int_t         * newl2g,
                                           pastix_int_t         * newnnz,
                                           MPI_Fint    * fortran_comm,
                                           pastix_int_t         * dof,
                                           int         * ierr)
{

  MPI_Comm        pastix_comm;
  (void)newnnz; (void)fortran_comm;

  pastix_comm = MPI_Comm_f2c(*fortran_comm);
  MALLOC_INTERN(*csc_data, 1, struct csc_data_);
  (*csc_data)->n      = 0;
  (*csc_data)->colptr = NULL;
  (*csc_data)->rows   = NULL;
  (*csc_data)->values = NULL;
  (*csc_data)->rhs    = NULL;
  (*csc_data)->nrhs   = *nrhs;
  (*csc_data)->perm   = NULL;
  (*csc_data)->l2g    = NULL;
  (*csc_data)->dof    = *dof;

  *ierr = cscd_redispatch(*n,    ia, ja, a, rhs,  *nrhs, l2g,
                          *newn, &((*csc_data)->colptr), &((*csc_data)->rows), &((*csc_data)->values),  &((*csc_data)->rhs), newl2g,
                          pastix_comm, (*csc_data)->dof);
}

void FORTRAN_CALL(cscd_redispatch_fortran_end)(csc_data_t ** csc_data,
                                               pastix_int_t         *dcolptr,
                                               pastix_int_t         *drow,
                                               pastix_float_t       *davals,
                                               pastix_float_t       *drhs)
{
  pastix_int_t nnz = 0;
  if ((*csc_data)->colptr != NULL) {
    nnz = (*csc_data)->colptr[(*csc_data)->n]-1;
    memcpy(dcolptr, (*csc_data)->colptr, (1+(*csc_data)->n)*sizeof(pastix_int_t));
    free((*csc_data)->colptr);
  }
  if ((*csc_data)->rows != NULL) {
    memcpy(drow,    (*csc_data)->rows,   nnz*sizeof(pastix_int_t));
    free((*csc_data)->rows);
  }
  if ((*csc_data)->values != NULL) {
    memcpy(davals,  (*csc_data)->values,
           (*csc_data)->dof*(*csc_data)->dof*nnz*sizeof(pastix_float_t));
    free((*csc_data)->values);
  }
  if ((*csc_data)->rhs != NULL) {
    memcpy(drhs, (*csc_data)->rhs, (*csc_data)->nrhs*(*csc_data)->n*sizeof(pastix_float_t));
    free((*csc_data)->rhs);
  }
  memFree_null(*csc_data);
}
