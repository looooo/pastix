/**
 *  PaStiX CSC management routines.
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
  File: pastix_fortran.c

  Interface to the PaStiX API functions.

 */
#include "common.h"
#include "z_cscd_utils.h"


#define FORTRAN_NAME(nu,nl,pl,pc)   \
  void nl pl;                       \
  void nu pl                        \
  { nl pc; }                        \
  void nl ## _ pl                   \
  { nl pc; }                        \
  void nl ## __ pl                  \
  { nl pc; }
/*
  Struct: csc_data_

  Contains the new CSCD

*/
struct csc_data_ {
  pastix_int_t        n;
  pastix_int_t       *colptr;
  pastix_int_t       *rows;
  pastix_complex64_t *values;
  pastix_complex64_t *rhs;
  pastix_int_t        nrhs;
  pastix_int_t       *perm;
  pastix_int_t       *l2g;
  pastix_int_t        dof;
};

/*
  Typedef: csc_data_t

  Type coresponding to the struct <csc_data_>
*/
typedef struct csc_data_ csc_data_t;

FORTRAN_NAME(Z_CSC_DISPATCH_FORTRAN,
             z_csc_dispatch_fortran,
             (csc_data_t         **csc_data,
              pastix_int_t        *gN,
              pastix_int_t        *gcolptr,
              pastix_int_t        *grow,
              pastix_complex64_t  *gavals,
              pastix_complex64_t  *grhs,
              pastix_int_t        *gperm,
              pastix_int_t        *ginvp,
              int                 *dispatch,
              pastix_int_t        *newn,
              pastix_int_t        *newnnz,
              MPI_Fint            *fortran_comm),
             (csc_data,
              gN,
              gcolptr,
              grow,
              gavals,
              grhs,
              gperm,
              ginvp,
              dispatch,
              newn,
              newnnz,
              fortran_comm))

void z_csc_dispatch_fortran(csc_data_t         **csc_data,
                            pastix_int_t        *gN,
                            pastix_int_t        *gcolptr,
                            pastix_int_t        *grow,
                            pastix_complex64_t  *gavals,
                            pastix_complex64_t  *grhs,
                            pastix_int_t        *gperm,
                            pastix_int_t        *ginvp,
                            int                 *dispatch,
                            pastix_int_t        *newn,
                            pastix_int_t        *newnnz,
                            MPI_Fint            *fortran_comm)
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

  z_csc_dispatch(*gN, gcolptr, grow, gavals, grhs, gperm, ginvp,
                 &((*csc_data)->n), &((*csc_data)->colptr), &((*csc_data)->rows), &((*csc_data)->values),
               &((*csc_data)->rhs), &((*csc_data)->perm),
               &((*csc_data)->l2g), *dispatch, pastix_comm);

  *newn   = (*csc_data)->n;
  *newnnz = (*csc_data)->colptr[(*csc_data)->n]-1;

}



FORTRAN_NAME(Z_CSC_DISPATCH_FORTRAN_END,
             z_csc_dispatch_fortran_end,
             (csc_data_t         **csc_data,
              pastix_int_t        *lcolptr,
              pastix_int_t        *lrow,
              pastix_complex64_t  *lavals,
              pastix_complex64_t  *lrhs,
              pastix_int_t        *lperm,
              pastix_int_t        *l2g),
             (csc_data,
              lcolptr,
              lrow,
              lavals,
              lrhs,
              lperm,
              l2g))
void z_csc_dispatch_fortran_end(csc_data_t         **csc_data,
                                pastix_int_t        *lcolptr,
                                pastix_int_t        *lrow,
                                pastix_complex64_t  *lavals,
                                pastix_complex64_t  *lrhs,
                                pastix_int_t       *lperm,
                                pastix_int_t        *l2g)
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
    memcpy(lavals,  (*csc_data)->values,   nnz*sizeof(pastix_complex64_t));
    free((*csc_data)->values);
  }
  if ((*csc_data)->rhs != NULL) {
    memcpy(lrhs, (*csc_data)->rhs, (*csc_data)->n*sizeof(pastix_complex64_t));
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

FORTRAN_NAME(Z_CSCD_REDISPATCH_FORTRAN,
             z_cscd_redispatch_fortran,
             (csc_data_t         **csc_data,
              pastix_int_t        *n,
              pastix_int_t        *ia,
              pastix_int_t        *ja,
              pastix_complex64_t  *a,
              pastix_complex64_t  *rhs,
              pastix_int_t        *nrhs,
              pastix_int_t        *l2g,
              pastix_int_t        *newn,
              pastix_int_t        *newl2g,
              pastix_int_t        *newnnz,
              MPI_Fint            *fortran_comm,
              pastix_int_t        *dof,
              int                 *ierr),
             (csc_data,
              n,
              ia,
              ja,
              a,
              rhs,
              nrhs,
              l2g,
              newn,
              newl2g,
              newnnz,
              fortran_comm,
              dof,
              ierr))
void z_cscd_redispatch_fortran(csc_data_t         **csc_data,
                               pastix_int_t        *n,
                               pastix_int_t        *ia,
                               pastix_int_t        *ja,
                               pastix_complex64_t  *a,
                               pastix_complex64_t  *rhs,
                               pastix_int_t        *nrhs,
                               pastix_int_t        *l2g,
                               pastix_int_t        *newn,
                               pastix_int_t        *newl2g,
                               pastix_int_t        *newnnz,
                               MPI_Fint            *fortran_comm,
                               pastix_int_t        *dof,
                               int                 *ierr)
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

  *ierr = z_cscd_redispatch(*n,    ia, ja, a, rhs,  *nrhs, l2g,
                            *newn, &((*csc_data)->colptr), &((*csc_data)->rows), &((*csc_data)->values),  &((*csc_data)->rhs), newl2g,
                          pastix_comm, (*csc_data)->dof);
}

FORTRAN_NAME(Z_CSCD_REDISPATCH_FORTRAN_END,
             z_cscd_redispatch_fortran_end,
             (csc_data_t         **csc_data,
              pastix_int_t        *dcolptr,
              pastix_int_t        *drow,
              pastix_complex64_t  *davals,
              pastix_complex64_t  *drhs),
             (csc_data,
              dcolptr,
              drow,
              davals,
              drhs))


void z_cscd_redispatch_fortran_end(csc_data_t         **csc_data,
                                   pastix_int_t        *dcolptr,
                                   pastix_int_t        *drow,
                                   pastix_complex64_t  *davals,
                                   pastix_complex64_t  *drhs)    
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
           (*csc_data)->dof*(*csc_data)->dof*nnz*sizeof(pastix_complex64_t));
    free((*csc_data)->values);
  }
  if ((*csc_data)->rhs != NULL) {
    memcpy(drhs, (*csc_data)->rhs, (*csc_data)->nrhs*(*csc_data)->n*sizeof(pastix_complex64_t));
    free((*csc_data)->rhs);
  }
  memFree_null(*csc_data);
}
