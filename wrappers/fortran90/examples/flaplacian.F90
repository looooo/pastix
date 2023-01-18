!
! @file flaplacian.F90
!
! Fortran 90 example using a laplacian matrix.
!
! @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.3.0
! @author Mathieu Faverge
! @date 2022-10-11
!
program flaplacian
  use iso_c_binding
  use spmf
  use pastixf
  implicit none

  integer(kind=spm_int_t),        dimension(:),   pointer             :: rowptr
  integer(kind=spm_int_t),        dimension(:),   pointer             :: colptr
  complex(kind=c_double_complex), dimension(:),   pointer             :: values
  complex(kind=c_double_complex), dimension(:,:), allocatable, target :: x0, x, b
  type(pastix_data_t),        pointer                                 :: pastix_data
  type(spmatrix_t),           pointer                                 :: spm
  type(spmatrix_t),           pointer                                 :: spm2
  integer(kind=pastix_int_t), target                                  :: iparm(iparm_size)
  real(kind=c_double),        target                                  :: dparm(dparm_size)
  integer(kind=pastix_int_t)                                          :: dim1, dim2, dim3, n, nnz
  integer(kind=pastix_int_t)                                          :: i, j, k, l, nrhs
  integer(c_int)                                                      :: info

  ! 1- Initialize the parameters and the solver
  ! (Done before any calls to spm to automatically intialize MPI if needed)
  call pastixInitParam( iparm, dparm )
  call pastixInit( pastix_data, MPI_COMM_WORLD, iparm, dparm )

  !
  ! Generate a 10x10x10 complex Laplacian
  !
  dim1 = 10
  dim2 = 10
  dim3 = 10
  n    = dim1 * dim2 * dim3
  nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)

  !
  ! Create the spm out of the internal data
  !
  allocate( spm )
  call spmInit( spm )
  spm%baseval = 1
  spm%mtxtype = SpmSymmetric
  spm%flttype = SpmComplex64
  spm%fmttype = SpmIJV
  spm%n       = n
  spm%nnz     = nnz
  spm%dof     = 1

  call spmUpdateComputedFields( spm )
  call spmAlloc( spm )
  call spmGetArray( spm, colptr=colptr, rowptr=rowptr, zvalues=values )

  l = 1
  do i=1,dim1
     do j=1,dim2
        do k=1,dim3
           rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           values(l) = (6., 0.)

           if (i == 1) then
              values(l) = values(l) - (1., 0.)
           end if
           if (i == dim1) then
              values(l) = values(l) - (1., 0.)
           end if
           if (j == 1) then
              values(l) = values(l) - (1., 0.)
           end if
           if (j == dim2) then
              values(l) = values(l) - (1., 0.)
           end if
           if (k == 1) then
              values(l) = values(l) - (1., 0.)
           end if
           if (k == dim3) then
              values(l) = values(l) - (1., 0.)
           end if

           values(l) = values(l) * 8.
           l = l + 1

           if (i < dim1) then
              rowptr(l) =  i    + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = (- 1.,  - 1.)
              l = l + 1
           end if
           if (j < dim2) then
              rowptr(l) = (i-1) + dim1 *  j    + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = (- 1., - 1.)
              l = l + 1
           end if
           if (k < dim3) then
              rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 *  k    + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = (-1., - 1. )
              l = l + 1
           end if
        end do
     end do
  end do

  if ( l .ne. nnz+1 ) then
     write(6,*) 'l ', l, " nnz ", nnz
  end if

  allocate( spm2 )
  call spmCheckAndCorrect( spm, spm2, info )
  if ( info .ne. 0 ) then
     call spmExit( spm )
     spm = spm2
  end if
  deallocate( spm2 )

  call spmPrintInfo( spm )

  !   2- The right hand side
  nrhs = 10
  allocate(x0(spm%nexp, nrhs))
  allocate(x( spm%nexp, nrhs))
  allocate(b( spm%nexp, nrhs))

  call spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm%nexp, b, spm%nexp, info )
  x = b

  !
  ! Solve the problem
  !

  ! 2- Analyze the problem
  call pastix_task_analyze( pastix_data, spm, info )

  ! 3- Factorize the matrix
  call pastix_task_numfact( pastix_data, spm, info )

  ! 4- Solve the problem
  call pastix_task_solve( pastix_data, spm%nexp, nrhs, x, spm%nexp, info )

  ! 5- Refine the solution
  call pastix_task_refine( pastix_data, spm%nexp, nrhs, b, spm%nexp, x, spm%nexp, info )

  ! 6- Destroy the C data structure
  call pastixFinalize( pastix_data )

  !
  ! Check the solution
  !
  call spmCheckAxb( dparm(DPARM_EPSILON_REFINEMENT), nrhs, spm, x0, spm%nexp, b, spm%nexp, x, spm%nexp, info )

  call spmExit( spm )
  deallocate(spm)
  deallocate(x0)
  deallocate(x)
  deallocate(b)

end program flaplacian
