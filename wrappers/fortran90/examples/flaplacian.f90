!
! @file flaplacian.f90
!
! Fortran 90 example using a laplacian matrix.
!
! @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
program fsimple
  use iso_c_binding
  use pastix_enums
  use spmf
  use pastixf
  implicit none

  integer(kind=pastix_int_t),     dimension(:), allocatable, target   :: rowptr
  integer(kind=pastix_int_t),     dimension(:), allocatable, target   :: colptr
  complex(kind=c_double_complex), dimension(:), allocatable, target   :: values
  complex(kind=c_double_complex), dimension(:,:), allocatable, target :: x0, x, b
  type(c_ptr)                                                         :: x0_ptr, x_ptr, b_ptr
  type(pastix_data_t),        pointer                                 :: pastix_data
  type(pastix_spm_t),         target                                  :: spm
  type(pastix_spm_t),         pointer                                 :: spm2
  integer(kind=pastix_int_t), target                                  :: iparm(iparm_size)
  real(kind=c_double),        target                                  :: dparm(dparm_size)
  integer(kind=pastix_int_t)                                          :: dim1, dim2, dim3, n, nnz
  integer(kind=pastix_int_t)                                          :: i, j, k, l, nrhs
  integer(c_int)                                                      :: info

  !
  ! Generate a 10x10x10 complex Laplacian
  !
  dim1 = 10
  dim2 = 10
  dim3 = 10
  n    = dim1 * dim2 * dim3
  nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)

  allocate(rowptr(nnz))
  allocate(colptr(nnz))
  allocate(values(nnz))

  l = 1
  do i=1,dim1
     do j=1,dim2
        do k=1,dim3
           rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           values(l) = 6.

           if (i == 1) then
              values(l) = values(l) - 1.
           end if
           if (i == dim1) then
              values(l) = values(l) - 1.
           end if
           if (j == 1) then
              values(l) = values(l) - 1.
           end if
           if (j == dim2) then
              values(l) = values(l) - 1.
           end if
           if (k == 1) then
              values(l) = values(l) - 1.
           end if
           if (k == dim3) then
              values(l) = values(l) - 1.
           end if

           values(l) = values(l) * 8.
           l = l + 1

           if (i < dim1) then
              rowptr(l) =  i    + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = - 1. - 1. * I
              l = l + 1
           end if
           if (j < dim2) then
              rowptr(l) = (i-1) + dim1 *  j    + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = - 1. - 1. * I
              l = l + 1
           end if
           if (k < dim3) then
              rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 *  k    + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = -1. - 1. * I
              l = l + 1
           end if
        end do
     end do
  end do

  if ( l .ne. nnz+1 ) then
     write(6,*) 'l ', l, " nnz ", nnz
  end if

  !
  ! Create the spm out of the internal data
  !
  call spmInit( spm )
  spm%mtxtype = PastixHermitian
  spm%flttype = PastixComplex64
  spm%fmttype = PastixIJV
  spm%n       = n
  spm%nnz     = nnz
  spm%dof     = 1
  spm%rowptr  = c_loc(rowptr)
  spm%colptr  = c_loc(colptr)
  spm%values  = c_loc(values)

  call spmUpdateComputedFields( spm )

  call spmCheckAndCorrect( spm, spm2 )
  if (.not. c_associated(c_loc(spm), c_loc(spm2))) then
     deallocate(rowptr)
     deallocate(colptr)
     deallocate(values)

     spm%rowptr = c_null_ptr
     spm%colptr = c_null_ptr
     spm%values = c_null_ptr

     call spmExit( spm )
     spm = spm2
  end if

  call spmPrintInfo( spm )

  !   2- The right hand side
  nrhs = 10
  allocate(x0(spm%n, nrhs))
  allocate(x( spm%n, nrhs))
  allocate(b( spm%n, nrhs))
  x0_ptr = c_loc(x0)
  x_ptr  = c_loc(x)
  b_ptr  = c_loc(b)

  call spmGenRHS( PastixRhsRndX, nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, info )
  x = b

  !
  ! Solve the problem
  !

  ! 1- Initialize the parameters and the solver
  call pastixInitParam( iparm, dparm )
  call pastixInit( pastix_data, 0, iparm, dparm )

  ! 2- Analyze the problem
  call pastix_task_analyze( pastix_data, spm, info )

  ! 3- Factorize the matrix
  call pastix_task_numfact( pastix_data, spm, info )

  ! 4- Solve the problem
  call pastix_task_solve( pastix_data, nrhs, x_ptr, spm%n, info )

  ! 5- Refine the solution
  call pastix_task_refine( pastix_data, spm%n, nrhs, b_ptr, spm%n, x_ptr, spm%n, info )

  ! 6- Destroy the C data structure
  call pastixFinalize( pastix_data )

  !
  ! Check the solution
  !
  call spmCheckAxb( nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, x_ptr, spm%n, info )

  call spmExit( spm )
  deallocate(x0)
  deallocate(x)
  deallocate(b)

end program fsimple
