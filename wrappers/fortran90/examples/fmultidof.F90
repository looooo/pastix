!
! @file fmultidof.F90
!
! Fortran 90 example using a matrix read with the spm driver.
!
! @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.4.0
! @author Mathieu Faverge
! @author Alycia Lisito
! @date 2024-07-05
!
program fmultidof
  use iso_c_binding
  use spmf
  use pastixf

  implicit none

  type(pastix_data_t),        pointer                      :: pastix_data
  type(pastix_order_t),       pointer                      :: order => null()
  type(spmatrix_t),           pointer                      :: spm
  !type(spmatrix_t),           pointer                      :: spm2
  integer(kind=pastix_int_t), target                       :: iparm(iparm_size)
  real(kind=c_double),        target                       :: dparm(dparm_size)
  real(kind=c_double)                                      :: normA
  integer(c_int)                                           :: info
  integer(kind=pastix_int_t)                               :: nrhs
  real(kind=c_double), dimension(:,:), allocatable, target :: x0, x, b

  ! 1- Initialize the parameters and the solver
  call pastixInitParam( iparm, dparm )
  call pastixInit( pastix_data, MPI_COMM_WORLD, iparm, dparm )

  !
  ! Initialize the problem
  !   1- The matrix
  ! See genLaplacian function for the parameters
  ! Here we generate a 10x101x0 laplacian with 5 degrees of freedom per unknown
  !
  allocate( spm )
  call spmReadDriver( SpmDriverLaplacian, "d:10:10:10:4.:.5:5", spm, info )
  if ( info .ne. 0 ) then
     error stop "spmReadDriver failed"
  end if

  !
  ! We do not check the matrix as it would expand the spm matrix, and
  ! the compressed version would not be used. If you reuse this
  ! testing, you have to make sure that the compressed form of the
  ! matrix is correct without calling the spmCheckAndCorrect. One way
  ! to do it, is to call the function on the matrix without the values
  ! to make sure its correct.
  !
  ! allocate( spm2 )
  ! call spmCheckAndCorrect( spm, spm2, info )
  ! if ( info .ne. 0 ) then
  !    call spmExit( spm )
  !    spm = spm2
  ! end if
  ! deallocate( spm2 )

  ! Print information about the matrix
  call spmPrintInfo( spm )

  ! Scale A for better stability with low-rank computations
  call spmNorm( SpmFrobeniusNorm, spm, normA )
  call spmScal( 1. / normA, spm )

  !
  ! Solve the problem
  !

  ! 2- Analyze the problem
  !   a- Split the analyze in substeps. The first three steps handle
  !      multi-dof spm, the last one does not yet
  call pastix_subtask_order( pastix_data, spm, order, info );
  if ( info .ne. 0 ) then
     error stop "pastix_subtask_order failed"
  end if

  call pastix_subtask_symbfact( pastix_data, info );
  if ( info .ne. 0 ) then
     error stop "pastix_subtask_symbfact failed"
  end if

  call pastix_subtask_reordering( pastix_data, info );
  if ( info .ne. 0 ) then
     error stop "pastix_subtask_reordering failed"
  end if

  !   b- Expand the matrix and the associated substructure
  call pastixExpand( pastix_data, spm )

  !   c- Finish the analyze step
  call pastix_subtask_blend( pastix_data, info )
  if ( info .ne. 0 ) then
     error stop "pastix_subtask_blend failed"
  end if

  ! 3- Factorize the matrix
  call pastix_task_numfact( pastix_data, spm, info )
  if ( info .ne. 0 ) then
     error stop "pastix_task_numfact failed"
  end if

  !
  ! We need to generate the right hand side once the spm has been
  ! expanded, spmGenRHS is not yet compatible with multi-dof
  !
  nrhs = 10
  allocate(x0(spm%nexp, nrhs))
  allocate(x( spm%nexp, nrhs))
  allocate(b( spm%nexp, nrhs))

  call spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm%nexp, b, spm%nexp, info )
  if ( info .ne. 0 ) then
     error stop "spmGenRHS failed"
  end if
  x = b

  ! 4- Solve the problem
  call pastix_task_solve( pastix_data, spm%nexp, nrhs, x, spm%nexp, info )
  if ( info .ne. 0 ) then
     error stop "pastix_task_solve failed"
  end if

  ! 5- Refine the solution
  call pastix_task_refine( pastix_data, spm%nexp, nrhs, b, spm%nexp, x, spm%nexp, info )
  if ( info .ne. 0 ) then
     error stop "pastix_task_refine failed"
  end if

  ! 6- Destroy the C data structure
  call pastixFinalize( pastix_data )

  !
  ! Check the solution
  !
  call spmCheckAxb( dparm(DPARM_EPSILON_REFINEMENT), nrhs, spm, x0, spm%nexp, b, spm%nexp, x, spm%nexp, info )
  if ( info .ne. 0 ) then
     error stop "spmCheckAxb failed"
  end if

  call spmExit( spm )
  deallocate(spm)
  deallocate(x0)
  deallocate(x)
  deallocate(b)

end program fmultidof
