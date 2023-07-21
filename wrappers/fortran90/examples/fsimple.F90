!
! @file fsimple.F90
!
! Fortran 90 example using a matrix read with the spm driver.
!
! @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.3.0
! @author Mathieu Faverge
! @date 2023-01-17
!
program fsimple
  use iso_c_binding
  use spmf
  use pastixf

  implicit none

  type(pastix_data_t),        pointer                      :: pastix_data
  type(spmatrix_t),           pointer                      :: spm
  type(spmatrix_t),           pointer                      :: spm2
  integer(kind=pastix_int_t), target                       :: iparm(iparm_size)
  real(kind=c_double),        target                       :: dparm(dparm_size)
  real(kind=c_double)                                      :: normA
  integer(c_int)                                           :: info
  integer(kind=pastix_int_t)                               :: nrhs
  real(kind=c_double), dimension(:,:), allocatable, target :: x0, x, b

  !
  ! Solve the problem
  !

  ! 1- Initialize the parameters and the solver
  call pastixInitParam( iparm, dparm )
  call pastixInit( pastix_data, MPI_COMM_WORLD, iparm, dparm )

  !
  ! Initialize the problem
  !   1- The matrix
  allocate( spm )
  call spmReadDriver( SpmDriverLaplacian, "d:10:10:10:2.", spm, info )

  allocate( spm2 )
  call spmCheckAndCorrect( spm, spm2, info )
  if ( info .ne. 0 ) then
     call spmExit( spm )
     spm = spm2
  end if
  deallocate( spm2 )

  call spmPrintInfo( spm )

  ! Scale A for better stability with low-rank computations
  call spmNorm( SpmFrobeniusNorm, spm, normA )
  call spmScal( 1. / normA, spm )

  !   2- The right hand side
  nrhs = 10
  allocate(x0(spm%nexp, nrhs))
  allocate(x( spm%nexp, nrhs))
  allocate(b( spm%nexp, nrhs))

  call spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm%nexp, b, spm%nexp, info )
  x = b

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

end program fsimple
