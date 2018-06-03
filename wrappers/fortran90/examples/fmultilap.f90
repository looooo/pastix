!
! @file flaplacian.f90
!
! Fortran 90 multiple rhs example.
!
! This example tries to solve the multi-rhs problem in three different manners:
!    1) The matrix is duplicated and solve as many times as the number
!        of RHS, both by the same subset of threads.
!    2) The matrix is duplicated per thread, solved in sequential
!       by each thread, and each thread solves a subset of RHS.
!    3) A single matrix is solved with all the available threads and
!       is used to solve all the RHS.
!
! @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
program flaplacian
  use iso_c_binding
  use pastix_enums
  use spmf
  use pastixf
  !$ use omp_lib
  implicit none

  type multilap_param
     logical :: checkmat    ! Enable/disable the check and correct on the spm matrix
     logical :: multirhs    ! Enable/disable multi-rhs solves
     logical :: output      ! Enable output information
     integer :: dim1        ! First dimension of the Laplacian matrix
     integer :: dim2        ! Second dimension of the Laplacian matrix
     integer :: dim3        ! Third dimension of the Laplacian matrix
     integer :: n           ! Size of the matrix
     integer :: nnz         ! Number of non zeroes in the matrix
     integer :: nb_outit    ! Number of outer iteration
     integer :: nb_mat      ! Number of matrices to factorize per iteration
     integer :: nb_rhs      ! Number of right hand side per matrix
     integer :: nb_solve    ! Number of solve steps per RHS
     integer :: nb_sys      ! Number of systems
     integer :: nb_thrd     ! Number of thread available
     integer :: nb_fact_omp ! Number of PaStiX instances for the factorization
     integer :: nb_fact_thd ! Number of threads per PaStiX instances for the factorization
     integer :: nb_solv_omp ! Number of PaStiX instances for the solve
     integer :: nb_solv_thd ! Number of threads per PaStiX instances for the solve
  end type multilap_param

  type rhs_arr
     integer(kind=pastix_int_t)                            :: idmat
     integer(kind=pastix_int_t)                            :: nrhs
     integer(kind=pastix_int_t)                            :: n
     integer(kind=pastix_int_t)                            :: size
     integer(kind=pastix_int_t)                            :: ori
     complex(kind=c_double_complex), dimension(:), pointer :: b0
     complex(kind=c_double_complex), dimension(:), pointer :: b
     complex(kind=c_double_complex), dimension(:), pointer :: x
  end type rhs_arr

  type sys_lin
     integer                                               :: nsys, pid
     type(pastix_data_t),                          pointer :: pastix_data
     integer(kind=pastix_int_t),     dimension(iparm_size) :: iparm
     real(kind=c_double),            dimension(dparm_size) :: dparm
     integer(kind=c_int),            dimension(:), pointer :: bindtab
     type(spmatrix_t),                             pointer :: spm
     type(spmatrix_t),                             pointer :: spm2
     integer(kind=pastix_int_t),     dimension(:), pointer :: rowptr
     integer(kind=pastix_int_t),     dimension(:), pointer :: colptr
     complex(kind=c_double_complex), dimension(:), pointer :: values
     type(rhs_arr),                  dimension(:), pointer :: rhs
  end type sys_lin

  type(sys_lin), dimension(:), allocatable, target  :: sla_lap
  type(sys_lin),                            pointer :: matrix
  type(rhs_arr),                            pointer :: rhs

  type(c_ptr)                                                       :: b0_ptr, x_ptr, b_ptr
  integer(kind=pastix_int_t),     dimension(iparm_size)             :: iparm
  real(kind=c_double),            dimension(dparm_size)             :: dparm
  integer                       , dimension(:), allocatable         :: ila_thrmn, ila_thrmx
  integer                       , dimension(:), allocatable         :: ila_thrsz
  integer(kind=pastix_int_t)                                        :: nrhs
  integer                                                           :: th, im, ir, i, j, k
  integer(c_int)                                                    :: info, ginfo = 0
  !
  integer, parameter :: MULTILAP_ANALYZE_TIME = 1
  integer, parameter :: MULTILAP_FACT_TIME    = 2
  integer, parameter :: MULTILAP_FACT_FLOPS   = 3
  integer, parameter :: MULTILAP_SOLV_TIME    = 4
  integer, parameter :: MULTILAP_MAX_STAT     = MULTILAP_SOLV_TIME
  !
  real(kind=c_double),  allocatable, dimension(:,:) :: dla_thread_stats
  real(kind=c_double),  dimension(MULTILAP_MAX_STAT)  :: dla_final_stats = 0.0

  type(multilap_param) :: params
  !
  ! Get the problem confirguration
  !
  call multilap_init( params )

  !
  ! Set some parameters in iparm
  !
  call pastixInitParam( iparm, dparm )

  iparm(IPARM_THREAD_NBR) = params%nb_fact_thd
  if ( iparm(IPARM_THREAD_NBR) == 1 ) then
     iparm(IPARM_SCHEDULER)  = PastixSchedSequential
  else
     iparm(IPARM_SCHEDULER)  = PastixSchedStatic
  end if
  iparm(IPARM_VERBOSE)       = PastixVerboseNot
  iparm(IPARM_FACTORIZATION) = PastixFactLLT

  ! Initialize the matrices structures
  allocate(sla_lap(params%nb_mat))
  do im = 1, params%nb_mat
     !
     matrix => sla_lap(im)
     !
     matrix%iparm(:) = iparm(:)
     matrix%dparm(:) = dparm(:)
     !

     ! Initialize the solve structures
     if ( params%multirhs ) then
        !
        ! multi-thread RHS is enabled
        ! Distribute the rhs solve in the same manner as the matrices factorizations
        !
        matrix%nsys = 1
        allocate(matrix%rhs( matrix%nsys ))

        rhs => matrix%rhs(1)
        rhs%nrhs  = params%nb_rhs
        rhs%n     = params%n
        rhs%size  = rhs%nrhs * rhs%n
        rhs%ori   = 0

        if ( params%output ) then
           ! Allocate a backup of b that will be destroyed by the check
           allocate(rhs%b0(rhs%size))
        end if
        allocate(rhs%b( rhs%size))
        allocate(rhs%x( rhs%size))

     else
        !
        ! multi-thread rhs is disabled
        ! Distribute the rhs solve among each sequential thread issued from the factorization
        !
        call split_parall( params%nb_rhs, params%nb_fact_thd )

        matrix%nsys = iparm(IPARM_THREAD_NBR)
        allocate(matrix%rhs( matrix%nsys ))

        do th = 1, iparm(IPARM_THREAD_NBR)

           rhs => matrix%rhs(th)

           rhs%nrhs = ila_thrsz(th)
           rhs%n    = params%n
           rhs%size = rhs%nrhs * rhs%n
           rhs%ori  = (ila_thrmn(th)-1) * rhs%n

           if ( rhs%nrhs .gt. 0 ) then
              if ( params%output ) then
                 ! Allocate a backup of b that will be destroyed by the check
                 allocate(rhs%b0(rhs%size))
              end if
              allocate(rhs%b( rhs%size))
              allocate(rhs%x( rhs%size))
           end if
        end do
     end if
  end do

  if (params%output) then
     write(6,*) '!--------------------------------------------------------------------!'
     write(6,*) '        Size of x = ', params%n
     write(6,*) '        Matrix      NRHS      Start         End'

     do im = 1, params%nb_mat
        matrix => sla_lap(im)
        do j = 1, matrix%nsys
           rhs => matrix%rhs(j)
           write(6,*) im, rhs%nrhs, rhs%ori + 1, rhs%ori + rhs%size
        end do
     end do
  end if

  !
  ! Outmost iteration
  !
  main_loop: do k = 1, params%nb_outit
     !
     write(6,*) '!====================================================================!'
     write(6,*) ' Outer iteration: ', k
     write(6,*) '!--------------------------------------------------------------------!'
     write(6,*) ' Nb of factorization performed in parallel      = ', params%nb_fact_omp
     write(6,*) ' Nb of threads used by PaStiX per factorization = ', params%nb_fact_thd

     !
     ! Distribute factorizations among parallel threads
     !
     call split_parall( params%nb_mat, params%nb_fact_omp )
     !
     allocate(dla_thread_stats(params%nb_mat, MULTILAP_MAX_STAT))
     dla_thread_stats(:,:) = 0.0

     !
     !$OMP PARALLEL NUM_THREADS(params%nb_fact_omp) DEFAULT(NONE) &
     !$OMP SHARED(params, sla_lap, k)                 &
     !$OMP SHARED(ila_thrsz, ila_thrmn, ila_thrmx)    &
     !$OMP SHARED(dla_thread_stats, iparm, dparm)     &
     !$OMP PRIVATE(th, im, ir, i, j )                 &
     !$OMP PRIVATE(matrix, rhs, x_ptr, b_ptr, b0_ptr) &
     !$OMP PRIVATE(info)

     !$OMP DO SCHEDULE(STATIC,1)
     fact_loop: do th = 1, params%nb_fact_omp
        do im = ila_thrmn(th), ila_thrmx(th)

           matrix => sla_lap(im)

           matrix%pid = th
           !
           ! Give each matrix a subset of cores
           !
           allocate( matrix%bindtab( iparm(IPARM_THREAD_NBR) ) )
           do i = 1, iparm(IPARM_THREAD_NBR)
              matrix%bindtab( i ) = (th-1) * iparm(IPARM_THREAD_NBR) + (i - 1)
           end do

           ! 0- Generate the SPM matrix
           call multilap_genOneMatrix( params, matrix, k, im )

           ! 1- Initialize the parameters and the solver
           call pastixInitWithAffinity( matrix%pastix_data, 0, &
                & matrix%iparm, matrix%dparm, &
                & matrix%bindtab )

           ! 2- Analyze the problem
           call pastix_task_analyze( matrix%pastix_data, matrix%spm, info )
           dla_thread_stats(th,MULTILAP_ANALYZE_TIME) = &
                & dla_thread_stats(th,MULTILAP_ANALYZE_TIME) + matrix%dparm(DPARM_ANALYZE_TIME)

           ! 3- Factorize the matrix
           call pastix_task_numfact( matrix%pastix_data, matrix%spm, info )

           dla_thread_stats(th,MULTILAP_FACT_TIME) = &
                & dla_thread_stats(th,MULTILAP_FACT_TIME) + matrix%dparm(DPARM_FACT_TIME)
           dla_thread_stats(th,MULTILAP_FACT_FLOPS) = &
                & dla_thread_stats(th,MULTILAP_FACT_FLOPS) + (matrix%dparm(DPARM_FACT_FLOPS)/(1024.0**3))

        end do
     end do fact_loop
     !$OMP END DO
     !$OMP END PARALLEL

     if ( params%output ) then
        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) '! Results per matrix'
        write(6,*) '!'
        do im = 1, params%nb_mat
           matrix => sla_lap(im)

           write(6,*) ' Matrix ', im, ' done by instance ', matrix%pid
           write(6,*) ' Time for analysis      ', matrix%dparm(DPARM_ANALYZE_TIME)
           write(6,*) ' Pred Time for fact     ', matrix%dparm(DPARM_PRED_FACT_TIME)
           write(6,*) ' Time for factorization ', matrix%dparm(DPARM_FACT_TIME)
           write(6,*) ' GFlops/s for fact      ', matrix%dparm(DPARM_FACT_FLOPS)/(1024.0**3)
        end do

        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) '! Results per PaStiX instance'
        write(6,*) '!'
        do th = 1, params%nb_fact_omp
           dla_thread_stats(th,MULTILAP_FACT_FLOPS) = &
                & dla_thread_stats(th,MULTILAP_FACT_FLOPS) / ila_thrsz(th)
           write(6,*) ' Thread ',th
           write(6,*) ' Time for analysis      ', dla_thread_stats(th,MULTILAP_ANALYZE_TIME)
           write(6,*) ' Time for factorization ', dla_thread_stats(th,MULTILAP_FACT_TIME)
           write(6,*) ' GFlops/s for fact      ', dla_thread_stats(th,MULTILAP_FACT_FLOPS)
        end do

        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) ' Outer iterate ', k
        write(6,*) ' Time for analysis      ', maxval(dla_thread_stats(:,MULTILAP_ANALYZE_TIME))
        write(6,*) ' Time for factorization ', maxval(dla_thread_stats(:,MULTILAP_FACT_TIME))
        write(6,*) ' GFlops/s for fact      ', sum(dla_thread_stats(:,MULTILAP_FACT_FLOPS)) / params%nb_fact_omp
     end if

     dla_final_stats(MULTILAP_ANALYZE_TIME) = dla_final_stats(MULTILAP_ANALYZE_TIME) + &
          & maxval(dla_thread_stats(:,MULTILAP_ANALYZE_TIME))
     dla_final_stats(MULTILAP_FACT_TIME) = dla_final_stats(MULTILAP_FACT_TIME) + &
          & maxval(dla_thread_stats(:,MULTILAP_FACT_TIME))
     dla_final_stats(MULTILAP_FACT_FLOPS) = dla_final_stats(MULTILAP_FACT_FLOPS) + &
          & (sum(dla_thread_stats(:,MULTILAP_FACT_FLOPS)) / params%nb_fact_omp)
     !
     deallocate(dla_thread_stats)

     !
     ! Distribute solves among parallel threads
     !
     allocate(dla_thread_stats( params%nb_solv_omp, MULTILAP_MAX_STAT ))
     dla_thread_stats(:,:) = 0.0

     write(6,*) '!--------------------------------------------------------------------!'
     write(6,*) ' Nb of OpenMP threads enrolled for solution      = ', params%nb_solv_omp
     write(6,*) ' Nb of pthreads used in PaStiX for solution      = ', params%nb_solv_thd

     call split_parall( params%nb_sys, params%nb_solv_omp )

     do i = 1, params%nb_solve

        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) ' Solve iteration nr. ', i

        ! Apply the forward permutation on b and x
        do im = 1, params%nb_mat
           matrix => sla_lap(im)
           do j = 1, matrix%nsys
              rhs => matrix%rhs(j)

              if ( .not. (rhs%nrhs .gt. 0)) then
                 cycle
              endif

              if (i==1) then
                 b0_ptr = c_loc(rhs%b0)
                 call spmGenRHS( SpmRhsRndX, rhs%nrhs, matrix%spm, &
                      & c_null_ptr, rhs%n, b0_ptr, rhs%n, info )
              end if

              ! Set the sequential scheduler to enable multiple solves in parallel with the same matrix
              if ( .not. params%multirhs ) then
                 matrix%iparm(IPARM_SCHEDULER) = PastixSchedSequential
              end if

              rhs%b(:) = rhs%b0(:)
              b_ptr = c_loc(rhs%b)

              ! Cannot be called in parallel with the same matrix for now
              call pastix_subtask_applyorder( matrix%pastix_data, SpmComplex64, PastixDirForward, &
                   &                          rhs%n, rhs%nrhs, b_ptr, rhs%n, info )

              rhs%x(:) = rhs%b(:)

           end do
        end do

        !
        !$OMP PARALLEL NUM_THREADS(params%nb_solv_omp) DEFAULT(NONE) &
        !$OMP SHARED(params, sla_lap, k, i)              &
        !$OMP SHARED(ila_thrsz, ila_thrmn, ila_thrmx)    &
        !$OMP SHARED(dla_thread_stats, iparm, dparm)     &
        !$OMP PRIVATE(th, im, ir, j )                    &
        !$OMP PRIVATE(matrix, rhs, x_ptr, b_ptr, b0_ptr) &
        !$OMP PRIVATE(info)

        !$OMP DO SCHEDULE(STATIC,1)
        solve_loop: do th = 1, params%nb_solv_omp
           solve_loop2: do j = ila_thrmn(th), ila_thrmx(th)

              if ( params%multirhs ) then
                 im = j
                 ir = 1
              else
                 im = ((j-1) / params%nb_fact_thd) + 1
                 ir = mod( j-1, params%nb_fact_thd ) + 1
              end if
              matrix => sla_lap(im)
              rhs    => matrix%rhs( ir )

              if ( .not. (rhs%nrhs .gt. 0)) then
                 cycle
              endif

              !   2- The right hand side
              x_ptr = c_loc(rhs%x)
              b_ptr = c_loc(rhs%b)

              ! 4- Solve the problem
              call pastix_subtask_solve( matrix%pastix_data, rhs%nrhs, &
                   & x_ptr, matrix%spm%n, info )

              dla_thread_stats(th, MULTILAP_SOLV_TIME) =  &
                   & dla_thread_stats(th,MULTILAP_SOLV_TIME) + matrix%dparm(DPARM_SOLV_TIME)

              ! 5- Refine the solution
              call pastix_subtask_refine(            &
                   & matrix%pastix_data,             &
                   & matrix%spm%n, rhs%nrhs,         &
                   & b_ptr, matrix%spm%n,            &
                   & x_ptr, matrix%spm%n, info )

              if ( .not. params%multirhs ) then
                 matrix%iparm(IPARM_SCHEDULER) = PastixSchedSequential
              end if

           end do solve_loop2
        end do solve_loop
        !$OMP END DO
        !$OMP END PARALLEL

        ! Apply the backward permutation on b and x
        do im = 1, params%nb_mat
           matrix => sla_lap(im)
           do j = 1, matrix%nsys
              rhs => matrix%rhs(j)

              if ( .not. (rhs%nrhs .gt. 0)) then
                 cycle
              endif

              x_ptr = c_loc(rhs%x)
              b_ptr = c_loc(rhs%b)

              ! Cannot be called in parallel with the same matrix for now
              call pastix_subtask_applyorder( matrix%pastix_data, SpmComplex64, PastixDirBackward, &
                   &                          rhs%n, rhs%nrhs, b_ptr, rhs%n, info )

              call pastix_subtask_applyorder( matrix%pastix_data, SpmComplex64, PastixDirBackward, &
                   &                          rhs%n, rhs%nrhs, x_ptr, rhs%n, info )

              ! Restore the initial scheduler
              if ( .not. params%multirhs ) then
                 matrix%iparm(IPARM_SCHEDULER) = iparm(IPARM_SCHEDULER)
              end if

           end do
        end do

        if (params%output) then
           do th = 1, params%nb_solv_omp
              do j = ila_thrmn(th), ila_thrmx(th)

                 if ( params%multirhs ) then
                    im = j
                    ir = 1
                 else
                    im = ((j-1) / params%nb_fact_thd) + 1
                    ir = mod( j-1, params%nb_fact_thd ) + 1
                 end if
                 matrix => sla_lap(im)
                 rhs    => matrix%rhs( ir )

                 if ( .not. (rhs%nrhs .gt. 0)) then
                    cycle
                 endif

                 write(6,*) '!--------------------------------------------------------------------!'
                 write(6,*) ' Check results for system ', im, ir

                 !
                 ! Check the solution
                 !
                 call spmPrintInfo( matrix%spm )

                 b_ptr = c_loc(rhs%b)
                 x_ptr = c_loc(rhs%x)

                 call spmCheckAxb( matrix%dparm(DPARM_EPSILON_REFINEMENT), rhs%nrhs, &
                      & matrix%spm,    &
                      & c_null_ptr, rhs%n, &
                      & b_ptr,      rhs%n, &
                      & x_ptr,      rhs%n, info )

                 ginfo = ginfo + info
              end do
           end do

           write(6,*) '!--------------------------------------------------------------------!'
           do th = 1, params%nb_fact_omp
              write(6,*) ' Thread ', th
              write(6,*) ' Time for solution      ', dla_thread_stats(th,MULTILAP_SOLV_TIME)
           end do
        end if
     end do ! End of iteration loop

     if (params%output) then
        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) ' Outer iterate ', k
        write(6,*) ' Time for solution      ', maxval(dla_thread_stats(:,MULTILAP_SOLV_TIME))
     end if

     dla_final_stats(MULTILAP_SOLV_TIME) = dla_final_stats(MULTILAP_SOLV_TIME) + &
          & maxval(dla_thread_stats(:,MULTILAP_SOLV_TIME))

     deallocate(dla_thread_stats)

     ! Free memory
     do im=1, params%nb_mat
        ! 6- Destroy the C data structure
        call spmExit( sla_lap(im)%spm )
     end do
  end do main_loop

  ! Destroy the matrices structures
  do im = 1, params%nb_mat
     !
     matrix => sla_lap(im)

     call pastixFinalize( matrix%pastix_data )

     do ir = 1, matrix%nsys
        rhs => matrix%rhs(ir)

        if ( rhs%nrhs .gt. 0 ) then
           if ( params%output ) then
              ! Allocate a backup of b that will be destroyed by the check
              deallocate(rhs%b0)
           end if
           deallocate(rhs%b)
           deallocate(rhs%x)
        end if
     end do
     deallocate( matrix%rhs )
     deallocate( matrix%spm )
     deallocate( matrix%bindtab )
  end do

  deallocate(sla_lap)
  if (allocated(ila_thrsz)) deallocate(ila_thrsz)
  if (allocated(ila_thrmn)) deallocate(ila_thrmn)
  if (allocated(ila_thrmx)) deallocate(ila_thrmx)

  write(6,*) '!====================================================================!'

  write(6,*) ' Overall final statistics'
  write(6,*) ' Time for analysis      ', dla_final_stats(MULTILAP_ANALYZE_TIME)
  write(6,*) ' Time for factorization ', dla_final_stats(MULTILAP_FACT_TIME)
  write(6,*) ' GFlops/s for fact      ', dla_final_stats(MULTILAP_FACT_FLOPS) / params%nb_outit
  write(6,*) ' Time for solution      ', dla_final_stats(MULTILAP_SOLV_TIME)
  write(6,*) '!====================================================================!'

  call exit(ginfo)

contains

  !
  ! @brief Initialize the testing parameters based on the input configuration file
  !
  ! @param[output] params
  !
  subroutine multilap_init( params )
    type(multilap_param), intent(out), target :: params
    integer                                   :: val, dim1, dim2, dim3
    integer                                   :: nbthrd_per_instance

    ! Read a dummy line.
    read( 5, * )

    ! Read if we enable/disable the verbose mode
    read( 5, * ) val
    params%output = .not. ( val .eq. 0 )

    ! Read if we enable/disable the check and correct
    read( 5, * ) val
    params%checkmat = .not. ( val .eq. 0 )

    ! Read if we enable/disable the multi-threaded solve
    read( 5, * ) val
    params%multirhs = .not. ( val .eq. 0 )

    ! Read the total number of threads
    read( 5, * ) params%nb_thrd

    if( params%nb_thrd < 1 ) then
       write( 6, fmt = 9999 ) 'nb_thrd', params%nb_thrd, params%nb_thrd, 1
       call exit(1)
    end if

    ! Read the number of pastix instances
    read( 5, * ) params%nb_fact_omp

    if( params%nb_fact_omp .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_fact_omp', params%nb_fact_omp, 1
       call exit(1)
    else if( params%nb_fact_omp .gt. params%nb_thrd) then
       write( 6, fmt = 9998 ) 'nb_fact_omp', params%nb_fact_omp, params%nb_thrd
       call exit(1)
    endif

    params%nb_fact_thd = params%nb_thrd / params%nb_fact_omp
    if ( (params%nb_fact_thd * params%nb_fact_omp) .ne. params%nb_thrd ) then
       write( 6, * ) 'nb_thrd (', params%nb_thrd, ')must be a multiple of nb_pastix (', params%nb_fact_omp, ')'
       call exit(1)
    end if

    ! Read the number of outer most iterations
    read( 5, * ) params%nb_outit

    if( params%nb_outit .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_outit', params%nb_outit, 1
       call exit(1)
    end if

    ! Read the number of matrices to factorize
    read( 5, * ) params%nb_mat

    if( params%nb_mat .lt. params%nb_fact_omp ) then
       write( 6, fmt = 9999 ) 'nb_mat', params%nb_mat, params%nb_fact_omp
       call exit(1)
    endif

    ! Read the number of RHS per iteration
    read( 5, * ) params%nb_rhs

    if( params%nb_rhs .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_rhs', params%nb_rhs, 1
       call exit(1)
    end if

    ! Read the number of solve per iteration
    read( 5, * ) params%nb_solve

    if( params%nb_solve .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_solve', params%nb_solve, 1
       call exit(1)
    end if

    ! Read the dimensions of the laplacian
    read( 5, * ) dim1

    if( dim1 .lt. 1 ) then
       write( 6, fmt = 9999 ) 'dim1', dim1, 1
       call exit(1)
    end if

    read( 5, * ) dim2

    if( dim2 .lt. 1 ) then
       write( 6, fmt = 9999 ) 'dim2', dim2, 1
       call exit(1)
    end if

    read( 5, * ) dim3

    if( dim3 .lt. 1 ) then
       write( 6, fmt = 9999 ) 'dim3', dim3, 1
       call exit(1)
    end if

    if ( params%multirhs ) then
       params%nb_sys      = params%nb_mat
       params%nb_solv_omp = params%nb_fact_omp
       params%nb_solv_thd = params%nb_fact_thd
    else
       params%nb_sys      = params%nb_mat * params%nb_fact_thd
       params%nb_solv_omp = params%nb_thrd
       params%nb_solv_thd = 1
    endif

    params%dim1 = dim1
    params%dim2 = dim2
    params%dim3 = dim3
    params%n    = dim1 * dim2 * dim3
    params%nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)

    write(6,*) '!--------------------------------------------------------------------!'
    write(6,*) '!           Multiple Laplacian testing configuration                 !'
    write(6,*) '!--------------------------------------------------------------------!'
    write(6,*) ' Nb of threads                               = ', params%nb_thrd
    write(6,*) ' Nb of PaStiX instances                      = ', params%nb_fact_omp
    write(6,*) ' Nb of outer iterations                      = ', params%nb_outit
    write(6,*) ' Nb of distinct matrices                     = ', params%nb_mat
    write(6,*) ' Nb of RHS to solve per matrix               = ', params%nb_rhs
    write(6,*) ' Nb of solve phase to perform per RHS        = ', params%nb_solve
    write(6,*) ' Size of each matrix                         = ', params%n, &
         &     '(', dim1, ' x ', dim2, ' x ', dim3, ')'
    write(6,*) ' Nbr of non zero entries per matrix          = ', params%nnz

    if ( params%multirhs ) then
       write(6,*) ' The multirhs mode is enabled'
    else
       write(6,*) ' The multirhs mode is disabled'
    endif

9997 format( A6, ' (', I6, ') is not a multiple of ', A6, ' (', I6 , ')' )
9998 format( ' Invalid input value: ', A8, '=', I6, '; must be <=', I6 )
9999 format( ' Invalid input value: ', A8, '=', I6, '; must be >=', I6 )

  end subroutine multilap_init

  !
  ! Generate a single laplacian spm matrix with the parameters stores
  ! in the matrix description
  !
  subroutine multilap_genOneMatrix( params, matrix, ib_out, ib )
    type(multilap_param),       intent(in),    target :: params
    type(sys_lin),              intent(inout), target :: matrix
    integer,                    intent(in)            :: ib, ib_out
    integer(kind=pastix_int_t)                        :: dim1, dim2, dim3, n, nnz
    integer                                           :: i, j, k, l

    !
    ! Laplacian dimensions
    !
    dim1   = params%dim1
    dim2   = params%dim2
    dim3   = params%dim3

    allocate(matrix%spm)
    allocate(matrix%rowptr(params%nnz))
    allocate(matrix%colptr(params%nnz))
    allocate(matrix%values(params%nnz))

    l = 1
    do i=1,dim1
       do j=1,dim2
          do k=1,dim3
             matrix%rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
             matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
             matrix%values(l) = 6.

             if (i == 1) then
                matrix%values(l) = matrix%values(l) - 1.
             end if
             if (i == dim1) then
                matrix%values(l) = matrix%values(l) - 1.
             end if
             if (j == 1) then
                matrix%values(l) = matrix%values(l) - 1.
             end if
             if (j == dim2) then
                matrix%values(l) = matrix%values(l) - 1.
             end if
             if (k == 1) then
                matrix%values(l) = matrix%values(l) - 1.
             end if
             if (k == dim3) then
                matrix%values(l) = matrix%values(l) - 1.
             end if

             matrix%values(l) = matrix%values(l) * 8.
             l = l + 1

             if (i < dim1) then
                matrix%rowptr(l) =  i    + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%values(l) = - 1. - 1. * i
                l = l + 1
             end if
             if (j < dim2) then
                matrix%rowptr(l) = (i-1) + dim1 *  j    + dim1 * dim2 * (k-1) + 1
                matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%values(l) = - 1. - 1. * i
                l = l + 1
             end if
             if (k < dim3) then
                matrix%rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 *  k    + 1
                matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%values(l) = -1. - 1. * i
                l = l + 1
             end if
          end do
       end do
    end do

    matrix%values(:) = matrix%values(:)*(dble(ib) * dble(ib_out) / 4.0)

    !
    ! Create the spm out of the internal data
    !
    call spmInit( matrix%spm )
    matrix%spm%mtxtype = SpmHermitian
    matrix%spm%flttype = SpmComplex64
    matrix%spm%fmttype = SpmIJV
    matrix%spm%n       = params%n
    matrix%spm%nnz     = params%nnz
    matrix%spm%dof     = 1
    matrix%spm%rowptr  = c_loc(matrix%rowptr)
    matrix%spm%colptr  = c_loc(matrix%colptr)
    matrix%spm%values  = c_loc(matrix%values)

    call spmUpdateComputedFields( matrix%spm )

    if (params%checkmat) then
       call spmCheckAndCorrect( matrix%spm, matrix%spm2 )
       if (.not. c_associated(c_loc(matrix%spm), c_loc(matrix%spm2))) then
          deallocate(matrix%rowptr)
          deallocate(matrix%colptr)
          deallocate(matrix%values)

          matrix%spm%rowptr = c_null_ptr
          matrix%spm%colptr = c_null_ptr
          matrix%spm%values = c_null_ptr

          call spmExit( matrix%spm )
          matrix%spm => matrix%spm2
       end if
    else
       call spmConvert(SpmCSC, matrix%spm, info)
    endif
  end subroutine multilap_genOneMatrix

  !
  subroutine split_parall(id_size, id_thr)

    integer, intent(in) :: id_size
    integer, intent(in) :: id_thr

    if (allocated(ila_thrsz)) deallocate(ila_thrsz)
    allocate(ila_thrsz(id_thr))
    if (allocated(ila_thrmn)) deallocate(ila_thrmn)
    allocate(ila_thrmn(id_thr))
    if (allocated(ila_thrmx)) deallocate(ila_thrmx)
    allocate(ila_thrmx(id_thr))

    ila_thrsz(:) = id_size / id_thr
    ila_thrsz(1:id_size-id_thr*ila_thrsz(1)) =&
         &  ila_thrsz(1:id_size-id_thr*ila_thrsz(1)) + 1

    ila_thrmn(1) = 1
    do i = 1, id_thr-1
       ila_thrmx(i)   = ila_thrmn(i) + ila_thrsz(i)-1
       ila_thrmn(i+1) = ila_thrmx(i) + 1
    end do
    ila_thrmx(id_thr) = id_size

  end subroutine split_parall

end program flaplacian


