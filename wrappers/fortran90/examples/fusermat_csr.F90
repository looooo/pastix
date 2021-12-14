!
! @file fusermat_csr.F90
!
! Fortran 90 example with a matrix defined by the user matrix.
! Additionaly, this example shows how to cheat with the solver to
! provide a CSR matrix without paying the cost of transposing it.
!
! @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.4.0
! @author Mathieu Faverge
! @author Tony Delarue
! @date 2024-07-05
!
program fusermat_csr
  use iso_c_binding
  use spmf
  use pastixf
#if defined(PASTIX_WITH_MPI)
  use mpi_f08
#endif
  implicit none

  type userMatCSR
     integer                                               :: n        ! Number of vertices
     integer                                               :: nnz      ! Number of non-zeroes in the graph
     integer(kind=pastix_int_t), dimension(:), allocatable :: rows     ! Array of size n+1 of indirections to rows for each vertex
     integer(kind=pastix_int_t), dimension(:), allocatable :: cols     ! Array of size nnz that corresponds to the global column numbering
     integer(kind=pastix_int_t), dimension(:), allocatable :: loc2glob ! Corresponding numbering from local to global
     real(kind=c_double),        dimension(:), allocatable :: vals     ! Values stored in the matrix
  end type userMatCSR

  type(pastix_data_t), pointer                      :: pastix_data

  type(userMatCSR),    pointer :: userMat  ! Local format of the matrix
  type(spmatrix_t),    pointer :: spm      ! Spm format of the matrix

  integer(kind=pastix_int_t), target                       :: iparm(IPARM_SIZE)
  real(kind=c_double),        target                       :: dparm(DPARM_SIZE)
  real(kind=c_double)                                      :: normA
  integer(c_int)                                           :: info
  integer(kind=pastix_int_t)                               :: nrhs
  real(kind=c_double), dimension(:,:), allocatable, target :: x, b

  logical,             parameter :: use_trans = .TRUE.
  real(kind=c_double), parameter :: alpha     = 10.
  real(kind=c_double), parameter :: beta      = 0.5

  !
  ! Solve the problem
  !

  ! 1- Initialize the parameters and the solver
  ! (Done before any calls to spm to automatically intialize MPI if needed)
  call pastixInitParam( iparm, dparm )

  ! Here you can decide to change the parameters from the defaults ones.
  ! A set of IPARM is listed below
  ! All the IPARM are listed here : https://solverstack.gitlabpages.inria.fr/pastix/group__pastix__api.html

  ! Verbose mode - Default: PastixVerboseNo
  ! Possible values : PastixVerboseNot, PastixVerboseNo, PastixVerboseYes
  iparm(IPARM_VERBOSE) = PastixVerboseNo

  ! Factorization mode - Default: PastixFactLU
  ! Possible values : PastixFactLLH, PastixFactLDLT, PastixFactLU, PastixFactLLT, PastixFactLDLH
  iparm(IPARM_FACTORIZATION) = PastixFactLU

  ! Refinement mode - Default: PastixRefineGMRES
  ! Possible values : PastixRefineGMRES, PastixRefineCG, PastixRefineSR, PastixRefineBiCGSTAB
  iparm(IPARM_REFINEMENT) = PastixRefineGMRES

  ! Scheduler mode - Default: PastixSchedDynamic
  ! Possible values : PastixSchedSequential, PastixSchedStatic, PastixSchedParsec, PastixSchedStarpu, PastixSchedDynamic
  iparm(IPARM_SCHEDULER) = PastixSchedDynamic

  call pastixInit( pastix_data, MPI_COMM_WORLD, iparm, dparm )

  ! Initialize the problem
  call create3DLaplacianMatrix( userMat, 4, 4, 4, alpha, beta )

  ! convert the internal matrix format to csc (see other examples to generate directly in the spm format)
  if ( use_trans ) then
     call convert2spmNoCopy( userMat, spm )
     iparm(IPARM_TRANSPOSE_SOLVE) = PastixTrans
  else
     call convert2spmCopy( userMat, spm )
     call cleanUserMat( userMat ) ! Can be done as soon as it is converted if we don't need it
  end if

  ! Scale A for better stability with low-rank computations
  call spmNorm( SpmFrobeniusNorm, spm, normA )
  call spmScal( 1. / normA, spm )

  ! 2- Analyze the problem
  call pastix_task_analyze( pastix_data, spm, info )
  if ( info .ne. 0 ) then
     error stop "pastix_task_analyze failed"
  end if

  ! 3- Factorize the matrix
  call pastix_task_numfact( pastix_data, spm, info )
  if ( info .ne. 0 ) then
     error stop "pastix_task_numfact failed"
  end if

  !   2- Generate 10 random right hand side
  !
  ! Here since the degree of freedom is 1, n and nexp are equal and
  ! can both be used, but it is a better practice to use nexp if
  ! dof changes
  !
  nrhs = 2
  allocate(x(spm%nexp, nrhs))
  allocate(b(spm%nexp, nrhs))

  call spmGenRHS( SpmRhsRndB, nrhs, spm, B=b, ldb=spm%nexp, info=info )
  if ( info .ne. 0 ) then
     error stop "spmGenRHS failed"
  end if
  x = b

  ! 4- Solve the problem
  call pastix_task_solve( pastix_data, spm%nexp, nrhs, x, spm%nexp, info )
  if ( info .ne. 0 ) then
     error stop "pastix_task_solve failed"
  end if

  ! 5- Refine the solution if needed (More optimized version of solver+refinement is presented in fmultilap.F90)
  call pastix_task_refine( pastix_data, spm%nexp, nrhs, b, spm%nexp, x, spm%nexp, info )
  if ( info .ne. 0 ) then
     error stop "pastix_rask_refine failed"
  end if

  !
  ! Check the solution
  !
  call spmCheckAxb( dparm(DPARM_EPSILON_REFINEMENT), nrhs, spm, &
       B=b, ldb=spm%nexp, X=x, ldx=spm%nexp, info=info )
  if ( info .ne. 0 ) then
     error stop "spmCheckAxB failed"
  end if

  ! Set fortran pointers to null before calling spmExit()
  call cleanSpm( spm, use_trans )
  if ( use_trans ) then
     call cleanUserMat( userMat ) ! Can be done as soon as it is converted if we don't need it
  end if

  deallocate(x)
  deallocate(b)

  ! 6- Destroy the C data structure (Should be last if used to call MPI_Finalize)
  call pastixFinalize( pastix_data )

contains

  !
  ! This routines creates a random sparse matrix with n vertices, and
  ! nnz non-zeroes. The values are set to corresponds to a Laplacian
  ! with alpha times the degree on the diagonal, and -beta on the edges.
  !
  subroutine create3DLaplacianMatrix( userMat, dim1, dim2, dim3, alpha, beta )

    implicit none

    ! Parameters
    type(userMatCSR),    intent(inout), pointer :: userMat
    integer,             intent(in)             :: dim1, dim2, dim3
    real(kind=c_double), intent(in)             :: alpha, beta

    ! Variables
    integer :: comm_rank, comm_size
    integer :: n, nnz
    integer :: ldim1, k1, k2
    integer :: i, j, k, l, di, dj, dk ! degrees on each dimension

    !
    ! Compute the local number of unknowns by spliting along the
    ! first dimension of the 3D grid if distributed
    !
#if defined(PASTIX_WITH_MPI)
    integer :: ierr
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, comm_rank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
#else
    comm_rank = 0
    comm_size = 1
#endif
    ldim1 = dim1 / comm_size;
    k1 = ldim1 *  comm_rank    + 1;
    k2 = ldim1 * (comm_rank+1);

    if (comm_rank < modulo( dim1, comm_size )) then
       k1 = k1 + comm_rank
       k2 = k2 + comm_rank
    end if
    if ((comm_rank+1) < modulo( dim1, comm_size )) then
       k2    = k2 + 1
       ldim1 = ldim1 + 1
    end if

    n = ldim1 * dim2 * dim3

    !
    ! Compute the local number of edges without forgetting the edges
    ! to connect to the neighbors
    !
    nnz = n ! Diagonal elements

    ! Add the beta edges - Sides of the squares
    nnz = nnz + (ldim1 - 1) *  dim2 * dim3;
    nnz = nnz + ( dim2 - 1) * ldim1 * dim3;
    nnz = nnz + ( dim3 - 1) * ldim1 * dim2;

    ! Add connexion edges
    if (comm_rank < min( dim1-1, comm_size-1 )) then
       nnz = nnz + dim2 * dim3
    end if

    ! Initialize the userMat structure
    allocate(userMat)
    userMat%n   = n
    userMat%nnz = nnz
    allocate( userMat%rows( n+1 ) )
    allocate( userMat%cols( nnz ) )
    allocate( userMat%vals( nnz ) )
    if (comm_size > 1) then
       allocate( userMat%loc2glob( n ) )
    end if

    ! Now let's build the matrix with the last dimension first
    n   = 1
    nnz = 1
    userMat%rows(n) = 1
    do i=1,dim3
       di = 0

       ! +1 after the first layer
       if ( i > 1 ) then
          di = di + 1;
       end if

       ! -1 on the last layer
       if ( i < dim3 ) then
          di = di + 1;
       end if

       do j=1,dim2
          dj = 0

          ! +1 after the first layer
          if ( j > 1 ) then
             dj = dj + 1;
          end if

          ! -1 on the last layer
          if ( j < dim2 ) then
             dj = dj + 1;
          end if

          l = (i-1) * dim1 * dim2 + (j-1) * dim1
          do k=k1,k2
             dk = 0

             ! +1 after the first layer
             if ( k > 1 ) then
                dk = dk + 1;
             end if

             ! -1 on the last layer
             if ( k < dim1 ) then
                dk = dk + 1;
             end if

             userMat%rows(n+1) = userMat%rows(n)
             if (comm_size > 1) then
                userMat%loc2glob( n ) = l + k
             end if

             ! Diagonal element
             userMat%rows( n+1 ) = userMat%rows( n+1 ) + 1
             userMat%cols( nnz ) = l + k
             userMat%vals( nnz ) = (di + dk + dj) * alpha;
             nnz = nnz + 1

             ! Connexion along dimension 1
             if ( k < dim1 ) then
                userMat%rows( n+1 ) = userMat%rows( n+1 ) + 1
                userMat%cols( nnz ) = l + k + 1
                userMat%vals( nnz ) = beta;
                nnz = nnz + 1
             end if

             ! Connexion along dimension 2
             if ( j < dim2 ) then
                userMat%rows( n+1 ) = userMat%rows( n+1 ) + 1
                userMat%cols( nnz ) = l + k + dim1
                userMat%vals( nnz ) = beta;
                nnz = nnz + 1
             end if

             ! Connexion along dimension 3
             if ( i < dim3 ) then
                userMat%rows( n+1 ) = userMat%rows( n+1 ) + 1
                userMat%cols( nnz ) = l + k + dim1 * dim2
                userMat%vals( nnz ) = beta;
                nnz = nnz + 1
             end if

             n = n + 1

          end do
       end do
    end do

  end subroutine create3DLaplacianMatrix

  ! Example using fortran allocated arrays
  ! If you are looking for a version with the arrays allocated
  ! thanks to the spm library, please refer to flaplacian.F90
  !
  subroutine convert2spmNoCopy( userMat, spm )
    implicit none

    ! Parameters
    type(userMatCSR), intent(inout), pointer :: userMat
    type(spmatrix_t), intent(out),   pointer :: spm

    allocate( spm )
    call spmInitDist( spm, MPI_COMM_WORLD )
    !
    ! Note that we cheat and use the symmetric property of the
    ! matrix to set it as CSC and avoid incorrrect pointers
    ! management bewteen C and Fortran, so we don't have to call
    ! check and correct that may deallocate fortran pointers.
    !
    spm%baseval = 1              ! Arrays are initialized 1-based
    spm%mtxtype = SpmSymmetric   ! Only one triangle of the matrix is given
    spm%flttype = SpmDouble      ! Values are stores in double
    spm%fmttype = SpmCSC         ! Format in CSC
    spm%n       = userMat%n      ! Local number of unknowns
    spm%nnz     = userMat%nnz    ! Local number of non zeroes
    spm%dof     = 1              ! Degree of freedom per unknown

    spm%colptr = c_loc( userMat%rows )
    spm%rowptr = c_loc( userMat%cols )
    spm%values = c_loc( userMat%vals )

    if (allocated( userMat%loc2glob )) then ! or if MPI_Comm_Size > 1
       spm%replicated = 0 ! Each node owns a portion of the matrix
       spm%loc2glob = c_loc( userMat%loc2glob )
    else
       spm%replicated = 1 ! The matrix is fully replicated on each node (or a single node is used)
    end if

    ! Compute the global information within the structure between the nodes involved
    call spmUpdateComputedFields( spm )

    ! We are sure that the matrix is correct, so we don't call
    ! CheckAndCorrect to not create C/Fortran conflict on the arrays:
    !
    !  - Using CSC format
    !  - No duplicated entries
    !  - entries for all diagonal values
    !

    ! Print basic information for debug
    call spmPrintInfo( spm )
  end subroutine convert2spmNoCopy

  ! Example using fortran allocated arrays
  ! If you are looking for a version with the arrays allocated
  ! thanks to the spm library, please refer to flaplacian.F90
  !
  subroutine convert2spmCopy( userMat, spm )
    implicit none

    ! Parameters
    type(userMatCSR), intent(inout), pointer :: userMat
    type(spmatrix_t), intent(out),   pointer :: spm

    ! Variables
    type(spmatrix_t),                      pointer :: spm2
    integer(kind=spm_int_t), dimension(:), pointer :: rowptr
    integer(kind=spm_int_t), dimension(:), pointer :: colptr
    real(kind=c_double),     dimension(:), pointer :: values
    integer(kind=spm_int_t), dimension(:), pointer :: loc2glob

    allocate( spm )
    call spmInit( spm )        ! By default using MPI_COMM_WORLD, switch to spmInitDist to select the communicator
    spm%baseval = 1            ! Arrays are initialized 1-based
    spm%mtxtype = SpmSymmetric ! Only one triangle of the matrix is given
    spm%flttype = SpmDouble    ! Values are stores in double
    spm%fmttype = SpmCSR       ! Format in CSC
    spm%n       = userMat%n    ! Local number of unknowns
    spm%nnz     = userMat%nnz  ! Local number of non zeroes
    spm%dof     = 1            ! Degree of freedom per unknown

    if (allocated( userMat%loc2glob )) then ! or if MPI_Comm_Size > 1
       spm%replicated = 0 ! Each node owns a portion of the matrix
    else
       spm%replicated = 1 ! The matrix is fully replicated on each node (or a single node is used)
    end if

    ! Compute the global information within the structure between the nodes involved
    call spmUpdateComputedFields( spm )

    call spmAlloc( spm ) ! Allocate the arrays through C calls

    ! get the arrays pointers in fortran
    call spmGetArray( spm, colptr=colptr, rowptr=rowptr, dvalues=values )

    colptr(:) = userMat%cols(:)
    rowptr(:) = userMat%rows(:)
    values(:) = userMat%vals(:)

    if (allocated( userMat%loc2glob )) then ! or if MPI_Comm_Size > 1
       call spmGetArray( spm, loc2glob = loc2glob )
       loc2glob(:) = userMat%loc2glob(:)
    end if

    ! We make sure that the matrix is correct
    !
    !  - Using CSC format
    !  - No duplicated entries
    !  - entries for all diagonal values
    !
    ! WARNING: this call can be costly in time and memory as it may
    ! double the memory to work on the matrix and it goes through
    ! all edges to check the correctness of the spm, thus being
    ! quite slow. If it can be avoided, it should.
    !
    ! When the matrix is correct and the only issue is the CSR
    ! storage, refer so siple_trans.c to solve the transpose
    ! problem.
    !
    allocate( spm2 )
    call spmCheckAndCorrect( spm, spm2, info )
    if ( info .ne. 0 ) then
       call spmExit( spm )
       spm = spm2
    end if
    deallocate( spm2 )

    ! Print basic information for debug
    call spmPrintInfo( spm )
  end subroutine convert2spmCopy

  subroutine cleanUserMat( userMat )
    implicit none

    type(userMatCSR), intent(inout), pointer :: userMat

    deallocate( userMat%cols )
    deallocate( userMat%rows )
    deallocate( userMat%vals )

    if (allocated( userMat%loc2glob )) then ! or if MPI_Comm_Size > 1
       deallocate( userMat%loc2glob )
    end if

    deallocate( userMat )
  end subroutine cleanUserMat

  subroutine cleanSpm( spm, fpointers )
    implicit none

    type(spmatrix_t), intent(inout), pointer :: spm
    logical,          intent(in)             :: fpointers

    if ( fpointers ) then
       spm%colptr = c_null_ptr
       spm%rowptr = c_null_ptr
       spm%values = c_null_ptr
       spm%loc2glob = c_null_ptr
    end if
    call spmExit( spm )
    deallocate(spm)

  end subroutine cleanSpm

end program fusermat_csr
