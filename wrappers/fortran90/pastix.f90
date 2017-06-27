module pastix_mod
  use iso_c_binding
  implicit none

  ! iparm enum
  integer, parameter :: iparm_verbose               = 0
  integer, parameter :: iparm_io_strategy           = iparm_verbose               + 1
  integer, parameter :: iparm_nnzeros               = iparm_io_strategy           + 1
  integer, parameter :: iparm_nnzeros_block_local   = iparm_nnzeros               + 1
  integer, parameter :: iparm_allocated_terms       = iparm_nnzeros_block_local   + 1
  integer, parameter :: iparm_produce_stats         = iparm_allocated_terms       + 1
  integer, parameter :: iparm_mc64                  = iparm_produce_stats         + 1
  integer, parameter :: iparm_ordering              = iparm_mc64                  + 1
  integer, parameter :: iparm_ordering_default      = iparm_ordering              + 1
  integer, parameter :: iparm_scotch_switch_level   = iparm_ordering_default      + 1
  integer, parameter :: iparm_scotch_cmin           = iparm_scotch_switch_level   + 1
  integer, parameter :: iparm_scotch_cmax           = iparm_scotch_cmin           + 1
  integer, parameter :: iparm_scotch_frat           = iparm_scotch_cmax           + 1
  integer, parameter :: iparm_metis_ctype           = iparm_scotch_frat           + 1
  integer, parameter :: iparm_metis_rtype           = iparm_metis_ctype           + 1
  integer, parameter :: iparm_metis_no2hop          = iparm_metis_rtype           + 1
  integer, parameter :: iparm_metis_nseps           = iparm_metis_no2hop          + 1
  integer, parameter :: iparm_metis_niter           = iparm_metis_nseps           + 1
  integer, parameter :: iparm_metis_ufactor         = iparm_metis_niter           + 1
  integer, parameter :: iparm_metis_compress        = iparm_metis_ufactor         + 1
  integer, parameter :: iparm_metis_ccorder         = iparm_metis_compress        + 1
  integer, parameter :: iparm_metis_pfactor         = iparm_metis_ccorder         + 1
  integer, parameter :: iparm_metis_seed            = iparm_metis_pfactor         + 1
  integer, parameter :: iparm_metis_dbglvl          = iparm_metis_seed            + 1
  integer, parameter :: iparm_sf_kass               = iparm_metis_dbglvl          + 1
  integer, parameter :: iparm_amalgamation_lvlblas  = iparm_sf_kass               + 1
  integer, parameter :: iparm_amalgamation_lvlcblk  = iparm_amalgamation_lvlblas  + 1
  integer, parameter :: iparm_reordering_split      = iparm_amalgamation_lvlcblk  + 1
  integer, parameter :: iparm_reordering_stop       = iparm_reordering_split      + 1
  integer, parameter :: iparm_min_blocksize         = iparm_reordering_stop       + 1
  integer, parameter :: iparm_max_blocksize         = iparm_min_blocksize         + 1
  integer, parameter :: iparm_distribution_level    = iparm_max_blocksize         + 1
  integer, parameter :: iparm_abs                   = iparm_distribution_level    + 1
  integer, parameter :: iparm_incomplete            = iparm_abs                   + 1
  integer, parameter :: iparm_level_of_fill         = iparm_incomplete            + 1
  integer, parameter :: iparm_factorization         = iparm_level_of_fill         + 1
  integer, parameter :: iparm_static_pivoting       = iparm_factorization         + 1
  integer, parameter :: iparm_inertia               = iparm_static_pivoting       + 1
  integer, parameter :: iparm_free_cscuser          = iparm_inertia               + 1
  integer, parameter :: iparm_refinement            = iparm_free_cscuser          + 1
  integer, parameter :: iparm_nbiter                = iparm_refinement            + 1
  integer, parameter :: iparm_itermax               = iparm_nbiter                + 1
  integer, parameter :: iparm_gmres_im              = iparm_itermax               + 1
  integer, parameter :: iparm_scheduler             = iparm_gmres_im              + 1
  integer, parameter :: iparm_thread_nbr            = iparm_scheduler             + 1
  integer, parameter :: iparm_autosplit_comm        = iparm_thread_nbr            + 1
  integer, parameter :: iparm_gpu_nbr               = iparm_autosplit_comm        + 1
  integer, parameter :: iparm_gpu_memory_percentage = iparm_gpu_nbr               + 1
  integer, parameter :: iparm_gpu_memory_block_size = iparm_gpu_memory_percentage + 1
  integer, parameter :: iparm_compress_min_width    = iparm_gpu_memory_block_size + 1
  integer, parameter :: iparm_compress_min_height   = iparm_compress_min_width    + 1
  integer, parameter :: iparm_compress_when         = iparm_compress_min_height   + 1
  integer, parameter :: iparm_compress_method       = iparm_compress_when         + 1
  integer, parameter :: iparm_thread_comm_mode      = iparm_compress_method       + 1
  integer, parameter :: iparm_modify_parameter      = iparm_thread_comm_mode      + 1
  integer, parameter :: iparm_start_task            = iparm_modify_parameter      + 1
  integer, parameter :: iparm_end_task              = iparm_start_task            + 1
  integer, parameter :: iparm_float_type            = iparm_end_task              + 1
  integer, parameter :: iparm_mtx_type              = iparm_float_type            + 1
  integer, parameter :: iparm_dof_nbr               = iparm_mtx_type              + 1
  integer, parameter :: iparm_size                  = iparm_dof_nbr               + 1

  ! dparm enum
  integer, parameter :: dparm_fill_in            = 0
  integer, parameter :: dparm_epsilon_refinement = dparm_fill_in            + 1
  integer, parameter :: dparm_relative_error     = dparm_epsilon_refinement + 1
  integer, parameter :: dparm_epsilon_magn_ctrl  = dparm_relative_error     + 1
  integer, parameter :: dparm_analyze_time       = dparm_epsilon_magn_ctrl  + 1
  integer, parameter :: dparm_pred_fact_time     = dparm_analyze_time       + 1
  integer, parameter :: dparm_fact_time          = dparm_pred_fact_time     + 1
  integer, parameter :: dparm_solv_time          = dparm_fact_time          + 1
  integer, parameter :: dparm_fact_flops         = dparm_solv_time          + 1
  integer, parameter :: dparm_fact_thflops       = dparm_fact_flops         + 1
  integer, parameter :: dparm_fact_rlflops       = dparm_fact_thflops       + 1
  integer, parameter :: dparm_solv_flops         = dparm_fact_rlflops       + 1
  integer, parameter :: dparm_solv_thflops       = dparm_solv_flops         + 1
  integer, parameter :: dparm_solv_rlflops       = dparm_solv_thflops       + 1
  integer, parameter :: dparm_refine_time        = dparm_solv_rlflops       + 1
  integer, parameter :: dparm_a_norm             = dparm_refine_time        + 1
  integer, parameter :: dparm_compress_tolerance = dparm_a_norm             + 1
  integer, parameter :: dparm_size               = dparm_compress_tolerance + 1

  integer, parameter :: PastixTaskInit       = 0
  integer, parameter :: PastixTaskOrdering   = 1
  integer, parameter :: PastixTaskSymbfact   = 2
  integer, parameter :: PastixTaskAnalyze    = 3
  integer, parameter :: PastixTaskNumfact    = 4
  integer, parameter :: PastixTaskSolve      = 5
  integer, parameter :: PastixTaskRefine     = 6
  integer, parameter :: PastixTaskClean      = 7

  ! class verbose:
  integer, parameter :: PastixVerboseNot   = 0
  integer, parameter :: PastixVerboseNo    = 1
  integer, parameter :: PastixVerboseYes   = 2

  ! class io:
  integer, parameter :: PastixIONo         = 0
  integer, parameter :: PastixIOLoad       = 1
  integer, parameter :: PastixIOSave       = 2
  integer, parameter :: PastixIOLoadGraph  = 4
  integer, parameter :: PastixIOSaveGraph  = 8
  integer, parameter :: PastixIOLoadCSC    = 16
  integer, parameter :: PastixIOSaveCSC    = 32

  ! class refine:
  integer, parameter :: PastixRefineGMRES    = 0
  integer, parameter :: PastixRefineCG       = 1
  integer, parameter :: PastixRefineSR       = 2
  integer, parameter :: PastixRefineBiCGSTAB = 3

  ! class coeftype:
  integer, parameter :: PastixPattern   = 0
  integer, parameter :: PastixFloat     = 2
  integer, parameter :: PastixDouble    = 3
  integer, parameter :: PastixComplex32 = 4
  integer, parameter :: PastixComplex64 = 5

  ! class fmttype:
  integer, parameter :: PastixCSC = 0
  integer, parameter :: PastixCSR = 1
  integer, parameter :: PastixIJV = 2

  ! class factotype:
  integer, parameter :: PastixFactLLT  = 0 ! Cholesky factorization
  integer, parameter :: PastixFactLDLT = 1 ! LDL^t factorization
  integer, parameter :: PastixFactLU   = 2 ! LU factorization
  integer, parameter :: PastixFactLDLH = 3 ! LDL^h factorization for complex matrices

  ! class scheduler:
  integer, parameter :: PastixSchedSequential = 0 ! Sequential
  integer, parameter :: PastixSchedStatic     = 1 ! Shared memory with static scheduler
  integer, parameter :: PastixSchedParsec     = 2 ! PaRSEC scheduler
  integer, parameter :: PastixSchedStarpu     = 3 ! StarPU scheduler
  integer, parameter :: PastixSchedDynamic    = 4 ! Shared memory with dynamic scheduler

  ! class order:
  integer, parameter :: PastixOrderScotch   = 0 ! Use Scotch ordering
  integer, parameter :: PastixOrderMetis    = 1 ! Use Metis ordering
  integer, parameter :: PastixOrderPersonal = 2 ! Apply user's permutation
  integer, parameter :: PastixOrderLoad     = 3 ! Load ordering from file
  integer, parameter :: PastixOrderPtScotch = 4 ! Use Pt-Scotch ordering
  integer, parameter :: PastixOrderParMetis = 5 ! Use ParMetis ordering

  ! class threadmode:
  integer, parameter :: PastixThreadMultiple = 1 ! All threads communicate
  integer, parameter :: PastixThreadFunneled = 2 ! One thread perform all the MPI Calls

  ! class error:
  integer, parameter :: PASTIX_SUCCESS            = 0  ! No error
  integer, parameter :: PASTIX_ERR_UNKNOWN        = 1  ! Unknown error
  integer, parameter :: PASTIX_ERR_ALLOC          = 2  ! Allocation error
  integer, parameter :: PASTIX_ERR_NOTIMPLEMENTED = 3  ! Not implemented feature
  integer, parameter :: PASTIX_ERR_OUTOFMEMORY    = 4  ! Not enough memory
  integer, parameter :: PASTIX_ERR_THREAD         = 5  ! Error with threads
  integer, parameter :: PASTIX_ERR_INTERNAL       = 6  ! Internal error
  integer, parameter :: PASTIX_ERR_BADPARAMETER   = 7  ! Bad parameters given
  integer, parameter :: PASTIX_ERR_FILE           = 8  ! Error in In/Out operations
  integer, parameter :: PASTIX_ERR_INTEGER_TYPE   = 9  ! Error with integer types
  integer, parameter :: PASTIX_ERR_IO             = 10 ! Error with input/output
  integer, parameter :: PASTIX_ERR_MPI            = 11 ! Error with MPI calls

  ! class compress_when:
  integer, parameter :: PastixCompressNever  = 0
  integer, parameter :: PastixCompressBegin  = 1
  integer, parameter :: PastixCompressEnd    = 2
  integer, parameter :: PastixCompressDuring = 3

  ! class compress_method:
  integer, parameter :: PastixCompressMethodSVD  = 0
  integer, parameter :: PastixCompressMethodRRQR = 1

  ! class driver:
  integer, parameter :: PastixDriverRSA        = 0 ! RSA Fortran driver
  integer, parameter :: PastixDriverHB         = 1 ! Harwell Boeing driver
  integer, parameter :: PastixDriverIJV        = 2 ! IJV Coordinate driver
  integer, parameter :: PastixDriverMM         = 3 ! Matrix Market C driver
  integer, parameter :: PastixDriverLaplacian  = 4 ! 3, 5, or 7 points Laplacian stencil generator
  integer, parameter :: PastixDriverXLaplacian = 5 ! 15-points Laplacian stencil generator
  integer, parameter :: PastixDriverGraph      = 6 ! Scotch Graph driver

  ! class rhstype:
  integer, parameter :: PastixRhsOne =  0
  integer, parameter :: PastixRhsI   =  1
  integer, parameter :: PastixRhsRndX = 2
  integer, parameter :: PastixRhsRndB = 3

  ! class layout:
  integer, parameter :: PastixRowMajor  = 101 ! Storage in row major order
  integer, parameter :: PastixColMajor  = 102 ! Storage in column major order

  ! class trans:
  integer, parameter :: PastixNoTrans   = 111 ! Use A
  integer, parameter :: PastixTrans     = 112 ! Use A^t
  integer, parameter :: PastixConjTrans = 113 ! Use conj(A^t)

  ! class mtxtype:
  integer, parameter :: PastixGeneral   = PastixNoTrans    ! The matrix is general
  integer, parameter :: PastixSymmetric = PastixTrans      ! The matrix is symmetric
  integer, parameter :: PastixHermitian = PastixConjTrans  ! The matrix is hermitian

  ! class uplo:
  integer, parameter :: PastixUpper      = 121 ! Use lower triangle of A
  integer, parameter :: PastixLower      = 122 ! Use upper triangle of A
  integer, parameter :: PastixUpperLower = 123 ! Use the full A

  ! class coefside:
  integer, parameter :: PastixLCoef  = 0 ! Coefficients of the lower triangular L are used
  integer, parameter :: PastixUCoef  = 1 ! Coefficients of the upper triangular U are used
  integer, parameter :: PastixLUCoef = 2 ! Coefficients of the upper/lower triangular U/L are used

  ! class diag:
  integer, parameter :: PastixNonUnit = 131 ! Diagonal is non unitary
  integer, parameter :: PastixUnit    = 132 ! Diagonal is unitary

  ! class side:
  integer, parameter :: PastixLeft  = 141 ! Apply operator on the left
  integer, parameter :: PastixRight = 142 ! Apply operator on the right

  ! class normtype:
  integer, parameter :: PastixOne       = 171 ! One norm:       max_j( sum_i( |a_{ij}| ) )
  integer, parameter :: PastixFrobenius = 174 ! Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) )
  integer, parameter :: PastixInf       = 175 ! Inifinite norm: max_i( sum_j( |a_{ij}| ) )
  integer, parameter :: PastixMax       = 177 ! Inifinite norm: max_{i,j}( | a_{ij} | )

  ! class direction:
  integer, parameter :: PastixForward  = 391 ! Forward direction
  integer, parameter :: PastixBackward = 392 ! Backward direction

  ! C structs converted to derived types.
  type, bind(c) :: pastix_spm_t
     integer(kind=c_int) :: mtxtype
     integer(c_int) :: flttype
     integer(c_int) :: fmttype
     integer(c_int) :: gN
     integer(c_int) :: n
     integer(c_int) :: gnnz
     integer(c_int) :: nnz
     integer(c_int) :: gNexp
     integer(c_int) :: nexp
     integer(c_int) :: gnnzexp
     integer(c_int) :: nnzexp
     integer(c_int) :: dof
     type(c_ptr) :: dofs
     integer(c_int) :: layout
     type(c_ptr) :: colptr
     type(c_ptr) :: rowptr
     type(c_ptr) :: loc2glob
     type(c_ptr) :: values
  end type pastix_spm_t

  ! Interfaces of the C functions.
  interface
     function pastix_c(pastix_data, pastix_comm, n, colptr, row, avals, perm, invp, &
          b, nrhs, iparm, dparm) &
          bind(c, name='pastix')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_c
       type(c_ptr) :: pastix_data
       integer(kind=c_int), value :: pastix_comm
       integer(c_int), value :: n
       type(c_ptr), value :: colptr
       type(c_ptr), value :: row
       type(c_ptr), value :: avals
       type(c_ptr), value :: perm
       type(c_ptr), value :: invp
       type(c_ptr), value :: b
       integer(c_int), value :: nrhs
       type(c_ptr), value :: iparm
       type(c_ptr), value :: dparm
     end function pastix_c
  end interface

  interface
     subroutine pastixInitParam_c(iparm, dparm) &
          bind(c, name='pastixInitParam')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: iparm
       type(c_ptr), value :: dparm
     end subroutine pastixInitParam_c
  end interface

  interface
     subroutine pastixInit_c(pastix_data, pastix_comm, iparm, dparm) &
          bind(c, name='pastixInit')
       use iso_c_binding
       implicit none
       type(c_ptr) :: pastix_data
       integer(kind=c_int), value :: pastix_comm
       type(c_ptr), value :: iparm
       type(c_ptr), value :: dparm
     end subroutine pastixInit_c
  end interface

  interface
     subroutine pastixFinalize_c(pastix_data) &
          bind(c, name='pastixFinalize')
       use iso_c_binding
       implicit none
       type(c_ptr) :: pastix_data
     end subroutine pastixFinalize_c
  end interface

  interface
     function pastix_task_analyze_c(pastix_data, spm) &
          bind(c, name='pastix_task_analyze')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_task_analyze_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_task_analyze_c
  end interface

  interface
     function pastix_task_numfact_c(pastix_data, spm) &
          bind(c, name='pastix_task_numfact')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_task_numfact_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_task_numfact_c
  end interface

  interface
     function pastix_task_solve_c(pastix_data, spm, nrhs, b, ldb) &
          bind(c, name='pastix_task_solve')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_task_solve_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
       integer(c_int), value :: nrhs
       type(c_ptr), value :: b
       integer(c_int), value :: ldb
     end function pastix_task_solve_c
  end interface

  interface
     function pastix_task_refine_c(pastix_data, x, rhsnbr, b) &
          bind(c, name='pastix_task_refine')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_task_refine_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: x
       integer(c_int), value :: rhsnbr
       type(c_ptr), value :: b
     end function pastix_task_refine_c
  end interface

  interface
     function pastix_subtask_order_c(pastix_data, spm, myorder) &
          bind(c, name='pastix_subtask_order')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_order_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
       type(c_ptr), value :: myorder
     end function pastix_subtask_order_c
  end interface

  interface
     function pastix_subtask_symbfact_c(pastix_data) &
          bind(c, name='pastix_subtask_symbfact')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_symbfact_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_symbfact_c
  end interface

  interface
     function pastix_subtask_reordering_c(pastix_data) &
          bind(c, name='pastix_subtask_reordering')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_reordering_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_reordering_c
  end interface

  interface
     function pastix_subtask_blend_c(pastix_data) &
          bind(c, name='pastix_subtask_blend')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_blend_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_blend_c
  end interface

  interface
     function pastix_subtask_spm2bcsc_c(pastix_data, spm) &
          bind(c, name='pastix_subtask_spm2bcsc')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_spm2bcsc_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_subtask_spm2bcsc_c
  end interface

  interface
     function pastix_subtask_bcsc2ctab_c(pastix_data, spm) &
          bind(c, name='pastix_subtask_bcsc2ctab')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_bcsc2ctab_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_subtask_bcsc2ctab_c
  end interface

  interface
     function pastix_subtask_sopalin_c(pastix_data, spm) &
          bind(c, name='pastix_subtask_sopalin')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_sopalin_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_subtask_sopalin_c
  end interface

  interface
     function pastix_subtask_applyorder_c(pastix_data, flttype, dir, m, n, b, ldb) &
          bind(c, name='pastix_subtask_applyorder')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_applyorder_c
       type(c_ptr), value :: pastix_data
       integer(c_int), value :: flttype
       integer(c_int), value :: dir
       integer(c_int), value :: m
       integer(c_int), value :: n
       type(c_ptr), value :: b
       integer(c_int), value :: ldb
     end function pastix_subtask_applyorder_c
  end interface

  interface
     function pastix_subtask_trsm_c(pastix_data, flttype, side, uplo, trans, diag, nrhs, b, &
          ldb) &
          bind(c, name='pastix_subtask_trsm')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_trsm_c
       type(c_ptr), value :: pastix_data
       integer(c_int), value :: flttype
       integer(c_int), value :: side
       integer(c_int), value :: uplo
       integer(c_int), value :: trans
       integer(c_int), value :: diag
       integer(c_int), value :: nrhs
       type(c_ptr), value :: b
       integer(c_int), value :: ldb
     end function pastix_subtask_trsm_c
  end interface

  interface
     function pastix_subtask_diag_c(pastix_data, flttype, nrhs, b, ldb) &
          bind(c, name='pastix_subtask_diag')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_diag_c
       type(c_ptr), value :: pastix_data
       integer(c_int), value :: flttype
       integer(c_int), value :: nrhs
       type(c_ptr), value :: b
       integer(c_int), value :: ldb
     end function pastix_subtask_diag_c
  end interface

  interface
     subroutine pastix_setSchurUnknownList_c(pastix_data, n, list) &
          bind(c, name='pastix_setSchurUnknownList')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pastix_data
       integer(c_int), value :: n
       type(c_ptr), value :: list
     end subroutine pastix_setSchurUnknownList_c
  end interface

  interface
     function pastix_getSchur_c(pastix_data, S, lds) &
          bind(c, name='pastix_getSchur')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_getSchur_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: S
       integer(c_int), value :: lds
     end function pastix_getSchur_c
  end interface

  interface
     subroutine pastix_getOptions_c(argc, argv, iparam, dparam, check, driver, filename) &
          bind(c, name='pastix_getOptions')
       use iso_c_binding
       implicit none
       integer(kind=c_int), value :: argc
       type(c_ptr) :: argv
       type(c_ptr), value :: iparam
       type(c_ptr), value :: dparam
       type(c_ptr), value :: check
       type(c_ptr), value :: driver
       type(c_ptr) :: filename
     end subroutine pastix_getOptions_c
  end interface

  interface
     subroutine spmInit_c(spm) &
          bind(c, name='spmInit')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmInit_c
  end interface

  interface
     subroutine spmExit_c(spm) &
          bind(c, name='spmExit')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmExit_c
  end interface

  interface
     function spmCopy_c(spm) &
          bind(c, name='spmCopy')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr) :: spmCopy_c
       type(c_ptr), value :: spm
     end function spmCopy_c
  end interface

  interface
     subroutine spmBase_c(spm, baseval) &
          bind(c, name='spmBase')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
       integer(kind=c_int), value :: baseval
     end subroutine spmBase_c
  end interface

  interface
     function spmFindBase_c(spm) &
          bind(c, name='spmFindBase')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(c_int) :: spmFindBase_c
       type(c_ptr), value :: spm
     end function spmFindBase_c
  end interface

  interface
     function spmConvert_c(ofmttype, ospm) &
          bind(c, name='spmConvert')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int) :: spmConvert_c
       integer(kind=c_int), value :: ofmttype
       type(c_ptr), value :: ospm
     end function spmConvert_c
  end interface

  interface
     subroutine spmUpdateComputedFields_c(spm) &
          bind(c, name='spmUpdateComputedFields')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmUpdateComputedFields_c
  end interface

  interface
     function spmNorm_c(ntype, spm) &
          bind(c, name='spmNorm')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       real(kind=c_double) :: spmNorm_c
       integer(c_int), value :: ntype
       type(c_ptr), value :: spm
     end function spmNorm_c
  end interface

  interface
     function spmMatVec_c(trans, alpha, spm, x, beta, y) &
          bind(c, name='spmMatVec')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int) :: spmMatVec_c
       integer(c_int), value :: trans
       type(c_ptr), value :: alpha
       type(c_ptr), value :: spm
       type(c_ptr), value :: x
       type(c_ptr), value :: beta
       type(c_ptr), value :: y
     end function spmMatVec_c
  end interface

  interface
     subroutine spmScalMatrix_c(alpha, spm) &
          bind(c, name='spmScalMatrix')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: spm
     end subroutine spmScalMatrix_c
  end interface

  interface
     subroutine spmScalVector_c(alpha, spm, x) &
          bind(c, name='spmScalVector')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       real(kind=c_double), value :: alpha
       type(c_ptr), value :: spm
       type(c_ptr), value :: x
     end subroutine spmScalVector_c
  end interface

  interface
     function spmSort_c(spm) &
          bind(c, name='spmSort')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int) :: spmSort_c
       type(c_ptr), value :: spm
     end function spmSort_c
  end interface

  interface
     function spmMergeDuplicate_c(spm) &
          bind(c, name='spmMergeDuplicate')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(c_int) :: spmMergeDuplicate_c
       type(c_ptr), value :: spm
     end function spmMergeDuplicate_c
  end interface

  interface
     function spmSymmetrize_c(spm) &
          bind(c, name='spmSymmetrize')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(c_int) :: spmSymmetrize_c
       type(c_ptr), value :: spm
     end function spmSymmetrize_c
  end interface

  interface
     function spmCheckAndCorrect_c(spm) &
          bind(c, name='spmCheckAndCorrect')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr) :: spmCheckAndCorrect_c
       type(c_ptr), value :: spm
     end function spmCheckAndCorrect_c
  end interface

  interface
     function spmGenRHS_c(type, nrhs, spm, x, ldx, b, ldb) &
          bind(c, name='spmGenRHS')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int) :: spmGenRHS_c
       integer(c_int), value :: type
       integer(kind=c_int), value :: nrhs
       type(c_ptr), value :: spm
       type(c_ptr), value :: x
       integer(kind=c_int), value :: ldx
       type(c_ptr), value :: b
       integer(kind=c_int), value :: ldb
     end function spmGenRHS_c
  end interface

  interface
     function spmCheckAxb_c(nrhs, spm, x0, ldx0, b, ldb, x, ldx) &
          bind(c, name='spmCheckAxb')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int) :: spmCheckAxb_c
       integer(kind=c_int), value :: nrhs
       type(c_ptr), value :: spm
       type(c_ptr), value :: x0
       integer(kind=c_int), value :: ldx0
       type(c_ptr), value :: b
       integer(kind=c_int), value :: ldb
       type(c_ptr), value :: x
       integer(kind=c_int), value :: ldx
     end function spmCheckAxb_c
  end interface

  interface
     function spmIntConvert_c(n, input) &
          bind(c, name='spmIntConvert')
       use iso_c_binding
       implicit none
       type(c_ptr) :: spmIntConvert_c
       integer(c_int), value :: n
       type(c_ptr), value :: input
     end function spmIntConvert_c
  end interface

  interface
     subroutine spmIntSort1Asc1_c(pbase, n) &
          bind(c, name='spmIntSort1Asc1')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pbase
       integer(c_int), value :: n
     end subroutine spmIntSort1Asc1_c
  end interface

  interface
     subroutine spmIntSort2Asc1_c(pbase, n) &
          bind(c, name='spmIntSort2Asc1')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pbase
       integer(c_int), value :: n
     end subroutine spmIntSort2Asc1_c
  end interface

  interface
     subroutine spmIntSort2Asc2_c(pbase, n) &
          bind(c, name='spmIntSort2Asc2')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pbase
       integer(c_int), value :: n
     end subroutine spmIntSort2Asc2_c
  end interface

  interface
     function spmReadDriver_c(driver, filename, spm, pastix_comm) &
          bind(c, name='spmReadDriver')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int) :: spmReadDriver_c
       integer(c_int), value :: driver
       type(c_ptr), value :: filename
       type(c_ptr), value :: spm
       integer(kind=c_int), value :: pastix_comm
     end function spmReadDriver_c
  end interface

  interface
     subroutine spm2Dense_c(spm) &
          bind(c, name='spm2Dense')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spm2Dense_c
  end interface

  interface
     function spmExpand_c(spm) &
          bind(c, name='spmExpand')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr) :: spmExpand_c
       type(c_ptr), value :: spm
     end function spmExpand_c
  end interface

  interface
     function spmDofExtend_c(spm, type, dof) &
          bind(c, name='spmDofExtend')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr) :: spmDofExtend_c
       type(c_ptr), value :: spm
       integer(kind=c_int), value :: type
       integer(kind=c_int), value :: dof
     end function spmDofExtend_c
  end interface

contains

  ! Wrappers of the C functions.
  subroutine pastix(pastix_data, pastix_comm, n, colptr, row, avals, perm, invp, &
       b, nrhs, iparm, dparm, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), pointer :: pastix_data
    integer(kind=c_int), intent(in) :: pastix_comm
    integer(c_int), intent(in) :: n
    integer(c_int), intent(inout), target :: colptr
    integer(c_int), intent(inout), target :: row
    type(c_ptr), intent(inout), target :: avals
    integer(c_int), intent(inout), target :: perm
    integer(c_int), intent(inout), target :: invp
    type(c_ptr), intent(inout), target :: b
    integer(c_int), intent(in) :: nrhs
    integer(c_int), intent(inout), target :: iparm
    real(kind=c_double), intent(inout), target :: dparm
    integer(kind=c_int), intent(out) :: info

    type(c_ptr) :: pastix_data_aux

    info = pastix_c(pastix_data_aux, pastix_comm, n, c_loc(colptr), c_loc(row), c_loc(avals), c_loc(perm), c_loc(invp), &
         c_loc(b), nrhs, c_loc(iparm), c_loc(dparm))
    call c_f_pointer(pastix_data_aux, pastix_data)
  end subroutine pastix

  subroutine pastixInitParam(iparm, dparm)
    use iso_c_binding
    implicit none
    integer(c_int), intent(inout), target :: iparm
    real(kind=c_double), intent(inout), target :: dparm

    call pastixInitParam_c(c_loc(iparm), c_loc(dparm))
  end subroutine pastixInitParam

  subroutine pastixInit(pastix_data, pastix_comm, iparm, dparm)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), pointer :: pastix_data
    integer(kind=c_int), intent(in) :: pastix_comm
    integer(c_int), intent(inout), target :: iparm
    real(kind=c_double), intent(inout), target :: dparm

    type(c_ptr) :: pastix_data_aux

    call pastixInit_c(pastix_data_aux, pastix_comm, c_loc(iparm), c_loc(dparm))
    call c_f_pointer(pastix_data_aux, pastix_data)
  end subroutine pastixInit

  subroutine pastixFinalize(pastix_data)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), pointer :: pastix_data

    type(c_ptr) :: pastix_data_aux

    call pastixFinalize_c(pastix_data_aux)
    call c_f_pointer(pastix_data_aux, pastix_data)
  end subroutine pastixFinalize

  subroutine pastix_task_analyze(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(out) :: info

    info = pastix_task_analyze_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_task_analyze

  subroutine pastix_task_numfact(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(out) :: info

    info = pastix_task_numfact_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_task_numfact

  subroutine pastix_task_solve(pastix_data, spm, nrhs, b, ldb, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(pastix_spm_t), intent(inout), target :: spm
    integer(c_int), intent(in) :: nrhs
    type(c_ptr), intent(inout), target :: b
    integer(c_int), intent(in) :: ldb
    integer(kind=c_int), intent(out) :: info

    info = pastix_task_solve_c(c_loc(pastix_data), c_loc(spm), nrhs, c_loc(b), ldb)
  end subroutine pastix_task_solve

  subroutine pastix_task_refine(pastix_data, x, rhsnbr, b, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(c_ptr), intent(inout), target :: x
    integer(c_int), intent(in) :: rhsnbr
    type(c_ptr), intent(inout), target :: b
    integer(kind=c_int), intent(out) :: info

    info = pastix_task_refine_c(c_loc(pastix_data), c_loc(x), rhsnbr, c_loc(b))
  end subroutine pastix_task_refine

  subroutine pastix_subtask_order(pastix_data, spm, myorder, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(pastix_spm_t), intent(inout), target :: spm
    type(c_ptr), intent(inout), target :: myorder
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_order_c(c_loc(pastix_data), c_loc(spm), c_loc(myorder))
  end subroutine pastix_subtask_order

  subroutine pastix_subtask_symbfact(pastix_data, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_symbfact_c(c_loc(pastix_data))
  end subroutine pastix_subtask_symbfact

  subroutine pastix_subtask_reordering(pastix_data, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_reordering_c(c_loc(pastix_data))
  end subroutine pastix_subtask_reordering

  subroutine pastix_subtask_blend(pastix_data, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_blend_c(c_loc(pastix_data))
  end subroutine pastix_subtask_blend

  subroutine pastix_subtask_spm2bcsc(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_spm2bcsc_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_subtask_spm2bcsc

  subroutine pastix_subtask_bcsc2ctab(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_bcsc2ctab_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_subtask_bcsc2ctab

  subroutine pastix_subtask_sopalin(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_sopalin_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_subtask_sopalin

  subroutine pastix_subtask_applyorder(pastix_data, flttype, dir, m, n, b, ldb, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    integer(c_int), intent(in) :: flttype
    integer(c_int), intent(in) :: dir
    integer(c_int), intent(in) :: m
    integer(c_int), intent(in) :: n
    type(c_ptr), intent(inout), target :: b
    integer(c_int), intent(in) :: ldb
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_applyorder_c(c_loc(pastix_data), flttype, dir, m, n, c_loc(b), ldb)
  end subroutine pastix_subtask_applyorder

  subroutine pastix_subtask_trsm(pastix_data, flttype, side, uplo, trans, diag, nrhs, b, &
       ldb, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    integer(c_int), intent(in) :: flttype
    integer(c_int), intent(in) :: side
    integer(c_int), intent(in) :: uplo
    integer(c_int), intent(in) :: trans
    integer(c_int), intent(in) :: diag
    integer(c_int), intent(in) :: nrhs
    type(c_ptr), intent(inout), target :: b
    integer(c_int), intent(in) :: ldb
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_trsm_c(c_loc(pastix_data), flttype, side, uplo, trans, diag, nrhs, c_loc(b), &
         ldb)
  end subroutine pastix_subtask_trsm

  subroutine pastix_subtask_diag(pastix_data, flttype, nrhs, b, ldb, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    integer(c_int), intent(in) :: flttype
    integer(c_int), intent(in) :: nrhs
    type(c_ptr), intent(inout), target :: b
    integer(c_int), intent(in) :: ldb
    integer(kind=c_int), intent(out) :: info

    info = pastix_subtask_diag_c(c_loc(pastix_data), flttype, nrhs, c_loc(b), ldb)
  end subroutine pastix_subtask_diag

  subroutine pastix_setSchurUnknownList(pastix_data, n, list)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    integer(c_int), intent(in) :: n
    integer(c_int), intent(inout), target :: list

    call pastix_setSchurUnknownList_c(c_loc(pastix_data), n, c_loc(list))
  end subroutine pastix_setSchurUnknownList

  subroutine pastix_getSchur(pastix_data, S, lds, info)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pastix_data
    type(c_ptr), intent(inout), target :: S
    integer(c_int), intent(in) :: lds
    integer(kind=c_int), intent(out) :: info

    info = pastix_getSchur_c(c_loc(pastix_data), c_loc(S), lds)
  end subroutine pastix_getSchur

  subroutine pastix_getOptions(argc, argv, iparam, dparam, check, driver, filename)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(in) :: argc
    character(kind=c_char), intent(inout), pointer :: argv
    integer(c_int), intent(inout), target :: iparam
    real(kind=c_double), intent(inout), target :: dparam
    integer(kind=c_int), intent(inout), target :: check
    integer(c_int), intent(inout), target :: driver
    character(kind=c_char), intent(inout), pointer :: filename

    type(c_ptr) :: argv_aux

    type(c_ptr) :: filename_aux

    call pastix_getOptions_c(argc, argv_aux, c_loc(iparam), c_loc(dparam), c_loc(check), c_loc(driver), filename_aux)
    call c_f_pointer(argv_aux, argv)
    call c_f_pointer(filename_aux, filename)
  end subroutine pastix_getOptions

  subroutine spmInit(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmInit_c(c_loc(spm))
  end subroutine spmInit

  subroutine spmExit(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmExit_c(c_loc(spm))
  end subroutine spmExit

  subroutine spmCopy(spm, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    type(pastix_spm_t), intent(out), pointer :: spmo

    call c_f_pointer(spmCopy_c(c_loc(spm)), spmo)
  end subroutine spmCopy

  subroutine spmBase(spm, baseval)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(in) :: baseval

    call spmBase_c(c_loc(spm), baseval)
  end subroutine spmBase

  subroutine spmFindBase(spm, value)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    integer(c_int), intent(out) :: value

    value = spmFindBase_c(c_loc(spm))
  end subroutine spmFindBase

  subroutine spmConvert(ofmttype, ospm, info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(in) :: ofmttype
    type(pastix_spm_t), intent(inout), target :: ospm
    integer(kind=c_int), intent(out) :: info

    info = spmConvert_c(ofmttype, c_loc(ospm))
  end subroutine spmConvert

  subroutine spmUpdateComputedFields(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmUpdateComputedFields_c(c_loc(spm))
  end subroutine spmUpdateComputedFields

  subroutine spmNorm(ntype, spm, value)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: ntype
    type(pastix_spm_t), intent(inout), target :: spm
    real(kind=c_double), intent(out) :: value

    value = spmNorm_c(ntype, c_loc(spm))
  end subroutine spmNorm

  subroutine spmMatVec(trans, alpha, spm, x, beta, y, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: trans
    type(c_ptr), intent(inout), target :: alpha
    type(pastix_spm_t), intent(inout), target :: spm
    type(c_ptr), intent(inout), target :: x
    type(c_ptr), intent(inout), target :: beta
    type(c_ptr), intent(inout), target :: y
    integer(kind=c_int), intent(out) :: info

    info = spmMatVec_c(trans, c_loc(alpha), c_loc(spm), c_loc(x), c_loc(beta), c_loc(y))
  end subroutine spmMatVec

  subroutine spmScalMatrix(alpha, spm)
    use iso_c_binding
    implicit none
    complex(kind=c_double_complex), intent(in) :: alpha
    type(pastix_spm_t), intent(inout), target :: spm

    call spmScalMatrix_c(alpha, c_loc(spm))
  end subroutine spmScalMatrix

  subroutine spmScalVector(alpha, spm, x)
    use iso_c_binding
    implicit none
    real(kind=c_double), intent(in) :: alpha
    type(pastix_spm_t), intent(inout), target :: spm
    type(c_ptr), intent(inout), target :: x

    call spmScalVector_c(alpha, c_loc(spm), c_loc(x))
  end subroutine spmScalVector

  subroutine spmSort(spm, info)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(out) :: info

    info = spmSort_c(c_loc(spm))
  end subroutine spmSort

  subroutine spmMergeDuplicate(spm, value)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    integer(c_int), intent(out) :: value

    value = spmMergeDuplicate_c(c_loc(spm))
  end subroutine spmMergeDuplicate

  subroutine spmSymmetrize(spm, value)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    integer(c_int), intent(out) :: value

    value = spmSymmetrize_c(c_loc(spm))
  end subroutine spmSymmetrize

  subroutine spmCheckAndCorrect(spm, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    type(pastix_spm_t), intent(out), pointer :: spmo

    call c_f_pointer(spmCheckAndCorrect_c(c_loc(spm)), spmo)
  end subroutine spmCheckAndCorrect

  subroutine spmGenRHS(type, nrhs, spm, x, ldx, b, ldb, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: type
    integer(kind=c_int), intent(in) :: nrhs
    type(pastix_spm_t), intent(inout), target :: spm
    type(c_ptr), intent(inout), target :: x
    integer(kind=c_int), intent(in) :: ldx
    type(c_ptr), intent(inout), target :: b
    integer(kind=c_int), intent(in) :: ldb
    integer(kind=c_int), intent(out) :: info

    info = spmGenRHS_c(type, nrhs, c_loc(spm), c_loc(x), ldx, c_loc(b), ldb)
  end subroutine spmGenRHS

  subroutine spmCheckAxb(nrhs, spm, x0, ldx0, b, ldb, x, ldx, info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(in) :: nrhs
    type(pastix_spm_t), intent(inout), target :: spm
    type(c_ptr), intent(inout), target :: x0
    integer(kind=c_int), intent(in) :: ldx0
    type(c_ptr), intent(inout), target :: b
    integer(kind=c_int), intent(in) :: ldb
    type(c_ptr), intent(inout), target :: x
    integer(kind=c_int), intent(in) :: ldx
    integer(kind=c_int), intent(out) :: info

    info = spmCheckAxb_c(nrhs, c_loc(spm), c_loc(x0), ldx0, c_loc(b), ldb, c_loc(x), ldx)
  end subroutine spmCheckAxb

  subroutine spmIntConvert(n, input, value)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n
    integer(kind=c_int), intent(inout), target :: input
    integer(c_int), intent(out), pointer :: value

    call c_f_pointer(spmIntConvert_c(n, c_loc(input)), value)
  end subroutine spmIntConvert

  subroutine spmIntSort1Asc1(pbase, n)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pbase
    integer(c_int), intent(in) :: n

    call spmIntSort1Asc1_c(c_loc(pbase), n)
  end subroutine spmIntSort1Asc1

  subroutine spmIntSort2Asc1(pbase, n)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pbase
    integer(c_int), intent(in) :: n

    call spmIntSort2Asc1_c(c_loc(pbase), n)
  end subroutine spmIntSort2Asc1

  subroutine spmIntSort2Asc2(pbase, n)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout), target :: pbase
    integer(c_int), intent(in) :: n

    call spmIntSort2Asc2_c(c_loc(pbase), n)
  end subroutine spmIntSort2Asc2

  subroutine spmReadDriver(driver, filename, spm, pastix_comm, info)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: driver
    character(kind=c_char), intent(inout), target :: filename
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(in) :: pastix_comm
    integer(kind=c_int), intent(out) :: info

    info = spmReadDriver_c(driver, c_loc(filename), c_loc(spm), pastix_comm)
  end subroutine spmReadDriver

  subroutine spm2Dense(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spm2Dense_c(c_loc(spm))
  end subroutine spm2Dense

  subroutine spmExpand(spm, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    type(pastix_spm_t), intent(out), pointer :: spmo

    call c_f_pointer(spmExpand_c(c_loc(spm)), spmo)
  end subroutine spmExpand

  subroutine spmDofExtend(spm, type, dof, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm
    integer(kind=c_int), intent(in) :: type
    integer(kind=c_int), intent(in) :: dof
    type(pastix_spm_t), intent(out), pointer :: spmo

    call c_f_pointer(spmDofExtend_c(c_loc(spm), type, dof), spmo)
  end subroutine spmDofExtend
  
end module pastix_mod
