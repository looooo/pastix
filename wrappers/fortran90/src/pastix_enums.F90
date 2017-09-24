module pastix_enums
  use iso_c_binding
  implicit none

  ! iparm enum
  integer, parameter :: iparm_verbose               = 1
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
  integer, parameter :: iparm_2dtasks_level         = iparm_max_blocksize         + 1
  integer, parameter :: iparm_2dtasks_width         = iparm_2dtasks_level         + 1
  integer, parameter :: iparm_abs                   = iparm_2dtasks_width         + 1
  integer, parameter :: iparm_incomplete            = iparm_abs                   + 1
  integer, parameter :: iparm_level_of_fill         = iparm_incomplete            + 1
  integer, parameter :: iparm_factorization         = iparm_level_of_fill         + 1
  integer, parameter :: iparm_static_pivoting       = iparm_factorization         + 1
  integer, parameter :: iparm_inertia               = iparm_static_pivoting       + 1
  integer, parameter :: iparm_free_cscuser          = iparm_inertia               + 1
  integer, parameter :: iparm_schur_fact_mode       = iparm_free_cscuser          + 1
  integer, parameter :: iparm_schur_solv_mode       = iparm_schur_fact_mode       + 1
  integer, parameter :: iparm_refinement            = iparm_schur_solv_mode       + 1
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
  integer, parameter :: dparm_fill_in            = 1
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

  ! class factmode:
  integer, parameter :: PastixFactModeLocal = 0
  integer, parameter :: PastixFactModeSchur = 1
  integer, parameter :: PastixFactModeBoth  = 2

  ! class solvmode:
  integer, parameter :: PastixSolvModeLocal     = 0
  integer, parameter :: PastixSolvModeInterface = 1
  integer, parameter :: PastixSolvModeSchur     = 2

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
  integer, parameter :: PastixRhsOne  = 0
  integer, parameter :: PastixRhsI    = 1
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
  integer, parameter :: pastix_int_t = PASTIX_INT_KIND

contains

  function pastix_getintsize()
    integer :: pastix_getintsize
    pastix_getintsize = PASTIX_INT_KIND
    return
  end function pastix_getintsize

end module pastix_enums
