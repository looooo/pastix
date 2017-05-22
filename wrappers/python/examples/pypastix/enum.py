import numpy as np
import ctypes

pastix_int = ctypes.c_int64

class iparm:
    verbose               = 0
    io_strategy           = verbose                + 1
    nnzeros               = io_strategy            + 1
    nnzeros_block_local   = nnzeros                + 1
    allocated_terms       = nnzeros_block_local    + 1
    produce_stats         = allocated_terms        + 1
    mc64                  = produce_stats          + 1
    ordering              = mc64                   + 1
    ordering_default      = ordering               + 1
    scotch_switch_level   = ordering_default       + 1
    scotch_cmin           = scotch_switch_level    + 1
    scotch_cmax           = scotch_cmin            + 1
    scotch_frat           = scotch_cmax            + 1
    metis_ctype           = scotch_frat            + 1
    metis_rtype           = metis_ctype            + 1
    metis_no2hop          = metis_rtype            + 1
    metis_nseps           = metis_no2hop           + 1
    metis_niter           = metis_nseps            + 1
    metis_ufactor         = metis_niter            + 1
    metis_compress        = metis_ufactor          + 1
    metis_ccorder         = metis_compress         + 1
    metis_pfactor         = metis_ccorder          + 1
    metis_seed            = metis_pfactor          + 1
    metis_dbglvl          = metis_seed             + 1
    sf_kass               = metis_dbglvl           + 1
    amalgamation_lvlblas  = sf_kass                + 1
    amalgamation_lvlcblk  = amalgamation_lvlblas   + 1
    reordering_split      = amalgamation_lvlcblk   + 1
    reordering_stop       = reordering_split       + 1
    min_blocksize         = reordering_stop        + 1
    max_blocksize         = min_blocksize          + 1
    distribution_level    = max_blocksize          + 1
    abs                   = distribution_level     + 1
    incomplete            = abs                    + 1
    level_of_fill         = incomplete             + 1
    factorization         = level_of_fill          + 1
    static_pivoting       = factorization          + 1
    inertia               = static_pivoting        + 1
    free_cscuser          = inertia                + 1
    refinement            = free_cscuser           + 1
    nbiter                = refinement             + 1
    itermax               = nbiter                 + 1
    gmres_im              = itermax                + 1
    scheduler             = gmres_im               + 1
    thread_nbr            = scheduler              + 1
    autosplit_comm        = thread_nbr             + 1
    gpu_nbr               = autosplit_comm         + 1
    gpu_memory_percentage = gpu_nbr                + 1
    gpu_memory_block_size = gpu_memory_percentage  + 1
    compress_min_width    = gpu_memory_block_size  + 1
    compress_min_height   = compress_min_width     + 1
    compress_when         = compress_min_height    + 1
    compress_method       = compress_when          + 1
    thread_comm_mode      = compress_method        + 1
    modify_parameter      = thread_comm_mode       + 1
    start_task            = modify_parameter       + 1
    end_task              = start_task             + 1
    baseval               = end_task               + 1
    float_type            = baseval                + 1
    mtx_type              = float_type             + 1
    dof_nbr               = mtx_type               + 1
    size                  = dof_nbr                + 1

class dparm:
    fill_in            = 0
    epsilon_refinement = fill_in            + 1
    relative_error     = epsilon_refinement + 1
    epsilon_magn_ctrl  = relative_error     + 1
    analyze_time       = epsilon_magn_ctrl  + 1
    pred_fact_time     = analyze_time       + 1
    fact_time          = pred_fact_time     + 1
    solv_time          = fact_time          + 1
    fact_flops         = solv_time          + 1
    fact_thflops       = fact_flops         + 1
    fact_rlflops       = fact_thflops       + 1
    solv_flops         = fact_rlflops       + 1
    solv_thflops       = solv_flops         + 1
    solv_rlflops       = solv_thflops       + 1
    refine_time        = solv_rlflops       + 1
    a_norm             = refine_time        + 1
    compress_tolerance = a_norm             + 1
    size               = compress_tolerance + 1

class task:
    Init       = 0 # Startup the library
    Ordering   = 1 # Ordering
    Symbfact   = 2 # Symbolic factorization
    Analyze    = 3 # Tasks mapping and scheduling
    Numfact    = 4 # Numerical factorization
    Solve      = 5 # Numerical solve
    Refine     = 6 # Numerical refinement
    Clean      = 7 # Clean


class verbose:
    Not        = 0 # Nothing
    No         = 1 # Default
    Yes        = 2 # Extended

class io:
    No         = 0
    Load       = 1
    Save       = 2
    LoadGraph  = 4
    SaveGraph  = 8
    LoadCSC    = 16
    SaveCSC    = 32

class refine:
    GMRES    = 0 # GMRES
    CG       = 1 # CG
    SR       = 2 # SR
    BiCGSTAB = 3 # BiCGstab

class coeftype:
    Pattern   = 0 # Pattern only, no values are stored
    Float     = 2 # Single precision real
    Double    = 3 # Double precision real
    Complex32 = 4 # Single precision complex
    Complex64 = 5 # Double precision complex

    @staticmethod
    def getptype ( dtype ):
        np_dict = {
            np.dtype('float32')    : coeftype.Float,
            np.dtype('float64')    : coeftype.Double,
            np.dtype('complex64')  : coeftype.Complex32,
            np.dtype('complex128') : coeftype.Complex64,
        }
        if dtype in np_dict:
            return np_dict[dtype]
        else:
            return -1

    @staticmethod
    def getnptype ( flttype ):
        np_dict = {
            coeftype.Float     : np.dtype('float32'),
            coeftype.Double    : np.dtype('float64'),
            coeftype.Complex32 : np.dtype('complex64'),
            coeftype.Complex64 : np.dtype('complex128')
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1

class fmttype:
    CSC = 0 # Compressed sparse column
    CSR = 1 # Compressed sparse row
    IJV = 2 # Coordinates

class factotype:
    LLT  = 0 # Cholesky factorization
    LDLT = 1 # LDL^t factorization
    LU   = 2 # LU factorization
    LDLH = 3 # LDL^h factorization for complex matrices

class rhstype:
    One =  0
    I   =  1
    RndX = 2
    RndB = 3

class normtype:
    One       = 171 # One norm:       max_j( sum_i( |a_{ij}| ) )
    Frobenius = 174 # Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) )
    Inf       = 175 # Inifinite norm: max_i( sum_j( |a_{ij}| ) )
    Max       = 177 # Inifinite norm: max_{i,j}( | a_{ij} | )

class order:
    RowMajor  = 101 # Storage in row major order
    ColMajor  = 102 # Storage in column major order

class trans:
    NoTrans   = 111 # Use A
    Trans     = 112 # Use A^t
    ConjTrans = 113 # Use conj(A^t)

class uplo:
    Upper      = 121 # Use lower triangle of A
    Lower      = 122 # Use upper triangle of A
    UpperLower = 123 # Use the full A

class diag:
    NonUnit = 131 # Diagonal is non unitary
    Unit    = 132 # Diagonal is unitary

class side:
    Left  = 141 # Apply operator on the left
    Right = 142 # Apply operator on the right

class coefside:
    LCoef = 0 # Coefficients of the lower triangular L are used
    UCoef = 1 # Coefficients of the upper triangular U are used

class direction:
    Forward  = 391 # Forward direction
    Backward = 392 # Backward direction

class mtxtype:
    General   = trans.NoTrans    # The matrix is general
    Symmetric = trans.Trans      # The matrix is symmetric
    Hermitian = trans.ConjTrans  # The matrix is hermitian

class driver:
    RSA        = 0  # RSA driver
    HB         = 1  # Harwell Boeing driver
    IJV        = 2  # IJV Coordinate driver
    MM         = 3  # Matrix Market driver
    Laplacian  = 4  # 3, 5, or 7 points Lapalacian stencil generator
    XLaplacian = 5  # 15-points Laplacian stencil generator
    Graph      = 6  # Scotch Graph driver
