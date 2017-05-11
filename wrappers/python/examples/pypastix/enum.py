import numpy as np
import ctypes

pastix_int = ctypes.c_int64

class iparm:
    factorization = 44
    size = 128

class dparm:
    size = 128

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
