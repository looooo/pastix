import numpy as np
from ctypes import *

pastix_c_int  = c_int64
pastix_np_int = np.dtype('int64')

class pastix_iparm:
    iparm_factorization = 35 # DANGER !!!
    iparm_size = 128

class pastix_dparm:
    dparm_size = 128

class pastix_task:
    PastixTaskInit       = 0 # Startup the library
    PastixTaskOrdering   = 1 # Ordering
    PastixTaskSymbfact   = 2 # Symbolic factorization
    PastixTaskAnalyze    = 3 # Tasks mapping and scheduling
    PastixTaskNumfact    = 4 # Numerical factorization
    PastixTaskSolve      = 5 # Numerical solve
    PastixTaskRefine     = 6 # Numerical refinement
    PastixTaskClean      = 7 # Clean


class pastix_verbose:
    PastixVerboseNot        = 0 # Nothing
    PastixVerboseNo         = 1 # Default
    PastixVerboseYes        = 2 # Extended

class pastix_io:
    PastixIONo         = 0
    PastixIOLoad       = 1
    PastixIOSave       = 2
    PastixIOLoadGraph  = 4
    PastixIOSaveGraph  = 8
    PastixIOLoadCSC    = 16
    PastixIOSaveCSC    = 32

class pastix_refine:
    PastixRefineGMRES    = 0 # GMRES
    PastixRefineCG       = 1 # CG
    PastixRefineSR       = 2 # SR
    PastixRefineBiCGSTAB = 3 # BiCGstab

class pastix_coeftype:
    PastixPattern   = 0 # Pattern only, no values are stored
    PastixFloat     = 2 # Single precision real
    PastixDouble    = 3 # Double precision real
    PastixComplex32 = 4 # Single precision complex
    PastixComplex64 = 5 # Double precision complex

    @staticmethod
    def get ( dtype ):
        np_dict = {
            np.dtype('float32')    : pastix_coeftype.PastixFloat,
            np.dtype('float64')    : pastix_coeftype.PastixDouble,
            np.dtype('complex64')  : pastix_coeftype.PastixComplex32,
            np.dtype('complex128') : pastix_coeftype.PastixComplex64,
        }
        if dtype in np_dict:
            return np_dict[dtype]
        else:
            return -1

    @staticmethod
    def getdtype ( flttype ):
        np_dict = {
            pastix_coeftype.PastixFloat     : np.dtype('float32'),
            pastix_coeftype.PastixDouble    : np.dtype('float64'),
            pastix_coeftype.PastixComplex32 : np.dtype('complex64'),
            pastix_coeftype.PastixComplex64 : np.dtype('complex128')
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1

class pastix_fmttype:
    PastixCSC = 0 # Compressed sparse column
    PastixCSR = 1 # Compressed sparse row
    PastixIJV = 2 # Coordinates

class pastix_factotype:
    PastixFactLLT  = 0 # Cholesky factorization
    PastixFactLDLT = 1 # LDL^t factorization
    PastixFactLU   = 2 # LU factorization
    PastixFactLDLH = 3 # LDL^h factorization for complex matrices

class pastix_rhstype:
    PastixRhsOne =  0
    PastixRhsI   =  1
    PastixRhsRndX = 2
    PastixRhsRndB = 3

class pastix_normtype:
    PastixOneNorm       = 171 # One norm:       max_j( sum_i( |a_{ij}| ) )
    PastixFrobeniusNorm = 174 # Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) )
    PastixInfNorm       = 175 # Inifinite norm: max_i( sum_j( |a_{ij}| ) )
    PastixMaxNorm       = 177 # Inifinite norm: max_{i,j}( | a_{ij} | )

class pastix_order:
    PastixRowMajor  = 101 # Storage in row major order
    PastixColMajor  = 102 # Storage in column major order

class pastix_trans:
    PastixNoTrans   = 111 # Use A
    PastixTrans     = 112 # Use A^t
    PastixConjTrans = 113 # Use conj(A^t)

class pastix_uplo:
    PastixUpper      = 121 # Use lower triangle of A
    PastixLower      = 122 # Use upper triangle of A
    PastixUpperLower = 123 # Use the full A

class pastix_diag:
    PastixNonUnit = 131 # Diagonal is non unitary
    PastixUnit    = 132 # Diagonal is unitary

class pastix_side:
    PastixLeft  = 141 # Apply operator on the left
    PastixRight = 142 # Apply operator on the right

class pastix_coefside:
    PastixLCoef = 0 # Coefficients of the lower triangular L are used
    PastixUCoef = 1 # Coefficients of the upper triangular U are used

class pastix_dir:
    PastixDirForward  = 391 # Forward direction
    PastixDirBackward = 392 # Backward direction

class pastix_mtxtype:
    PastixGeneral   = pastix_trans.PastixNoTrans    # The matrix is general
    PastixSymmetric = pastix_trans.PastixTrans      # The matrix is symmetric
    PastixHermitian = pastix_trans.PastixConjTrans  # The matrix is hermitian

class pastix_driver:
    PastixDriverRSA        = 0  # RSA driver
    PastixDriverHB         = 1  # Harwell Boeing driver
    PastixDriverIJV        = 2  # IJV Coordinate driver
    PastixDriverMM         = 3  # Matrix Market driver
    PastixDriverLaplacian  = 4  # 3, 5, or 7 points Lapalacian stencil generator
    PastixDriverXLaplacian = 5  # 15-points Laplacian stencil generator
    PastixDriverGraph      = 6  # Scotch Graph driver
