"""
 @file spm.py

 SPM python wrapper

 @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @date 2017-05-04

"""
from ctypes import *
import numpy as np
import scipy.sparse as sps

from .enum    import *
from .__spm__ import *

class spm():

    dtype = None

    def __init__( self, A=None, mtxtype_=mtxtype.General, driver=None, filename="" ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        if mtxtype_ == mtxtype.SymPosDef:
            mtxtype_ = mtxtype.Symmetric
        if mtxtype_ == mtxtype.HerPosDef:
            mtxtype_ = mtxtype.Hermitian

        self.spm_c = pypastix_spm_t( mtxtype_,
                                     coeftype.Double,
                                     fmttype.CSC,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     1, None,
                                     layout.ColMajor,
                                     None, None, None, None )
        self.id_ptr = pointer( self.spm_c )
        if A is not None:
            self.fromsps( A, mtxtype_ )
        elif driver is not None:
            self.fromdriver( driver, filename )

    def fromsps( self, A, mtxtype_=mtxtype.General ):
        """
        Initialize the SPM wrapper by loading the libraries
        """

        if not sps.isspmatrix(A):
            raise TypeError( "A must be of type scipy.sparse matrix" )

        # Assume A is already in Scipy sparse format
        self.dtype = A.dtype
        flttype = coeftype.getptype( A.dtype )
        if flttype == -1:
            raise TypeError( "Invalid data type. Must be part of (f4, f8, c8 or c16)" )

        A = sps.csc_matrix( A )
        A.sort_indices()

        # Pointer variables
        self.py_colptr = np.array( A.indptr[:],  dtype=pastix_int )
        self.py_rowptr = np.array( A.indices[:], dtype=pastix_int )
        self.py_values = np.array( A.data[:] )

        if mtxtype_ == mtxtype.SymPosDef:
            mtxtype_ = mtxtype.Symmetric
        if mtxtype_ == mtxtype.HerPosDef:
            mtxtype_ = mtxtype.Hermitian

        self.spm_c.mtxtype  = mtxtype_
        self.spm_c.flttype  = flttype
        self.spm_c.fmttype  = fmttype.CSC
        self.spm_c.n        = A.shape[0]
        self.spm_c.nnz      = A.getnnz()
        self.spm_c.dof      = 1
        self.spm_c.dofs     = None
        self.spm_c.layout   = layout.ColMajor
        self.spm_c.colptr   = self.py_colptr.ctypes.data_as(POINTER(pastix_int))
        self.spm_c.rowptr   = self.py_rowptr.ctypes.data_as(POINTER(pastix_int))
        self.spm_c.loc2glob = None
        self.spm_c.values   = self.py_values.ctypes.data_as(c_void_p)

        self.id_ptr = pointer( self.spm_c )

        self.updateComputedField()
        self.checkAndCorrect()

    def tosps( self ):
        """
        Return a Scipy sparse matrix
        """
        n      = int( self.spm_c.n )
        nnz    = int( self.spm_c.nnz )
        cflt   = coeftype.getctype(  self.spm_c.flttype )
        nflt   = coeftype.getnptype( self.spm_c.flttype )
        colptr = self.py_colptr.copy()
        rowptr = self.py_rowptr.copy()
        values = self.py_values.copy()

        baseval = colptr[0]
        colptr = colptr - baseval
        rowptr = rowptr - baseval

        return sps.csc_matrix((values, rowptr, colptr), shape=(n, n))

    def fromdriver( self, driver=driver.Laplacian, filename="10:10:10" ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        if filename == "":
            raise ValueError("filename must be prodived")

        pyspm_spmReadDriver( driver, filename.encode('utf-8'), self.id_ptr, c_int(0) )

        self.dtype = coeftype.getnptype( self.spm_c.flttype )
        self.checkAndCorrect()

    def scale( self, alpha ):
        return pyspm_spmScalMatrix( float(alpha), self.id_ptr )

    def norm( self, ntype=normtype.Frobenius ):
        return pyspm_spmNorm( ntype, self.id_ptr )

    def base( self, baseval ):
        pyspm_spmBase( self.id_ptr, baseval )

    def findBase( self ):
        return pyspm_spmFindBase( self.id_ptr )

    def updateComputedField( self ):
        pyspm_spmUpdateComputedFields( self.id_ptr )

    def printInfo( self ):
        pyspm_spmPrintInfo( self.id_ptr )

    def printSpm( self ):
        pyspm_spmPrint( self.id_ptr )

    def checkAndCorrect( self ):
        spm1 = self.id_ptr
        spm2 = pyspm_spmCheckAndCorrect( self.id_ptr )
        if (( spm1.contents.fmttype == spm2.contents.fmttype ) and
            ( spm1.contents.nnzexp  == spm2.contents.nnzexp  ) ):
            return
        self.spm_c = cast( spm2, POINTER(pypastix_spm_t) ).contents
        self.id_ptr = pointer( self.spm_c )

        self.py_colptr = np.frombuffer( (pastix_int * (n+1)).from_address( cast(self.spm_c.colptr, c_void_p).value ), pastix_int ).copy()
        self.py_rowptr = np.frombuffer( (pastix_int *  nnz ).from_address( cast(self.spm_c.rowptr, c_void_p).value ), pastix_int ).copy()
        self.py_values = np.frombuffer( (cflt       *  nnz ).from_address( self.spm_c.values ), nflt ).copy()

    def __checkVector( self, n, nrhs, x ):
        if x.dtype != self.dtype:
            raise TypeError( "Vectors must use the same arithmetic as the spm" )

        if x.ndim > 2:
            raise TypeError( "Vectors must be of dimension 1 or 2" )

        if x.shape[0] < n:
            raise TypeError( "Vectors must be of size at least ", n )

        if (x.ndim == 1 and nrhs > 1) or (x.ndim>1 and x.shape[1] < nrhs):
            raise TypeError( "At least nrhs vectors must be stored in the vector" )

    def checkAxb( self, x0, b, x, nrhs=-1 ):
        # if libspm == None:
        #     raise EnvironmentError( "SPM Instance badly instanciated" )

        b = np.array(b, self.dtype)
        x = np.asarray(x, self.dtype)

        n = self.spm_c.n
        if nrhs == -1:
            if x.ndim == 1:
                nrhs = 1
            else:
                nrhs = x.shape[1]

        if x0 is not None:
            x0 = np.asarray(x0, self.dtype)
            self.__checkVector( n, nrhs, x0 )
            ldx0 = x0.shape[0]
            x0ptr = x0.ctypes.data_as(c_void_p)
        else:
            ldx0 = 1
            x0ptr = None
        self.__checkVector( n, nrhs, b )
        self.__checkVector( n, nrhs, x )

        ldb  = b.shape[0]
        ldx  = x.shape[0]

        pyspm_spmCheckAxb( nrhs, self.id_ptr,
                           x0ptr, ldx0,
                           b.ctypes.data_as(c_void_p), ldb,
                           x.ctypes.data_as(c_void_p), ldx )

    def genRHS( self, rhstype=rhstype.One, nrhs=1 ):
        # if libspm == None:
        #     raise EnvironmentError( "SPM Instance badly instanciated" )

        n = self.spm_c.n
        b = np.zeros(n, self.dtype)
        x = np.zeros(n, self.dtype)

        ldb = b.shape[0]
        ldx = x.shape[0]

        self.__checkVector( n, nrhs, x )
        self.__checkVector( n, nrhs, b )

        pyspm_spmGenRHS( rhstype, nrhs, self.id_ptr,
                         x.ctypes.data_as(c_void_p), ldx,
                         b.ctypes.data_as(c_void_p), ldb )
        return x, b
