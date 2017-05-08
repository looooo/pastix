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
import scipy.sparse as spp

from . import libspm
from .enum import *

class spm():

    dtype = None

    class c_spm(Structure):
        _fields_ = [("mtxtype", c_int                ),
                    ("flttype", c_int                ),
                    ("fmttype", c_int                ),
                    ("gN",      pastix_c_int         ),
                    ("n",       pastix_c_int         ),
                    ("gnnz",    pastix_c_int         ),
                    ("nnz",     pastix_c_int         ),
                    ("gNexp",   pastix_c_int         ),
                    ("nexp",    pastix_c_int         ),
                    ("gnnzexp", pastix_c_int         ),
                    ("nnzexp",  pastix_c_int         ),
                    ("dof",     pastix_c_int         ),
                    ("dofs",    POINTER(pastix_c_int)),
                    ("layout",  c_int                ),
                    ("colptr",  POINTER(pastix_c_int)),
                    ("rowptr",  POINTER(pastix_c_int)),
                    ("loc2glob",POINTER(pastix_c_int)),
                    ("values",  c_void_p             ) ]

    def __init__( self ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        self.spm_c = self.c_spm( pastix_mtxtype.PastixGeneral,
                                 pastix_coeftype.PastixDouble,
                                 pastix_fmttype.PastixCSC,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, None,
                                 pastix_order.PastixColMajor,
                                 None, None, None, None)
        self.id_ptr = pointer( self.spm_c )

    def fromspp( self, A, mtxtype=pastix_mtxtype.PastixGeneral ):
        """
        Initialize the SPM wrapper by loading the libraries
        """

        # Assume A is already in Scipy sparse format
        self.dtype = A.dtype
        flttype = pastix_coeftype.get( A.dtype )
        print( "Floating point arithmetic is", flttype )
        if flttype == -1:
            raise TypeError( "Invalid data type. Must be part of (f4, f8, c8 or c16)" )

        A = spp.csc_matrix( A )
        A.sort_indices()

        # Pointer variables
        self.py_colptr    = np.array( A.indptr[:],  dtype=pastix_np_int )
        self.py_rowptr    = np.array( A.indices[:], dtype=pastix_np_int )
        self.py_values    = np.array( A.data[:] )

        self.spm_c.mtxtype= mtxtype
        self.spm_c.flttype= flttype
        self.spm_c.n      = A.shape[0]
        self.spm_c.n      = A.shape[0]
        self.spm_c.nnz    = A.getnnz()
        self.spm_c.colptr = self.py_colptr.ctypes.data_as(POINTER(pastix_c_int))
        self.spm_c.rowptr = self.py_rowptr.ctypes.data_as(POINTER(pastix_c_int))
        self.spm_c.values = self.py_values.ctypes.data_as(c_void_p)

        self.id_ptr = pointer( self.spm_c )

        self.updateComputedField()
        self.checkAndCorrect()

    def fromdriver( self, driver=pastix_driver.PastixDriverLaplacian, filename="10:10:10" ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        print(filename)
        libspm.spmReadDriver.argtypes = [c_int, c_char_p, POINTER(self.c_spm), c_int]
        libspm.spmReadDriver( driver, filename.encode('utf-8'), self.id_ptr, c_int(0) )

        # Assume A is already in Scipy sparse format
        print(self.spm_c.flttype)
        self.dtype = pastix_coeftype.getdtype( self.spm_c.flttype )

    def updateComputedField( self ):
        libspm.spmUpdateComputedFields.argtypes = [POINTER(self.c_spm)]
        libspm.spmUpdateComputedFields( self.id_ptr )

    def printInfo( self ):
        libspm.spmPrintInfo.argtypes = [POINTER(self.c_spm), c_void_p]
        libspm.spmPrintInfo( self.id_ptr, None )

    #def print( self ):
    #    libspm.spmPrint.argtypes = [POINTER(self.c_spm), c_void_p]
    #    libspm.spmPrint( self.id_ptr, None )

    def checkAndCorrect( self ):
        libspm.spmCheckAndCorrect.argtypes = [POINTER(self.c_spm)]
        libspm.spmCheckAndCorrect.restype = POINTER(self.c_spm)
        self.id_ptr = libspm.spmCheckAndCorrect( self.id_ptr )

    def __checkVector( self, n, nrhs, x ):
        if x.dtype != self.dtype:
            raise TypeError( "Vectors must use the same arithmetic as the spm" )

        if x.ndim > 2:
            raise TypeError( "Vectors must be of dimension 1 or 2" )

        if x.shape[0] < n:
            raise TypeError( "Vectors must be of dimension at least ", n )

        if x.shape[1] < nrhs:
            raise TypeError( "At least nrhs vectors must be stored in the vector" )

    def checkAxb( self, x0, b, x, nrhs=-1 ):
        if libspm == None:
            raise EnvironmentError( "SPM Instance badly instanciated" )

        n = self.spm_c.n
        if nrhs == -1:
            if x.ndim == 1:
                nrhs = 1
            else:
                nrhs = x.shape[1]

        self.__checkVector( n, nrhs, x0 )
        self.__checkVector( n, nrhs, b )
        self.__checkVector( n, nrhs, x )

        ldx0 = x0.shape[0]
        ldb  = b.shape[0]
        ldx  = x.shape[0]

        libspm.spmCheckAxb.argtypes = [ c_int, POINTER(self.c_spm), c_void_p, c_int, c_void_p, c_int, c_void_p, c_int]
        libspm.spmCheckAxb( nrhs, self.id_ptr,
                                 x0.ctypes.data_as(c_void_p), ldx0,
                                 b.ctypes.data_as(c_void_p),  ldb,
                                 x.ctypes.data_as(c_void_p),  ldx )
