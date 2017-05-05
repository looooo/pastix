"""
 @file pyspm.py

 SPM python wrapper

 @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @date 2017-05-04

"""
from ctypes import *
from ctypes.util import find_library
from pastix_api import *
import numpy as np
import scipy.sparse as spp

class spm():

    libspm = None
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

    def __load_spm(self):
        """
        Load the SPM library
        """

        self.libspm = find_library('pastix_spm')
        if self.libspm == None:
            raise EnvironmentError("Could not find shared library: pastix_spm")
        __spm_loaded = 1
        return cdll.LoadLibrary(self.libspm)

    def __init__( self ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        if self.libspm == None:
            self.libspm = self.__load_spm();

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

        self.libspm.spmUpdateComputedFields.argtypes = [POINTER(self.c_spm)]
        self.libspm.spmUpdateComputedFields( self.id_ptr )

    def fromdriver( self, driver=pastix_driver.PastixDriverLaplacian, filename="10:10:10" ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        print(filename)
        self.libspm.spmReadDriver.argtypes = [c_int, c_char_p, POINTER(self.c_spm), c_int]
        self.libspm.spmReadDriver( driver, filename.encode('utf-8'), self.id_ptr, c_int(0) )

        # Assume A is already in Scipy sparse format
        print(self.spm_c.flttype)
        self.dtype = pastix_coeftype.getdtype( self.spm_c.flttype )

    def printInfo( self ):
        self.libspm.spmPrintInfo.argtypes = [POINTER(self.c_spm), c_void_p]
        self.libspm.spmPrintInfo( self.id_ptr, None )

    def checkAxb( self, nrhs, x0, ldx0, b, ldb, x, ldx ):
        if self.libspm == None:
            raise EnvironmentError( "SPM Instance badly instanciated" )

        if x0.dtype != self.dtype:
            raise TypeError( "b must use the same arithmetic as the spm" )

        if b.dtype != self.dtype:
            raise TypeError( "b must use the same arithmetic as the spm" )

        if x.dtype != self.dtype:
            raise TypeError( "x must use the same arithmetic as the spm" )

        self.libspm.spmCheckAxb.argtypes = [ c_int, POINTER(self.c_spm), c_void_p, c_int, c_void_p, c_int, c_void_p, c_int]
        self.libspm.spmCheckAxb( nrhs, self.id_ptr,
                                 x0.ctypes.data_as(c_void_p), ldx0,
                                 b.ctypes.data_as(c_void_p),  ldb,
                                 x.ctypes.data_as(c_void_p),  ldx )
