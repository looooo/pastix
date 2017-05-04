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

    def __init__( self, A, mtxtype=pastix_mtxtype.PastixGeneral ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        if self.libspm == None:
            self.libspm = self.__load_spm();

        # Assume A is already in Scipy sparse format
        self.dtype = A.dtype
        flttype = pastix_coeftype.get( A.dtype )
        print( "Floating point arithmetic is", flttype )
        if flttype == -1:
            raise TypeError( "Invalid data type. Must be part of (f4, f8, c8 or c16)" )

        A  = spp.csc_matrix( A )
        A.sort_indices()
        n   = A.shape[0]
        nnz = A.getnnz()

        # Pointer variables
        self.py_dofs      = None
        self.py_colptr    = np.array( A.indptr[:],  dtype=pastix_np_int )
        self.py_rowptr    = np.array( A.indices[:], dtype=pastix_np_int )
        self.py_values    = np.array( A.data[:] )
        self.py_loc2glob  = None

        self.spm_c = self.c_spm( mtxtype, flttype, pastix_fmttype.PastixCSC,
                                 0, n, 0, nnz, 0, 0, 0, 0,
                                 1, None,
                                 pastix_order.PastixColMajor,
                                 self.py_colptr.ctypes.data_as(POINTER(pastix_c_int)),
                                 self.py_rowptr.ctypes.data_as(POINTER(pastix_c_int)),
                                 None,
                                 self.py_values.ctypes.data_as(c_void_p))

        self.id_ptr = pointer( self.spm_c )

        self.libspm.spmUpdateComputedFields.argtypes = [POINTER(self.c_spm)]
        self.libspm.spmUpdateComputedFields( self.id_ptr )

    def __load_spm(self):
        """
        Load the SPM library
        """

        self.libspm = find_library('pastix_spm')
        if self.libspm == None:
            raise EnvironmentError("Could not find shared library: pastix_spm")
        __spm_loaded = 1
        return cdll.LoadLibrary(self.libspm)


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
