"""
 @file pypastix.py

 PaStiX python wrapper

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
from pyspm import *
import os
import numpy as np
import scipy.sparse as spp

class pastix():

    libpastix = None
    _c_int = c_int64
    _dtype = np.dtype('int64')

    @staticmethod
    def __init__( ):
        """Initialize the PaStiX wrapper by loading the libraries"""
        if pastix.libpastix == None:
            pastix.libpastix = pastix.__load_pastix();

    @staticmethod
    def __load_pastix( ):
        """
        Load the PASTIX library
        """

        pastix.libpastix = find_library('pastix')
        if pastix.libpastix == None:
            raise EnvironmentError("Could not find shared library: pastix")
        return cdll.LoadLibrary(pastix.libpastix)

    @staticmethod
    def initParam( iparm, dparm ):
        if pastix.libpastix == None:
            pastix.__init__()

        pastix.libpastix.pastixInitParam.argtypes = [POINTER(pastix_c_int), POINTER(c_double)]
        pastix.libpastix.pastixInitParam( iparm.ctypes.data_as(POINTER(pastix_c_int)),
                                          dparm.ctypes.data_as(POINTER(c_double)) )

    @staticmethod
    def init( iparm, dparm ):
        if pastix.libpastix == None:
            pastix.__init__()

        pastix_data = c_void_p();
        pastix.libpastix.pastixInit.argtypes = [POINTER(c_void_p), c_int, POINTER(pastix_c_int), POINTER(c_double)]
        pastix.libpastix.pastixInit( pointer( pastix_data ), c_int(0),
                                     iparm.ctypes.data_as(POINTER(pastix_c_int)),
                                     dparm.ctypes.data_as(POINTER(c_double)) )

        return pastix_data

    @staticmethod
    def finalize( pastix_data, iparm, dparm ):

        pastix.libpastix.pastixFinalize.argtypes = [POINTER(c_void_p), c_int, POINTER(pastix_c_int), POINTER(c_double)]
        pastix.libpastix.pastixFinalize( pointer( pastix_data ), c_int(0),
                                         iparm.ctypes.data_as(POINTER(pastix_c_int)),
                                         dparm.ctypes.data_as(POINTER(c_double)) )

    #
    # Pastix Tasks
    #
    @staticmethod
    def analyze( pastix_data, spm ):
        if pastix.libpastix == None:
            pastix.__init__()

        pastix.libpastix.pastix_task_analyze.argtypes = [c_void_p, POINTER(spm.c_spm)]
        pastix.libpastix.pastix_task_analyze( pastix_data, spm.id_ptr )

    @staticmethod
    def numfact( pastix_data, spm ):
        if pastix.libpastix == None:
            pastix.__init__()

        pastix.libpastix.pastix_task_numfact.argtypes = [c_void_p, POINTER(spm.c_spm)]
        pastix.libpastix.pastix_task_numfact( pastix_data, spm.id_ptr )

    @staticmethod
    def solve( pastix_data, spm, nrhs, x, ldx ):
        if pastix.libpastix == None:
            pastix.__init__()

        if x.dtype != spm.dtype:
            raise TypeError( "x must use the same arithmetic as the spm" )

        pastix.libpastix.pastix_task_solve.argtypes = [c_void_p, POINTER(spm.c_spm), pastix_c_int, c_void_p, pastix_c_int]
        pastix.libpastix.pastix_task_solve( pastix_data, spm.id_ptr, nrhs,
                                            x.ctypes.data_as(c_void_p), ldx )


    @staticmethod
    def refine( pastix_data, nrhs, b, ldb, x, ldx ):
        if pastix.libpastix == None:
            pastix.__init__()

        if b.dtype != spm.dtype:
            raise TypeError( "b must use the same arithmetic as the spm" )

        if x.dtype != spm.dtype:
            raise TypeError( "x must use the same arithmetic as the spm" )

        pastix.libpastix.pastix_task_refine.argtypes = [c_void_p, c_void_p, pastix_c_int, c_void_p]
        pastix.libpastix.pastix_task_refine( pastix_data,
                                             x.ctypes.data_as(c_void_p), nrhs,
                                             b.ctypes.data_as(c_void_p) )



    #
    # Pastix Schur
    #
    @staticmethod
    def setSchurUnknownList( pastix_data, n, list ):
        if pastix.libpastix == None:
            pastix.__init__()

        if _dtype != list.dtype:
            raise TypeError( "list must use the same integer type as the pastix library" )

        pastix.libpastix.pastix_setSchurUnknownList.argtypes = [c_void_p, pastix_c_int, POINTER(pastix_c_int)]
        pastix.libpastix.pastix_setSchurUnknownList( pastix_data, n,
                                                     list.ctypes.data_as(POINTER(pastix_c_int)) )

    @staticmethod
    def getSchur( pastix_data, S ):
        if pastix.libpastix == None:
            pastix.__init__()

        if _dtype != list.dtype:
            raise TypeError( "list must use the same integer type as the pastix library" )

        pastix.libpastix.pastix_getSchur.argtypes = [c_void_p, c_void_p, pastix_c_int ]
        pastix.libpastix.pastix_getSchur.restype  = c_void_p
        S = np.array( pastix.libpastix.pastix_getSchur( pastix_data,
                                                        S.ctypes.data_as(c_void_p),
                                                        S.shape[0] ), order='F' ).reshape( S.shape[0], S.shape[1] )

    #
    # Pastix Schur
    #
    @staticmethod
    def subtask_applyorder( pastix_data, dir, m, n, b, ldb ):
        if pastix.libpastix == None:
            pastix.__init__()

        flttype = pastix_coeftype.get( b.dtype )

        pastix.libpastix.pastix_pastix_subtask_applyorder.argtypes = [c_void_p, c_int, c_int, pastix_c_int, pastix_c_int, c_void_p, pastix_c_int]
        pastix.libpastix.pastix_pastix_subtask_applyorder( pastix_data, flttype, dir, m, n,
                                                           b.ctypes.data_as(c_void_p), ldb )

    @staticmethod
    def subtask_trsm( pastix_data, side, uplo, trans, diag, nrhs, b, ldb ):
        if pastix.libpastix == None:
            pastix.__init__()

        flttype = pastix_coeftype.get( b.dtype )

        pastix.libpastix.pastix_pastix_subtask_trsm.argtypes = [c_void_p, c_int, c_int, c_int, c_int, c_int, pastix_c_int, c_void_p, pastix_c_int]
        pastix.libpastix.pastix_pastix_subtask_trsm( pastix_data, flttype, side, uplo, trans, diag, nrhs,
                                                     b.ctypes.data_as(c_void_p), ldb )

    @staticmethod
    def subtask_diag( pastix_data, nrhs, b, ldb ):
        if pastix.libpastix == None:
            pastix.__init__()

        flttype = pastix_coeftype.get( b.dtype )

        pastix.libpastix.pastix_pastix_subtask_diag.argtypes = [c_void_p, c_int, pastix_c_int, c_void_p, pastix_c_int]
        pastix.libpastix.pastix_pastix_subtask_diag( pastix_data, flttype, nrhs,
                                                     b.ctypes.data_as(c_void_p), ldb )

