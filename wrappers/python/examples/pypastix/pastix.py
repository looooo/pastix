"""
 @file pastix.py

 PaStiX python wrapper

 @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @date 2017-05-04

"""
from ctypes import *
import numpy as np
from . import libpastix
from .enum import *
from .spm import spm

def initParam():
    iparm_ = np.zeros( iparm.size, dtype=pastix_int)
    dparm_ = np.zeros( dparm.size, dtype='float64' )
    libpastix.pastixInitParam.argtypes = [POINTER(pastix_int), POINTER(c_double)]
    libpastix.pastixInitParam( iparm_.ctypes.data_as(POINTER(pastix_int)),
                               dparm_.ctypes.data_as(POINTER(c_double)) )
    return iparm_, dparm_

def init( iparm, dparm ):
    pastix_data = c_void_p();
    libpastix.pastixInit.argtypes = [POINTER(c_void_p), c_int, POINTER(pastix_int), POINTER(c_double)]
    libpastix.pastixInit( pointer( pastix_data ), c_int(0),
                          iparm.ctypes.data_as(POINTER(pastix_int)),
                          dparm.ctypes.data_as(POINTER(c_double)) )

    return pastix_data

def finalize( pastix_data, iparm, dparm ):
    libpastix.pastixFinalize.argtypes = [POINTER(c_void_p)]
    libpastix.pastixFinalize( pointer( pastix_data ) )

def __getnrhs(nrhs, x):
    if nrhs == -1:
        if x.ndim == 1:
            nrhs = 1
        else:
            nrhs = x.shape[1]
    return nrhs

#
# Pastix Tasks
#
def task_analyze( pastix_data, spm ):
    libpastix.pastix_task_analyze.argtypes = [c_void_p, POINTER(spm.c_spm)]
    libpastix.pastix_task_analyze( pastix_data, spm.id_ptr )

def task_numfact( pastix_data, spm ):
    libpastix.pastix_task_numfact.argtypes = [c_void_p, POINTER(spm.c_spm)]
    libpastix.pastix_task_numfact( pastix_data, spm.id_ptr )

def task_solve( pastix_data, x, nrhs=-1 ):

    n    = x.shape[0]
    nrhs = __getnrhs( nrhs, x )

    x = np.asarray(x, spm.dtype)

    libpastix.pastix_task_solve.argtypes = [c_void_p, pastix_int, c_void_p, pastix_int]
    libpastix.pastix_task_solve( pastix_data, nrhs, x.ctypes.data_as(c_void_p), n )

def task_refine( pastix_data, b, x, nrhs=-1 ):

    nrhs = __getnrhs( nrhs, x )
    b = np.asarray( b, spm.dtype )
    x = np.asarray( x, spm.dtype )

    if b.dtype != x.dtype:
        raise TypeError( "b and x must use the same arithmetic" )

    libpastix.pastix_task_refine.argtypes = [c_void_p, c_void_p, pastix_int, c_void_p]
    libpastix.pastix_task_refine( pastix_data,
                                  x.ctypes.data_as(c_void_p), nrhs,
                                  b.ctypes.data_as(c_void_p) )

#
# Numerical solve subtasks
#
def subtask_applyorder( pastix_data, direction, b ):
    flttype = coeftype.getptype( b.dtype )

    m = b.shape[0]
    if b.ndim == 1:
        n = 1
    else:
        n = b.shape[1]
    ldb = m

    libpastix.pastix_subtask_applyorder.argtypes = [c_void_p, c_int, c_int, pastix_int, pastix_int, c_void_p, pastix_int]
    libpastix.pastix_subtask_applyorder( pastix_data, flttype, direction, m, n,
                                         b.ctypes.data_as(c_void_p), ldb )

def subtask_trsm( pastix_data, side, uplo, trans, diag, b, nrhs=-1 ):
    flttype = coeftype.getptype( b.dtype )

    n    = b.shape[0]
    nrhs = __getnrhs( nrhs, b )
    ldb  = n

    libpastix.pastix_subtask_trsm.argtypes = [c_void_p, c_int, c_int, c_int, c_int, c_int, pastix_int, c_void_p, pastix_int]
    libpastix.pastix_subtask_trsm( pastix_data, flttype, side, uplo, trans, diag, nrhs,
                                   b.ctypes.data_as(c_void_p), ldb )

def subtask_diag( pastix_data, b ):
    flttype = coeftype.getptype( b.dtype )

    ldb  = b.shape[0]
    nrhs = __getnrhs( nrhs, b )

    libpastix.pastix_subtask_diag.argtypes = [c_void_p, c_int, pastix_int, c_void_p, pastix_int]
    libpastix.pastix_subtask_diag( pastix_data, flttype, nrhs,
                                   b.ctypes.data_as(c_void_p), ldb )

#
# Schur complement manipulation routines.
#
def setSchurUnknownList( pastix_data, schur_list ):
    n = schur_list.shape[0]
    schur_list = np.array(schur_list, dtype=pastix_int)

    libpastix.pastix_setSchurUnknownList.argtypes = [c_void_p, pastix_int, POINTER(pastix_int)]
    libpastix.pastix_setSchurUnknownList( pastix_data, n,
                                          schur_list.ctypes.data_as(POINTER(pastix_int)) )

def getSchur( pastix_data, S ):
    libpastix.pastix_getSchur.argtypes = [c_void_p, c_void_p, pastix_int ]
    libpastix.pastix_getSchur( pastix_data,
                               S.ctypes.data_as(c_void_p),
                               S.shape[0] )
