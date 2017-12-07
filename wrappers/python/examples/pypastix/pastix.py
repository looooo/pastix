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
from .enum import *
from .spm import spm
from .__pastix__ import *

def initParam():
    iparm_ = np.zeros( iparm.size, dtype=pastix_int)
    dparm_ = np.zeros( dparm.size, dtype='float64' )
    pypastix_pastixInitParam( iparm_, dparm_ )
    return iparm_, dparm_

def init( iparm, dparm ):
    pastix_data = c_void_p();
    pypastix_pastixInit( pastix_data, c_int(0), iparm, dparm )
    return pastix_data

def finalize( pastix_data, iparm, dparm ):
    pypastix_pastixFinalize( pastix_data )

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
    pypastix_pastix_task_analyze( pastix_data, spm.id_ptr )

def task_numfact( pastix_data, spm ):
    pypastix_pastix_task_numfact( pastix_data, spm.id_ptr )

def task_solve( pastix_data, x, nrhs=-1 ):

    n    = x.shape[0]
    nrhs = __getnrhs( nrhs, x )

    x = np.asarray(x, spm.dtype)

    pypastix_pastix_task_solve( pastix_data, nrhs, x.ctypes.data_as(c_void_p), n )

def task_refine( pastix_data, b, x, nrhs=-1 ):

    nrhs = __getnrhs( nrhs, x )
    b = np.asarray( b, spm.dtype )
    x = np.asarray( x, spm.dtype )

    if b.dtype != x.dtype:
        raise TypeError( "b and x must use the same arithmetic" )

    pypastix_pastix_task_refine( pastix_data,
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

    pypastix_pastix_subtask_applyorder( pastix_data, flttype, direction, m, n,
                                         b.ctypes.data_as(c_void_p), ldb )

def subtask_trsm( pastix_data, side, uplo, trans, diag, b, nrhs=-1 ):
    flttype = coeftype.getptype( b.dtype )

    n    = b.shape[0]
    nrhs = __getnrhs( nrhs, b )
    ldb  = n

    pypastix_pastix_subtask_trsm( pastix_data, flttype, side, uplo, trans, diag, nrhs,
                                  b.ctypes.data_as(c_void_p), ldb )

def subtask_diag( pastix_data, b ):
    flttype = coeftype.getptype( b.dtype )

    ldb  = b.shape[0]
    nrhs = __getnrhs( nrhs, b )

    pypastix_pastix_subtask_diag( pastix_data, flttype, nrhs,
                                  b.ctypes.data_as(c_void_p), ldb )

#
# Schur complement manipulation routines.
#
def setSchurUnknownList( pastix_data, schur_list ):
    n = schur_list.shape[0]
    schur_list = np.array(schur_list, dtype=pastix_int)
    pypastix_pastixSetSchurUnknownList( pastix_data, n,
                                        schur_list.ctypes.data_as(POINTER(pastix_int)) )

def getSchur( pastix_data, S ):
    pypastix_pastixGetSchur( pastix_data,
                             S.ctypes.data_as(c_void_p),
                             S.shape[0] )
