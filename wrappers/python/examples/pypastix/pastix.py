"""
 @file pastix.py

 PaStiX python wrapper

 @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.4.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2024-07-05

"""
from ctypes import *
import numpy as np
import spm
from .enum       import *
from .enum       import __pastix_int__
from .__pastix__ import *

def initParam():
    iparm_ = np.zeros( iparm.size, dtype=__pastix_int__)
    dparm_ = np.zeros( dparm.size, dtype='float64' )
    pypastix_pastixInitParam( iparm_, dparm_ )
    return iparm_, dparm_

def init( iparm, dparm, comm=pypastix_default_comm ):
    pastix_data = c_void_p()
    pypastix_pastixInit( pastix_data, comm, iparm, dparm )
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
# Order structure
#
def orderGet( pastix_data ):
    return pypastix_pastixOrderGet( pastix_data )

def orderGetPermtab( pastix_data, order = None ):
    if not order:
        order = orderGet( pastix_data )

    pyorder = cast( order, POINTER(pypastix_order_t) ).contents

    perm = np.frombuffer( (__pastix_int__ * (pyorder.vertnbr)).from_address( cast(pyorder.permtab, c_void_p).value ), dtype=__pastix_int__ )
    return perm

def orderGetPeritab( pastix_data, order = None ):
    if not order:
        order = orderGet( pastix_data )

    pyorder = cast( order, POINTER(pypastix_order_t) ).contents

    invp = np.frombuffer( (__pastix_int__ * (pyorder.vertnbr)).from_address( cast(pyorder.peritab, c_void_p).value ), dtype=__pastix_int__ )
    return invp

#
# Pastix Tasks
#
def task_analyze( pastix_data, spm ):
    pypastix_pastix_task_analyze( pastix_data, spm.id_ptr )

def task_numfact( pastix_data, spm ):
    pypastix_pastix_task_numfact( pastix_data, spm.id_ptr )

def task_solve( pastix_data, spm, x, nrhs=-1 ):

    m    = x.shape[0]
    nrhs = __getnrhs( nrhs, x )

    x = np.asarray(x, spm.dtype)

    pypastix_pastix_task_solve( pastix_data, m, nrhs, x.ctypes.data_as(c_void_p), m )

def task_refine( pastix_data, spm, b, x, nrhs=-1 ):

    n    = x.shape[0]
    nrhs = __getnrhs( nrhs, x )
    b    = np.asarray( b, spm.dtype )
    ldb  = b.shape[0]
    x    = np.asarray( x, spm.dtype )
    ldx  = x.shape[0]

    if b.dtype != x.dtype:
        raise TypeError( "b and x must use the same arithmetic" )

    pypastix_pastix_task_refine( pastix_data, n, nrhs,
                                 b.ctypes.data_as(c_void_p), ldb,
                                 x.ctypes.data_as(c_void_p), ldx )

#
# Numerical solve subtasks
#
def rhsInit():
    rhs = c_void_p();
    pypastix_pastixRhsInit( rhs )
    return rhs

def rhsFinalize( rhs ):
    pypastix_pastixRhsFinalize( rhs )

def subtask_applyorder( pastix_data, direction, b, Bp ):
    m = b.shape[0]
    if b.ndim == 1:
        n = 1
    else:
        n = b.shape[1]
    ldb = m

    pypastix_pastix_subtask_applyorder( pastix_data, direction, m, n,
                                        b.ctypes.data_as(c_void_p), ldb, Bp )

def subtask_trsm( pastix_data, side, uplo, trans, diag, b ):
    pypastix_pastix_subtask_trsm( pastix_data, side, uplo, trans, diag, b )

def subtask_diag( pastix_data, b ):
    pypastix_pastix_subtask_diag( pastix_data, b )

#
# Schur complement manipulation routines.
#
def setSchurUnknownList( pastix_data, schur_list ):
    n = schur_list.shape[0]
    schur_list = np.array(schur_list, dtype=__pastix_int__)
    pypastix_pastixSetSchurUnknownList( pastix_data, n, schur_list )

def getSchur( pastix_data, S ):
    pypastix_pastixGetSchur( pastix_data,
                             S.ctypes.data_as(c_void_p),
                             S.shape[0] )

def rhsSchurGet( pastix_data, rhsB, B ):
    m = B.shape[0]
    if B.ndim == 1:
        n = 1
    else:
        n = B.shape[1]
    ldb = m
    pypastix_pastixRhsSchurGet( pastix_data, m, n, rhsB,
                                B.ctypes.data_as(c_void_p), ldb )

def rhsSchurSet( pastix_data, B, rhsB ):
    m = B.shape[0]
    if B.ndim == 1:
        n = 1
    else:
        n = B.shape[1]
    ldb = m
    pypastix_pastixRhsSchurSet( pastix_data, m, n,
                                B.ctypes.data_as(c_void_p), ldb, rhsB )

