"""

 @file __pastix__.py

 PaStiX python wrapper

 @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-05-04

This file has been automatically generated with gen_wrappers.py

"""
from ctypes import *
import numpy as np


from . import libpastix
from .enum import pastix_int
from .__spm__ import pypastix_spm_t

class pypastix_order_t(Structure):
    _fields_ = [("baseval",  pastix_int         ),
                ("vertnbr",  pastix_int         ),
                ("cblknbr",  pastix_int         ),
                ("permtab",  POINTER(pastix_int)),
                ("peritab",  POINTER(pastix_int)),
                ("rangtab",  POINTER(pastix_int)),
                ("treetab",  POINTER(pastix_int)) ]

def pypastix_pastixOrderInit( ordeptr, baseval, vertnbr, cblknbr, perm, invp,
                              rang, tree ):
    libpastix.pastixOrderInit.argtypes = [ c_void_p, pastix_int, pastix_int,
                                           pastix_int, POINTER(pastix_int),
                                           POINTER(pastix_int),
                                           POINTER(pastix_int),
                                           POINTER(pastix_int) ]
    libpastix.pastixOrderInit.restype = c_int
    return libpastix.pastixOrderInit( ordeptr, baseval, vertnbr, cblknbr,
                                      perm.ctypes.data_as( POINTER(pastix_int) ),
                                      invp.ctypes.data_as( POINTER(pastix_int) ),
                                      rang.ctypes.data_as( POINTER(pastix_int) ),
                                      tree.ctypes.data_as( POINTER(pastix_int) ) )

def pypastix_pastixOrderAlloc( ordeptr, vertnbr, cblknbr ):
    libpastix.pastixOrderAlloc.argtypes = [ c_void_p, pastix_int, pastix_int ]
    libpastix.pastixOrderAlloc.restype = c_int
    return libpastix.pastixOrderAlloc( ordeptr, vertnbr, cblknbr )

def pypastix_pastixOrderExit( ordeptr ):
    libpastix.pastixOrderExit.argtypes = [ c_void_p ]
    libpastix.pastixOrderExit( ordeptr )

def pypastix_pastixOrderBase( ordeptr, baseval ):
    libpastix.pastixOrderBase.argtypes = [ c_void_p, pastix_int ]
    libpastix.pastixOrderBase( ordeptr, baseval )

def pypastix_pastixOrderCheck( ordeptr ):
    libpastix.pastixOrderCheck.argtypes = [ c_void_p ]
    libpastix.pastixOrderCheck.restype = c_int
    return libpastix.pastixOrderCheck( ordeptr )

def pypastix_pastixOrderCopy( ordedst, ordesrc ):
    libpastix.pastixOrderCopy.argtypes = [ c_void_p, c_void_p ]
    libpastix.pastixOrderCopy.restype = c_int
    return libpastix.pastixOrderCopy( ordedst, ordesrc )

def pypastix_pastixOrderGet( pastix_data ):
    libpastix.pastixOrderGet.argtypes = [ c_void_p ]
    libpastix.pastixOrderGet.restype = c_void_p
    return libpastix.pastixOrderGet( pastix_data )

def pypastix_pastixOrderLoad( pastix_data, ordeptr ):
    libpastix.pastixOrderLoad.argtypes = [ c_void_p, c_void_p ]
    libpastix.pastixOrderLoad.restype = c_int
    return libpastix.pastixOrderLoad( pastix_data, ordeptr )

def pypastix_pastixOrderSave( pastix_data, ordeptr ):
    libpastix.pastixOrderSave.argtypes = [ c_void_p, c_void_p ]
    libpastix.pastixOrderSave.restype = c_int
    return libpastix.pastixOrderSave( pastix_data, ordeptr )

def pypastix_pastix( pastix_data, pastix_comm, n, colptr, row, avals, perm,
                     invp, b, nrhs, iparm, dparm ):
    libpastix.pastix.argtypes = [ c_void_p, c_int, pastix_int,
                                  POINTER(pastix_int), POINTER(pastix_int),
                                  c_void_p, POINTER(pastix_int),
                                  POINTER(pastix_int), c_void_p, pastix_int,
                                  POINTER(pastix_int), POINTER(c_double) ]
    libpastix.pastix.restype = c_int
    return libpastix.pastix( pointer( pastix_data ), pastix_comm, n,
                             colptr.ctypes.data_as( POINTER(pastix_int) ),
                             row.ctypes.data_as( POINTER(pastix_int) ), avals,
                             perm.ctypes.data_as( POINTER(pastix_int) ),
                             invp.ctypes.data_as( POINTER(pastix_int) ), b,
                             nrhs, iparm.ctypes.data_as( POINTER(pastix_int) ),
                             dparm.ctypes.data_as( POINTER(c_double) ) )

def pypastix_pastixInitParam( iparm, dparm ):
    libpastix.pastixInitParam.argtypes = [ POINTER(pastix_int),
                                           POINTER(c_double) ]
    libpastix.pastixInitParam( iparm.ctypes.data_as( POINTER(pastix_int) ),
                               dparm.ctypes.data_as( POINTER(c_double) ) )

def pypastix_pastixInit( pastix_data, pastix_comm, iparm, dparm ):
    libpastix.pastixInit.argtypes = [ c_void_p, c_int, POINTER(pastix_int),
                                      POINTER(c_double) ]
    libpastix.pastixInit( pointer( pastix_data ), pastix_comm,
                          iparm.ctypes.data_as( POINTER(pastix_int) ),
                          dparm.ctypes.data_as( POINTER(c_double) ) )

def pypastix_pastixInitWithAffinity( pastix_data, pastix_comm, iparm, dparm,
                                     bindtab ):
    libpastix.pastixInitWithAffinity.argtypes = [ c_void_p, c_int,
                                                  POINTER(pastix_int),
                                                  POINTER(c_double), c_int_p ]
    libpastix.pastixInitWithAffinity( pointer( pastix_data ), pastix_comm,
                                      iparm.ctypes.data_as( POINTER(pastix_int) ),
                                      dparm.ctypes.data_as( POINTER(c_double) ),
                                      bindtab.ctypes.data_as( c_int_p ) )

def pypastix_pastixFinalize( pastix_data ):
    libpastix.pastixFinalize.argtypes = [ c_void_p ]
    libpastix.pastixFinalize( pointer( pastix_data ) )

def pypastix_pastix_task_analyze( pastix_data, spm ):
    libpastix.pastix_task_analyze.argtypes = [ c_void_p,
                                               POINTER(pypastix_spm_t) ]
    libpastix.pastix_task_analyze.restype = c_int
    return libpastix.pastix_task_analyze( pastix_data, spm )

def pypastix_pastix_task_numfact( pastix_data, spm ):
    libpastix.pastix_task_numfact.argtypes = [ c_void_p,
                                               POINTER(pypastix_spm_t) ]
    libpastix.pastix_task_numfact.restype = c_int
    return libpastix.pastix_task_numfact( pastix_data, spm )

def pypastix_pastix_task_solve( pastix_data, nrhs, b, ldb ):
    libpastix.pastix_task_solve.argtypes = [ c_void_p, pastix_int, c_void_p,
                                             pastix_int ]
    libpastix.pastix_task_solve.restype = c_int
    return libpastix.pastix_task_solve( pastix_data, nrhs, b, ldb )

def pypastix_pastix_task_refine( pastix_data, x, nrhs, b ):
    libpastix.pastix_task_refine.argtypes = [ c_void_p, c_void_p, pastix_int,
                                              c_void_p ]
    libpastix.pastix_task_refine.restype = c_int
    return libpastix.pastix_task_refine( pastix_data, x, nrhs, b )

def pypastix_pastix_subtask_order( pastix_data, spm, myorder ):
    libpastix.pastix_subtask_order.argtypes = [ c_void_p,
                                                POINTER(pypastix_spm_t),
                                                c_void_p ]
    libpastix.pastix_subtask_order.restype = c_int
    return libpastix.pastix_subtask_order( pastix_data, spm, myorder )

def pypastix_pastix_subtask_symbfact( pastix_data ):
    libpastix.pastix_subtask_symbfact.argtypes = [ c_void_p ]
    libpastix.pastix_subtask_symbfact.restype = c_int
    return libpastix.pastix_subtask_symbfact( pastix_data )

def pypastix_pastix_subtask_reordering( pastix_data ):
    libpastix.pastix_subtask_reordering.argtypes = [ c_void_p ]
    libpastix.pastix_subtask_reordering.restype = c_int
    return libpastix.pastix_subtask_reordering( pastix_data )

def pypastix_pastix_subtask_blend( pastix_data ):
    libpastix.pastix_subtask_blend.argtypes = [ c_void_p ]
    libpastix.pastix_subtask_blend.restype = c_int
    return libpastix.pastix_subtask_blend( pastix_data )

def pypastix_pastix_subtask_spm2bcsc( pastix_data, spm ):
    libpastix.pastix_subtask_spm2bcsc.argtypes = [ c_void_p,
                                                   POINTER(pypastix_spm_t) ]
    libpastix.pastix_subtask_spm2bcsc.restype = c_int
    return libpastix.pastix_subtask_spm2bcsc( pastix_data, spm )

def pypastix_pastix_subtask_bcsc2ctab( pastix_data ):
    libpastix.pastix_subtask_bcsc2ctab.argtypes = [ c_void_p ]
    libpastix.pastix_subtask_bcsc2ctab.restype = c_int
    return libpastix.pastix_subtask_bcsc2ctab( pastix_data )

def pypastix_pastix_subtask_sopalin( pastix_data ):
    libpastix.pastix_subtask_sopalin.argtypes = [ c_void_p ]
    libpastix.pastix_subtask_sopalin.restype = c_int
    return libpastix.pastix_subtask_sopalin( pastix_data )

def pypastix_pastix_subtask_applyorder( pastix_data, flttype, dir, m, n, b, ldb ):
    libpastix.pastix_subtask_applyorder.argtypes = [ c_void_p, c_int, c_int,
                                                     pastix_int, pastix_int,
                                                     c_void_p, pastix_int ]
    libpastix.pastix_subtask_applyorder.restype = c_int
    return libpastix.pastix_subtask_applyorder( pastix_data, flttype, dir, m, n,
                                                b, ldb )

def pypastix_pastix_subtask_trsm( pastix_data, flttype, side, uplo, trans, diag,
                                  nrhs, b, ldb ):
    libpastix.pastix_subtask_trsm.argtypes = [ c_void_p, c_int, c_int, c_int,
                                               c_int, c_int, pastix_int,
                                               c_void_p, pastix_int ]
    libpastix.pastix_subtask_trsm.restype = c_int
    return libpastix.pastix_subtask_trsm( pastix_data, flttype, side, uplo,
                                          trans, diag, nrhs, b, ldb )

def pypastix_pastix_subtask_diag( pastix_data, flttype, nrhs, b, ldb ):
    libpastix.pastix_subtask_diag.argtypes = [ c_void_p, c_int, pastix_int,
                                               c_void_p, pastix_int ]
    libpastix.pastix_subtask_diag.restype = c_int
    return libpastix.pastix_subtask_diag( pastix_data, flttype, nrhs, b, ldb )

def pypastix_pastixSetSchurUnknownList( pastix_data, n, list ):
    libpastix.pastixSetSchurUnknownList.argtypes = [ c_void_p, pastix_int,
                                                     POINTER(pastix_int) ]
    libpastix.pastixSetSchurUnknownList( pastix_data, n, list )

def pypastix_pastixGetSchur( pastix_data, S, lds ):
    libpastix.pastixGetSchur.argtypes = [ c_void_p, c_void_p, pastix_int ]
    libpastix.pastixGetSchur.restype = c_int
    return libpastix.pastixGetSchur( pastix_data, S, lds )

def pypastix_pastixGetDiag( pastix_data, D, incD ):
    libpastix.pastixGetDiag.argtypes = [ c_void_p, c_void_p, pastix_int ]
    libpastix.pastixGetDiag.restype = c_int
    return libpastix.pastixGetDiag( pastix_data, D, incD )

def pypastix_pastixGetOptions( argc, argv, iparm, dparm, check, driver,
                               filename ):
    libpastix.pastixGetOptions.argtypes = [ c_int, c_char_p,
                                            POINTER(pastix_int),
                                            POINTER(c_double), c_int_p, c_int_p,
                                            c_char_p ]
    libpastix.pastixGetOptions( argc, pointer( argv ),
                                iparm.ctypes.data_as( POINTER(pastix_int) ),
                                dparm.ctypes.data_as( POINTER(c_double) ),
                                check, driver, pointer( filename ) )


