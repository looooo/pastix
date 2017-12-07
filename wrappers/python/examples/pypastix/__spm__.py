"""

 @file __spm__.py

 SPM python wrapper

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

from . import libspm
from .enum import pastix_int

class pypastix_spm_t(Structure):
    _fields_ = [("mtxtype",   c_int              ),
                ("flttype",   c_int              ),
                ("fmttype",   c_int              ),
                ("gN",        pastix_int         ),
                ("n",         pastix_int         ),
                ("gnnz",      pastix_int         ),
                ("nnz",       pastix_int         ),
                ("gNexp",     pastix_int         ),
                ("nexp",      pastix_int         ),
                ("gnnzexp",   pastix_int         ),
                ("nnzexp",    pastix_int         ),
                ("dof",       pastix_int         ),
                ("dofs",      POINTER(pastix_int)),
                ("layout",    c_int              ),
                ("colptr",    POINTER(pastix_int)),
                ("rowptr",    POINTER(pastix_int)),
                ("loc2glob",  POINTER(pastix_int)),
                ("values",    c_void_p           ) ]

def pyspm_spmInit( spm ):
    libspm.spmInit.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmInit( spm )

def pyspm_spmExit( spm ):
    libspm.spmExit.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmExit( spm )

def pyspm_spmCopy( spm ):
    libspm.spmCopy.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmCopy.restype = POINTER(pypastix_spm_t)
    return libspm.spmCopy( spm )

def pyspm_spmBase( spm, baseval ):
    libspm.spmBase.argtypes = [ POINTER(pypastix_spm_t), c_int ]
    libspm.spmBase( spm, baseval )

def pyspm_spmFindBase( spm ):
    libspm.spmFindBase.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmFindBase.restype = pastix_int
    return libspm.spmFindBase( spm )

def pyspm_spmConvert( ofmttype, ospm ):
    libspm.spmConvert.argtypes = [ c_int, POINTER(pypastix_spm_t) ]
    libspm.spmConvert.restype = c_int
    return libspm.spmConvert( ofmttype, ospm )

def pyspm_spmUpdateComputedFields( spm ):
    libspm.spmUpdateComputedFields.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmUpdateComputedFields( spm )

def pyspm_spmGenFakeValues( spm ):
    libspm.spmGenFakeValues.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmGenFakeValues( spm )

def pyspm_spmNorm( ntype, spm ):
    libspm.spmNorm.argtypes = [ c_int, POINTER(pypastix_spm_t) ]
    libspm.spmNorm.restype = c_double
    return libspm.spmNorm( ntype, spm )

def pyspm_spmMatVec( trans, alpha, spm, x, beta, y ):
    libspm.spmMatVec.argtypes = [ c_int, c_void_p, POINTER(pypastix_spm_t),
                                  c_void_p, c_void_p, c_void_p ]
    libspm.spmMatVec.restype = c_int
    return libspm.spmMatVec( trans, alpha, spm, x, beta, y )

def pyspm_spmScalMatrix( alpha, spm ):
    libspm.spmScalMatrix.argtypes = [ c_double, POINTER(pypastix_spm_t) ]
    libspm.spmScalMatrix( alpha, spm )

def pyspm_spmScalVector( alpha, spm, x ):
    libspm.spmScalVector.argtypes = [ c_double, POINTER(pypastix_spm_t),
                                      c_void_p ]
    libspm.spmScalVector( alpha, spm, x )

def pyspm_spmSort( spm ):
    libspm.spmSort.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmSort.restype = c_int
    return libspm.spmSort( spm )

def pyspm_spmMergeDuplicate( spm ):
    libspm.spmMergeDuplicate.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmMergeDuplicate.restype = pastix_int
    return libspm.spmMergeDuplicate( spm )

def pyspm_spmSymmetrize( spm ):
    libspm.spmSymmetrize.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmSymmetrize.restype = pastix_int
    return libspm.spmSymmetrize( spm )

def pyspm_spmCheckAndCorrect( spm ):
    libspm.spmCheckAndCorrect.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmCheckAndCorrect.restype = POINTER(pypastix_spm_t)
    return libspm.spmCheckAndCorrect( spm )

def pyspm_spmGenRHS( type, nrhs, spm, x, ldx, b, ldb ):
    libspm.spmGenRHS.argtypes = [ c_int, pastix_int, POINTER(pypastix_spm_t),
                                  c_void_p, pastix_int, c_void_p, pastix_int ]
    libspm.spmGenRHS.restype = c_int
    return libspm.spmGenRHS( type, nrhs, spm, x, ldx, b, ldb )

def pyspm_spmCheckAxb( nrhs, spm, x0, ldx0, b, ldb, x, ldx ):
    libspm.spmCheckAxb.argtypes = [ pastix_int, POINTER(pypastix_spm_t),
                                    c_void_p, pastix_int, c_void_p, pastix_int,
                                    c_void_p, pastix_int ]
    libspm.spmCheckAxb.restype = c_int
    return libspm.spmCheckAxb( nrhs, spm, x0, ldx0, b, ldb, x, ldx )

def pyspm_spmIntConvert( n, input ):
    libspm.spmIntConvert.argtypes = [ pastix_int, c_int_p ]
    libspm.spmIntConvert.restype = POINTER(pastix_int)
    return libspm.spmIntConvert( n, input )

def pyspm_spmLoad( spm ):
    libspm.spmLoad.argtypes = [ POINTER(pypastix_spm_t), c_void_p ]
    libspm.spmLoad.restype = c_int
    return libspm.spmLoad( spm, None )

def pyspm_spmSave( spm ):
    libspm.spmSave.argtypes = [ POINTER(pypastix_spm_t), c_void_p ]
    libspm.spmSave.restype = c_int
    return libspm.spmSave( spm, None )

def pyspm_spmReadDriver( driver, filename, spm, pastix_comm ):
    libspm.spmReadDriver.argtypes = [ c_int, c_char_p, POINTER(pypastix_spm_t),
                                      c_int ]
    libspm.spmReadDriver.restype = c_int
    return libspm.spmReadDriver( driver, filename, spm, pastix_comm )

def pyspm_spm2Dense( spm ):
    libspm.spm2Dense.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spm2Dense( spm )

def pyspm_spmPrint( spm ):
    libspm.spmPrint.argtypes = [ POINTER(pypastix_spm_t), c_void_p ]
    libspm.spmPrint( spm, None )

def pyspm_spmPrintInfo( spm ):
    libspm.spmPrintInfo.argtypes = [ POINTER(pypastix_spm_t), c_void_p ]
    libspm.spmPrintInfo( spm, None )

def pyspm_spmExpand( spm ):
    libspm.spmExpand.argtypes = [ POINTER(pypastix_spm_t) ]
    libspm.spmExpand.restype = POINTER(pypastix_spm_t)
    return libspm.spmExpand( spm )

def pyspm_spmDofExtend( spm, type, dof ):
    libspm.spmDofExtend.argtypes = [ POINTER(pypastix_spm_t), c_int, c_int ]
    libspm.spmDofExtend.restype = POINTER(pypastix_spm_t)
    return libspm.spmDofExtend( spm, type, dof )


