#!/usr/bin/env python
"""
Wrapper Python
==============

 @file wrappers/pastix_python.py

 PaStiX Python wrapper variables

 @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.0
 @author Mathieu Faverge
 @author Tony Delarue
 @date 2022-09-27

"""
filename_prefix = "wrappers/python/examples/pypastix/"

enums_python_coeftype='''
    @staticmethod
    def getptype ( dtype ):
        np_dict = {
            np.dtype('float32')    : coeftype.Float,
            np.dtype('float64')    : coeftype.Double,
            np.dtype('complex64')  : coeftype.Complex32,
            np.dtype('complex128') : coeftype.Complex64,
        }
        if dtype in np_dict:
            return np_dict[dtype]
        else:
            return -1

    @staticmethod
    def getctype ( flttype ):
        class c_float_complex(Structure):
            _fields_ = [("r",c_float),("i", c_float)]
        class c_double_complex(Structure):
            _fields_ = [("r",c_double),("i", c_double)]

        np_dict = {
            coeftype.Float     : c_float,
            coeftype.Double    : c_double,
            coeftype.Complex32 : c_float_complex,
            coeftype.Complex64 : c_double_complex
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1

    @staticmethod
    def getnptype ( flttype ):
        np_dict = {
            coeftype.Float     : np.dtype('float32'),
            coeftype.Double    : np.dtype('float64'),
            coeftype.Complex32 : np.dtype('complex64'),
            coeftype.Complex64 : np.dtype('complex128')
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1
'''

enums = {
    'filename'    : filename_prefix + 'enum.py.in',
    'description' : "PaStiX python wrapper to define enums and datatypes",
    'header'      : """
# Start with __ to prevent broadcast to file importing enum
__pastix_int__ = @PASTIX_PYTHON_INTEGER@
__pastix_mpi_enabled__ = @PASTIX_PYTHON_MPI_ENABLED@

""",
    'footer'      : "",
    'enums'       : { 'coeftype' : enums_python_coeftype,
                      'mtxtype'  : "    SymPosDef = trans.ConjTrans + 1\n    HerPosDef = trans.ConjTrans + 2\n" }
}

common = {
    'filename'    : filename_prefix + '__pastix__.py',
    'description' : "PaStiX python wrapper",
    'header'      : """
from . import libpastix
from .enum import __pastix_int__
from .enum import __pastix_mpi_enabled__
from spm import pyspm_spmatrix_t
import spm

if __pastix_mpi_enabled__:
    from mpi4py import MPI
    if MPI._sizeof(MPI.Comm) == sizeof(c_long):
        pypastix_mpi_comm = c_long
    elif MPI._sizeof(MPI.Comm) == sizeof(c_int):
        pypastix_mpi_comm = c_int
    else:
        pypastix_mpi_comm = c_void_p

    pypastix_default_comm = MPI.COMM_WORLD

    def pypastix_convert_comm( comm ):
        comm_ptr = MPI._addressof(comm)
        return pypastix_mpi_comm.from_address(comm_ptr)
else:
    pypastix_mpi_comm = c_int

    pypastix_default_comm = 0

    def pypastix_convert_comm( comm ):
        return c_int(comm)

""",
    'footer'      : "",
    'enums'       : {}
}

# set indentation in the python file
indent="    "
iindent=4

# translation_table of types
types_dict = {
    "int":            ("c_int"),
    "int8_t":         ("c_int8"),
    "seed_t":                 ("c_ulonglong"),
    "unsigned long long int": ("c_ulonglong"),
    "spm_coeftype_t": ("c_int"),
    "spm_dir_t":      ("c_int"),
    "spm_trans_t":    ("c_int"),
    "spm_uplo_t":     ("c_int"),
    "spm_diag_t":     ("c_int"),
    "spm_side_t":     ("c_int"),
    "spm_driver_t":   ("c_int"),
    "spm_fmttype_t":  ("c_int"),
    "spm_layout_t":   ("c_int"),
    "spm_normtype_t": ("c_int"),
    "spm_rhstype_t":  ("c_int"),
    "spm_mtxtype_t":  ("c_int"),
    "spm_int_t":      ("__spm_int__"),
    "spmatrix_t":     ("pyspm_spmatrix_t"),
    "size_t":         ("c_size_t"),
    "char":           ("c_char"),
    "double":         ("c_double"),
    "float":          ("c_float"),
    "spm_complex64_t":("c_double_complex"),
    "spm_complex32_t":("c_float_complex"),
    "void":           ("c_void"),
    "MPI_Comm":       ("pypastix_mpi_comm"),
    "FILE":           ("c_void"),
    "pastix_coeftype_t": ("c_int"),
    "pastix_dir_t":      ("c_int"),
    "pastix_trans_t":    ("c_int"),
    "pastix_uplo_t":     ("c_int"),
    "pastix_diag_t":     ("c_int"),
    "pastix_side_t":     ("c_int"),
    "pastix_driver_t":   ("c_int"),
    "pastix_fmttype_t":  ("c_int"),
    "pastix_layout_t":   ("c_int"),
    "pastix_normtype_t": ("c_int"),
    "pastix_rhstype_t":  ("c_int"),
    "pastix_mtxtype_t":  ("c_int"),
    "pastix_int_t":      ("__pastix_int__"),
    "pastix_data_t":     ("c_void"),
    "pastix_rhs_t":      ("c_void_p"),
    "pastix_order_t":    ("c_void"),
    "pastix_graph_t":    ("c_void"),
    "PASTIX_Comm":       ("pypastix_mpi_comm"),
}
