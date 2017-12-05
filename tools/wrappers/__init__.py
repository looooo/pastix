"""
Wrappers
========

 @file __init__.py

 PaStiX wrapper generators module intialization

 @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Mathieu Faverge
 @date 2017-05-04

"""

# translation_table of types
types_dict = {
    "int":               ("integer(kind=c_int)"),
    "pastix_coeftype_t": ("integer(c_int)"),
    "pastix_dir_t":      ("integer(c_int)"),
    "pastix_trans_t":    ("integer(c_int)"),
    "pastix_uplo_t":     ("integer(c_int)"),
    "pastix_diag_t":     ("integer(c_int)"),
    "pastix_side_t":     ("integer(c_int)"),
    "pastix_driver_t":   ("integer(c_int)"),
    "pastix_fmttype_t":  ("integer(c_int)"),
    "pastix_layout_t":   ("integer(c_int)"),
    "pastix_normtype_t": ("integer(c_int)"),
    "pastix_rhstype_t":  ("integer(c_int)"),
    "pastix_mtxtype_t":  ("integer(c_int)"),
    "pastix_data_t":     ("type(pastix_data_t)"),
    "pastix_spm_t":      ("type(pastix_spm_t)"),
    "pastix_int_t":      ("integer(kind=pastix_int_t)"),
    "pastix_order_t":    ("type(pastix_order_t)"),
    "size_t":            ("integer(kind=c_size_t)"),
    "char":              ("character(kind=c_char)"),
    "double":            ("real(kind=c_double)"),
    "float":             ("real(kind=c_float)"),
    "pastix_complex64_t":("complex(kind=c_double_complex)"),
    "pastix_complex32_t":("complex(kind=c_float_complex)"),
    "void":              ("type(c_ptr)"),
    "MPI_Comm":          ("integer(kind=c_int)"),
    "FILE":              ("type(c_ptr)"),
}

# translation_table with names of auxiliary variables
return_variables_dict = {
    "double":            ("value"),
    "float":             ("value"),
    "pastix_int_t":      ("value"),
    "pastix_spm_t":      ("spmo"),
}

# global list used to determine derived types
derived_types = ['pastix_int_t']

# name arrays which will be translated to assumed-size arrays, e.g. pA(*)
arrays_names_2D = ["pA", "pB", "pC", "pAB", "pQ", "pX", "pAs"]
arrays_names_1D = ["ipiv", "values", "work"]

__all__ = [ 'types_dict', 'return_variables_dict', 'derived_types', 'arrays_names_1D', 'arrays_names_2D' ]

from .wrap_python  import *
from .wrap_fortran import *
