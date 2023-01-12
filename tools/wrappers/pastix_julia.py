#!/usr/bin/env python
"""
Wrapper Julia
=============

 @file wrappers/pastix_julia.py

 PaStiX Julia wrapper variables

 @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.0
 @author Mathieu Faverge
 @author Tony Delarue
 @date 2022-09-27

"""
filename_prefix = "wrappers/julia/PaStiX/src/"

enums = {
    'filename'    : filename_prefix + 'pastix_enums.jl.in',
    'description' : "PaStiX julia wrapper to define enums and datatypes",
    'header'      : """
Pastix_int_t = @PASTIX_JULIA_INTEGER@
pastix_mpi_enabled = @PASTIX_JULIA_MPI_ENABLED@
const Pastix_data_t = Cvoid
const Pastix_rhs_t = Ptr{Cvoid}
const Pastix_graph_t = Ptr{Cvoid}

using spm
if pastix_mpi_enabled
    using MPI
end

function __get_mpi_type__()
    if !pastix_mpi_enabled
        return Cint
    elseif sizeof(MPI.MPI_Comm) == sizeof(Clong)
        return Clong
    elseif sizeof(MPI.MPI_Comm) == sizeof(Cint)
        return Cint
    end
    return Cvoid
end

""",
    'footer'      : "",
    'enums'       : {}
}

common = {
    'filename'    : filename_prefix + 'PaStiX.jl',
    'description' : "PaStiX julia wrapper",
    'header'      : """
module PaStiX
using CBinding
using Libdl
include(\"pastix_enums.jl\")

function pastix_library_path()
    x = Libdl.dlext
    return \"libpastix.$x\"
end
libpastix = pastix_library_path()

""",
    'footer'      : "end #module\n",
    'enums'       : {}
}

# set indentation in the python file
indent="    "
iindent=4

# translation_table of types
types_dict = {
    "int":            ("Cint"),
    "int8_t":         ("Int8"),
    "seed_t":                 ("Culonglong"),
    "unsigned long long int": ("Culonglong"),
    "spm_coeftype_t": ("spm.spm_coeftype_t"),
    "spm_dir_t":      ("spm.spm_dir_t"),
    "spm_trans_t":    ("spm.spm_trans_t"),
    "spm_uplo_t":     ("spm.spm_uplo_t"),
    "spm_diag_t":     ("spm.spm_diag_t"),
    "spm_side_t":     ("spm.spm_side_t"),
    "spm_driver_t":   ("spm.spm_driver_t"),
    "spm_fmttype_t":  ("spm.spm_fmttype_t"),
    "spm_layout_t":   ("spm.spm_layout_t"),
    "spm_normtype_t": ("spm.spm_normtype_t"),
    "spm_rhstype_t":  ("spm.spm_rhstype_t"),
    "spm_mtxtype_t":  ("spm.spm_mtxtype_t"),
    "spm_int_t":      ("spm.spm_int_t"),
    "spmatrix_t":     ("spm.spmatrix_t"),
    "size_t":         ("Csize_t"),
    "char":           ("Cchar"),
    "double":         ("Cdouble"),
    "float":          ("Cfloat"),
    "spm_complex64_t":("ComplexF64"),
    "spm_complex32_t":("ComplexF32"),
    "void":           ("Cvoid"),
    "MPI_Comm":       ("__get_mpi_type__()"),
    "FILE":           ("Cvoid"),
    "pastix_coeftype_t": ("spm.spm_coeftype_t"),
    "pastix_dir_t":      ("spm.spm_dir_t"),
    "pastix_trans_t":    ("Pastix_trans_t"),
    "pastix_uplo_t":     ("Cint"),
    "pastix_diag_t":     ("Pastix_diag_t"),
    "pastix_side_t":     ("Cint"),
    "pastix_driver_t":   ("Pastix_driver_t"),
    "pastix_fmttype_t":  ("Pastix_fmttype_t"),
    "pastix_layout_t":   ("Pastix_layout_t"),
    "pastix_normtype_t": ("Pastix_normtype_t"),
    "pastix_rhstype_t":  ("Pastix_rhstype_t"),
    "pastix_mtxtype_t":  ("Pastix_mtxtype_t"),
    "pastix_int_t":      ("Pastix_int_t"),
    "pastix_data_t":     ("Pastix_data_t"),
    "pastix_rhs_t":      ("Pastix_rhs_t"),
    "pastix_ordering_t": ("Pastix_ordering_t"),
    "pastix_order_t":    ("Pastix_order_t"),
    "pastix_graph_t":    ("Pastix_graph_t"),
    "PASTIX_Comm":       ("__get_mpi_type__()"),
}
