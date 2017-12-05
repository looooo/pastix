#!/usr/bin/env python
import os
import re
import argparse

indent="    "

class wrap_python:

    @staticmethod
    def header( f ):
        modname = os.path.basename( f['filename'] )
        modname = re.sub(r"\.in", "", modname)
        header = '''
"""

 @file '''+modname+'''

 PaStiX python definition of enums and datatypes

 @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-05-04

This file has been automatically generated with gen_wrappers.py

"""
import numpy as np
import ctypes

''' + f['header'] + "\n"
        return header;

    @staticmethod
    def footer( f ):
        filename = os.path.basename( f['filename'] )
        modname = re.sub(r".f90", "", filename, flags=re.IGNORECASE)
        footer = f['footer']
        return footer

    @staticmethod
    def enum( f, enum ):
        """Generate an interface for an enum.
           Translate it into constants."""

        ename  = enum[0]
        params = enum[1]

        # initialize a string with the fortran interface
        py_interface = "class " + ename + ":\n"

        # loop over the arguments of the enum to get max param length
        length=0
        for param in params:
            # Convert IPARM/DPARM to lower case
            if param[0][1:5] == "PARM":
                param[0] = re.sub(r"[ID]PARM_", "", param[0])
                param[0] = param[0].lower()

            # Remove Pastix from everything
            param[0] = re.sub(r"Pastix", "", param[0])

            if ename == "error":
                param[0] = re.sub(r"PASTIX_", "", param[0])
                param[0] = re.sub(r"ERR_", "", param[0])
            elif ename == "fact_mode" or ename == "factotype" or ename == "solv_mode":
                param[0] = re.sub(r"Fact", "", param[0])
                param[0] = re.sub(r"Solv", "", param[0])
                param[0] = re.sub(r"Mode", "", param[0])
            elif ename == "scheduler":
                param[0] = re.sub(r"Sched", "", param[0])
            elif ename == "threadmode":
                param[0] = re.sub(r"Thread", "", param[0])
            elif ename[0:8] == "compress":
                param[0] = re.sub(r"Compress", "", param[0])
                param[0] = re.sub(r"When", "", param[0])
                param[0] = re.sub(r"Method", "", param[0])
            elif ename == "rhstype":
                param[0] = re.sub(r"Rhs", "", param[0])
            elif ename == "trans":
                param[0] = param[0]
            elif ename == "mtxtype":
                param[1] = re.sub(r"Pastix", "trans.", param[1])
            elif ename == "normtype":
                param[0] = re.sub(r"Norm", "", param[0])
            else:
                param[0] = re.sub(ename, "", param[0], flags=re.IGNORECASE)
            length = max( length, len(param[0]))
        fmt="%-"+ str(length) + "s"

        # loop over the arguments of the enum
        for param in params:
            name  = param[0]
            value = str(param[1])

            py_interface += indent + format(fmt % name) + " = " + value + "\n"

        if ename in f['enums']:
            py_interface += f['enums'][ename]

        return py_interface

    @staticmethod
    def struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        f_interface = ""

        name = struct[0][2]
        f_interface += indent + "type, bind(c) :: " + name + "\n"
        # loop over the arguments of the enum
        for j in range(1,len(struct)):
            f_interface += indent + indent + iso_c_interface_type(struct[j], True)
            f_interface += "\n"

        f_interface += indent + "end type " + name + "\n"

        return f_interface

    @staticmethod
    def function(function):
        """Generate an interface for a function."""

        # is it a function or a subroutine
        if (function[0][0] == "void"):
            is_function = False
        else:
            is_function = True

        c_symbol = function[0][2]
        f_symbol = c_symbol + "_c"

        used_derived_types = set([])
        for arg in function:
            type_name = arg[0]
            if (type_name in derived_types):
                used_derived_types.add(type_name)

        # initialize a string with the fortran interface
        f_interface = ""
        f_interface += indent + "interface\n"

        if (is_function):
            f_interface += indent + indent + "function "
        else:
            f_interface += indent + indent + "subroutine "

        f_interface += f_symbol + "("

        if (is_function):
            initial_indent = len(indent + indent + "function " + f_symbol + "(") * " "
        else:
            initial_indent = len(indent + indent + "subroutine " + f_symbol + "(") * " "

        # loop over the arguments to compose the first line
        for j in range(1,len(function)):
            if (j != 1):
                f_interface += ", "
            if (j%9 == 0):
                f_interface += "&\n" + initial_indent

            f_interface += function[j][2]

        f_interface += ") &\n"
        f_interface += indent + indent + "  " + "bind(c, name='" + c_symbol +"')\n"

        # add common header
        f_interface += indent + 2*indent + "use iso_c_binding\n"
        # import derived types
        for derived_type in used_derived_types:
            f_interface += indent + 2*indent + "import " + derived_type +"\n"
        f_interface += indent + 2*indent + "implicit none\n"


        # add the return value of the function
        if (is_function):
            f_interface +=  indent + 2*indent + iso_c_interface_type(function[0], True) + "_c"
            f_interface += "\n"

        # loop over the arguments to describe them
        for j in range(1,len(function)):
            f_interface += indent + 2*indent + iso_c_interface_type(function[j], False)
            f_interface += "\n"

        if (is_function):
            f_interface += indent + indent + "end function\n"
        else:
            f_interface += indent + indent + "end subroutine\n"

        f_interface += indent + "end interface\n"

        return f_interface

    @staticmethod
    def wrapper(function):
        """Generate a wrapper for a function.
           void functions in C will be called as subroutines,
           functions in C will be turned to subroutines by appending
           the return value as the last argument."""

        # is it a function or a subroutine
        if (function[0][0] == "void"):
            is_function = False
        else:
            is_function = True

        c_symbol = function[0][2]
        f_symbol = c_symbol + "_c"

        if (is_function):
            initial_indent_signature = len(indent + "subroutine " + c_symbol + "(") * " "
            initial_indent_call      = len(indent + indent + "info = " + f_symbol + "(") * " "
        else:
            initial_indent_signature = len(indent + "subroutine " + c_symbol + "(") * " "
            initial_indent_call      = len(indent + indent + "call " + f_symbol + "(") * " "

        # loop over the arguments to compose the first line and call line
        signature_line = ""
        call_line = ""
        double_pointers = []
        for j in range(1,len(function)):
            if (j != 1):
                signature_line += ", "
                call_line += ", "

            # do not make the argument list too long
            if (j%9 == 0):
                call_line      += "&\n" + initial_indent_call
                signature_line += "&\n" + initial_indent_signature

            # pointers
            arg_type    = function[j][0]
            arg_pointer = function[j][1]
            arg_name    = function[j][2]

            signature_line += arg_name
            if (arg_pointer == "**"):
                aux_name = arg_name + "_aux"
                call_line += aux_name
                double_pointers.append(arg_name)
            elif (arg_pointer == "*"):
                call_line += "c_loc(" + arg_name + ")"
            else:
                call_line += arg_name

        contains_derived_types = False
        for arg in function:
            if (arg[0] in derived_types):
                contains_derived_types = True

        # initialize a string with the fortran interface
        f_wrapper = ""
        f_wrapper += indent + "subroutine "
        f_wrapper += c_symbol + "("

        # add the info argument at the end
        f_wrapper += signature_line
        if (is_function):
            if (len(function) > 1):
                f_wrapper += ", "

            return_type = function[0][0]
            return_pointer = function[0][1]
            if (return_type == "int"):
                return_var = "info"
            else:
                return_var = return_variables_dict[return_type]

            f_wrapper += return_var

        f_wrapper += ")\n"

        # add common header
        f_wrapper += indent + indent + "use iso_c_binding\n"
        f_wrapper += indent + indent + "implicit none\n"

        # loop over the arguments to describe them
        for j in range(1,len(function)):
            f_wrapper += indent + indent + iso_c_wrapper_type(function[j]) + "\n"

        # add the return info value of the function
        if (is_function):
            if (function[0][1] == "*"):
                f_target = ", pointer"
            else:
                f_target = ""

            f_wrapper += indent + indent + types_dict[return_type] + ", intent(out)" + f_target + " :: " + return_var + "\n"

        f_wrapper += "\n"

        # loop over potential double pointers and generate auxiliary variables for them
        for double_pointer in double_pointers:
            aux_name = double_pointer + "_aux"
            f_wrapper += indent + indent + "type(c_ptr) :: " + aux_name + "\n"
            f_wrapper += "\n"

        if (is_function):
            f_return = return_var
            f_return += " = "
        else:
            f_return = "call "

        # generate the call to the C function
        if (is_function and return_pointer == "*"):
            f_wrapper += indent + indent + "call c_f_pointer(" + f_symbol + "(" + call_line + "), " + return_var + ")\n"
        else:
            f_wrapper += indent + indent + f_return + f_symbol + "(" + call_line + ")\n"

        # loop over potential double pointers and translate them to Fortran pointers
        for double_pointer in double_pointers:
            aux_name = double_pointer + "_aux"
            f_wrapper += indent + indent + "call c_f_pointer(" + aux_name + ", " + double_pointer + ")\n"

        f_wrapper += indent + "end subroutine\n"

        return f_wrapper
