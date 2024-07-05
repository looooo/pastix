#!/usr/bin/env python
"""
Wrapper Python
==============

 @file wrappers/wrap_python.py

 PaStiX generator for the python wrapper

 @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.2
 @author Mathieu Faverge
 @author Tony Delarue
 @author Esragul Korkmaz
 @author Selmane Lebdaoui
 @date 2023-07-21

"""
import os
import re
import argparse
import time
from . import *
from .pastix_python import *

def function_prepare_arg( function, arg ):
    """Generate a declaration for a variable in the Fortran wrapper."""

    py_type = types_dict[arg['type']]
    if arg['pointer'] > 0:
        if py_type == "c_void":
            py_type = "c_void_p"
        elif py_type == "c_char":
            py_type = "c_char_p"
        elif py_type == "c_int":
            py_type = "c_int_p"
        else:
            py_type = "POINTER("+py_type+")"

    py_name = arg['name']
    py_call = py_name

    if arg['type'] == 'FILE':
        py_name = ""
        py_call = "None"

    # detect array argument
    if (arg['pointer'] > 0) and (arg['type'] != "void"):
        vname = re.sub(r"^opt_", "", py_name, flags=re.IGNORECASE)
        if vname in arrays_names_2D:
            py_call = py_name + ".ctypes.data_as( " + py_type + " )"
        elif vname in arrays_names_1D:
            py_call = py_name + ".ctypes.data_as( " + py_type + " )"

    if arg['pointer'] > 1:
        py_call = "pointer( " + py_call + " )"

    if arg['pointer'] > 0 and (arg['type'] == "pastix_rhs_t"):
        py_call = "pointer( " + py_call + " )"

    # Call to communicators
    if arg['type'] == "MPI_Comm":
        py_call = "pypastix_convert_comm( " + py_call + " )"

    arg['py_name'] = py_name
    arg['py_type'] = py_type
    arg['py_call'] = py_call

def iso_c_interface_type(arg, return_value, args_list, args_size):
    """Generate a declaration for a variable in the interface."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    f_type = types_dict[arg[0]]
    if is_pointer:
        if f_type == "c_void":
            f_type = "c_void_p"
        elif f_type == "c_char":
            f_type = "c_char_p"
        elif f_type == "c_int":
            f_type = "c_int_p"
        else:
            f_type = "POINTER("+f_type+")"

    if (not return_value and arg[1] != "**"):
        f_pointer = "value"
    else:
        f_pointer = ""

    f_name = format("\"%s\", " % arg[2] )

    args_size[0] = max(args_size[0], len(f_name))
    args_size[1] = max(args_size[1], len(f_type))
    args_list.append( [ f_name, f_type ] );

class wrap_python:

    @staticmethod
    def write_header( f ):
        filename = os.path.basename( f['filename'] )
        filename = re.sub(r"\.in", "", filename)
        header = '''"""

 @file ''' + filename + '''

 ''' + f['description'] + '''

 @copyright 2017-''' + time.strftime( "%Y" ) + ''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.2
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Tony Delarue
 @author Selmane Lebdaoui
 @date ''' + time.strftime( "%Y-%m-%d" ) + '''

 This file has been automatically generated with gen_wrappers.py

 @ingroup wrap_python

"""
from ctypes import *
import numpy as np
'''
        if f['header'] != "":
            header += f['header']

        return header

    @staticmethod
    def write_footer( f ):
        return f['footer']

    @staticmethod
    def write_enum( f, enum ):
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
            param[0] = re.sub(r"Spm", "", param[0])

            if ename == "error":
                param[0] = re.sub(r"PASTIX_", "", param[0])
                param[0] = re.sub(r"SPM_", "", param[0])
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
                param[1] = re.sub(r"Spm", "trans.", param[1])
            elif ename == "normtype":
                param[0] = re.sub(r"Norm", "", param[0])
            elif ename == "ordering":
                param[0] = re.sub(r"Order", "", param[0])
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

        py_interface += "\n"
        return py_interface

    @staticmethod
    def write_struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        py_interface = ""

        s = 0
        name = struct[0][2]
        name = re.sub(r"pastix_", "", name)
        py_interface +=  "class pypastix_" + name + "(Structure):\n"

        s = iindent
        py_interface +=  s*" " + "_fields_ = ["
        headline = (s+len("_fields_ = ["))*" "

        slist = []
        ssize = [ 0, 0 ]

        # loop over the arguments of the enum
        for j in range(1,len(struct)):
            iso_c_interface_type(struct[j], True, slist, ssize)

        s += iindent
        fmt = "(%-"+ str(ssize[0]) + "s %-"+ str(ssize[1]) +"s)"

        # loop over the arguments of the struct
        for j in range(0,len(slist)):
            if (j > 0):
                py_interface += ",\n" + headline

            py_interface += format( fmt % (slist[j][0], slist[j][1]) )

        py_interface += " ]\n\n"

        return py_interface

    @staticmethod
    def write_function( function ):
        """Generate an interface for a function."""

        return_type    = function['rettype']['type']
        return_pointer = function['rettype']['pointer']

        c_symbol = function['name']
        if "pastix" in c_symbol:
            libname = "libpastix"
            prefix  = "pypastix_"
        elif "spm" in c_symbol:
            libname = "libspm"
            prefix  = "pyspm_"
        else:
            print("ERROR: function name without pastix nor spm", c_symbol)
            return

        # Function declaration
        func_line = "def " + prefix + c_symbol + "("
        func_line_length = len(func_line)
        func_line_indent = len(func_line)

        # Arguments line
        args_line = iindent*" " + libname + "." + c_symbol + ".argtypes = ["
        args_line_length = len(args_line)
        args_line_indent = len(args_line)

        if function['is_function']:
            call_line = iindent*" " + "return " + libname + "." + c_symbol + "("

            py_type = types_dict[return_type]
            if return_pointer == 1:
                if py_type == "c_void":
                    py_type = "c_void_p"
                else:
                    py_type = "POINTER("+py_type+")"
            elif return_pointer > 1:
                py_type = "c_void_p"

            retv_line = iindent*" " + libname + "." + c_symbol + ".restype = " + py_type
        else:
            call_line = iindent*" " + libname + "." + c_symbol + "("
            retv_line = ""
        call_line_length = len(call_line)
        call_line_indent = len(call_line)

        j = 0
        for arg in function['args']:
            function_prepare_arg( function, arg )

            isfile = (arg['type'] == "FILE")
            if j > 0:
                if not isfile:
                    func_line += ","
                    func_line_length += 2
                args_line += ","
                args_line_length += 2
                call_line += ","
                call_line_length += 2

            # func_line
            l = len(arg['py_name'])
            if not isfile:
                if ((func_line_length + l) > 78):
                    func_line_length = func_line_indent
                    func_line += "\n" + func_line_indent*" "
                func_line += " " + arg['py_name']
                func_line_length += l

            # args_line
            l = len(arg['py_type'])
            if ((args_line_length + l) > 78):
                args_line_length = args_line_indent
                args_line += "\n" + args_line_indent*" "
            args_line += " " + arg['py_type']
            args_line_length += l

            # call_line
            l = len(arg['py_call'])
            if ((call_line_length + l) > 78):
                call_line_length = call_line_indent
                call_line += "\n" + call_line_indent*" "
            call_line += " " + arg['py_call']
            call_line_length += l

            j += 1

        str_function  = func_line + " ):\n"
        str_function += args_line + " ]\n"
        if len(retv_line) > 0:
            str_function += retv_line + "\n"
        str_function += call_line + " )\n\n"

        return str_function

    @staticmethod
    def write_file( data, enum_list, struct_list, function_list ):
        """
        Generate a single python file. It will contains:
        enums, structs and interfaces of all C functions
        """

        modulefile = open( data['filename'], "w" )

        header_str = wrap_python.write_header( data )
        modulefile.write( header_str )

        # enums
        if (enum_list and len(enum_list) > 0):
            for enum in enum_list:
                enum_cpy = gen_enum_copy( enum )
                enum_str = wrap_python.write_enum( data, enum_cpy )
                modulefile.write( enum_str )

        # derived types
        if (struct_list and len(struct_list) > 0):
            for struct in struct_list:
                struct_str = wrap_python.write_struct( struct )
                modulefile.write( struct_str )

        # functions
        if (function_list and len(function_list) > 0):
            for function in function_list:
                function_str = wrap_python.write_function( function )
                modulefile.write( function_str )

        footer_str = wrap_python.write_footer( data )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write( enum_list, struct_list, function_list ):
        f = wrap_python.write_file( enums, enum_list, None, None )
        print( "Exported file: " + f )

        f = wrap_python.write_file( common, None, struct_list, function_list )
        print( "Exported file: " + f )
