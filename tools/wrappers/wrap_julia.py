#!/usr/bin/env python
"""
Wrapper Julia
==============

 @file wrappers/wrap_julia.py

 PaStiX generator for the  wrapper

 @copyright 2019-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.1
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date 2023-07-21

"""
import os
import re
import argparse
import time
from . import *
from .pastix_julia import *

def function_prepare_arg( function, arg, return_value ):
    """Generate a declaration for a variable in the interface."""

    jtype = types_dict[arg['type']]
    if arg['pointer'] > 0:
        if jtype == "Cvoid":
            jtype = "Ptr{Cvoid}"
        elif jtype == "Cchar":
            jtype = "Cstring"
        elif jtype == "Cint":
            jtype = "Ptr{Cint}"
        else:
            jtype = "Ptr{"+jtype+"}"
    if arg['pointer'] > 1:
        jtype = "Ptr{Cvoid}"

    arg['jtype'] = jtype
    arg['jname'] = format("%s::" % arg['name'] )

    # Update the maximum length
    sizes = function['sizes']
    sizes['type'] = max( sizes['type'], len(arg['jtype']) )
    sizes['name'] = max( sizes['name'], len(arg['jname']) )

def iso_c_interface_type(arg, return_value, args_list, args_size):
    """Generate a declaration for a variable in the interface."""

    is_double_ptr = False
    if (arg[1] == "*" or arg[1] == "**" ):
        is_pointer = True
    elif (arg[1] == "**"):
        is_double_ptr = True
    else:
        is_pointer = False

    f_type = types_dict[arg[0]]
    if is_pointer:
        if f_type == "Cvoid":
            f_type = "Ptr{Cvoid}"
        elif f_type == "Cchar":
            f_type = "Cstring"
        elif f_type == "Cint":
            f_type = "Ptr{Cint}"
        else:
            f_type = "Ptr{"+f_type+"}"
    if  is_double_ptr:
        f_type = "Ptr{" + f_type + "}"

    if (not return_value and arg[1] != "**"):
        f_pointer = "value"
    else:
        f_pointer = ""

    f_name = format("%s::" % arg[2] )

    args_size[0] = max(args_size[0], len(f_name))
    args_size[1] = max(args_size[1], len(f_type))
    args_list.append( [ f_name, f_type ] );

class wrap_julia:

    @staticmethod
    def write_header( f ):
        filename = os.path.basename( f['filename'] )
        filename = re.sub(r"\.in", "", filename)
        header = '''#=

 @file ''' + filename + '''

 ''' + f['description'] + '''

 @copyright 2020-''' + time.strftime( "%Y" ) + ''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.1
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date ''' + time.strftime( "%Y-%m-%d" ) + '''

 This file has been automatically generated with gen_wrappers.py

 @ingroup wrap_julia

=#
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
        bib = ""
        Bib = ""
        BIB = ""
        if ("SPM" in f['description']):
            bib = "spm"
            Bib = "Spm"
            BIB = "SPM"
        elif ("PaStiX" in f['description']):
            bib = "pastix"
            Bib = "Pastix"
            BIB = "PASTIX"
        jl_interface = "@cenum " + Bib + "_" + ename + "_t " + "{\n"

        # loop over the arguments of the enum to get max param length
        # And modify the names first
        length=0
        for param in params:
            # Remove Pastix from everything
            param[0] = re.sub(r"Pastix", "", param[0])
            param[0] = re.sub(r"PASTIX_", "", param[0])

            if ename == "mtxtype":
                param[1] = re.sub(r"trans.", "Spm", param[1])
                param[1] = re.sub(r"Pastix", "", param[1])
                param[1] = param[1].lower()
            # if ename == "verbose":
            #     param[0] = re.sub(r"Verbose", "", param[0])
            length = max( length, len(param[0]))
        fmt="%-"+ str(length) + "s"

        # Increment for index array enums
        inc = 0
        if ename[1:5] == "parm":
            inc=1

        # loop over the arguments of the enum
        suffix=""
        for param in params:
            name = param[0].lower()
            if isinstance(param[1],int):
                if name[1:10] == "parm_size":
                    value = str(param[1])
                else:
                    value = str(param[1] + inc)
            else:
                value = str(param[1])
            if(ename == "error" and  name=="SUCCESS"):
                jl_interface += indent + format(fmt % name) + indent + " = " + value + ",\n"
            else :
                jl_interface += indent + format(fmt % (name + suffix)) + " = " + value + ",\n"

        jl_interface+="}\n\n"
        return jl_interface

    @staticmethod
    def write_struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        py_interface = ""

        s = 0
        name = struct[0][2]
        name = re.sub(r"pastix_", "", name)
        py_interface +=  "@cstruct " + "Pastix" + "_" + name + " {\n"
        s = iindent
        py_interface +=  s*" "
        headline = s*" "

        slist = []
        ssize = [ 0, 0 ]

        # loop over the arguments of the enum
        for j in range(1,len(struct)):
            iso_c_interface_type(struct[j], True, slist, ssize)

        s += iindent

        # loop over the arguments of the struct
        for j in range(0,len(slist)):
            if (j > 0):
                py_interface += "\n" + headline
            py_interface += format(slist[j][0] + slist[j][1])

        py_interface += "\n}\n\n"
        return py_interface

    @staticmethod
    def write_function( function ):
        """Generate an interface for a function."""

        return_type    = function['rettype']['type']
        return_pointer = function['rettype']['pointer']

        sizes = { 'type' : 0, 'name' : 0 }
        function['sizes']  = sizes
        for arg in function['args']:
            function_prepare_arg( function, arg, False )

        c_symbol = function['name']
        if "pastix" in c_symbol:
            libname = "libpastix"
            prefix  = ""
        elif "spm" in c_symbol:
            libname = "libspm"
            prefix  = ""
        else:
            print("ERROR: function name without pastix nor spm")
            return
        cbinding_line = "@cbindings " + libname + " begin\n"
        func_line = indent + "@cextern " + c_symbol + "( "

        # Print the argument of the function
        j = 0
        for arg in function['args']:
            if j >= 1:
                func_line += ", "
            func_line += arg['jname'] + arg['jtype']
            j += 1

        # Print the return type
        ret_val = types_dict[return_type]
        if return_pointer == 1:
            ret_val = "Ptr{"+ret_val+"}"
        elif return_pointer > 1:
            ret_val = "Ptr{Ptr{"+ret_val+"}}"
        func_line += " )::" + ret_val

        # Final print
        str_function = cbinding_line + func_line + "\nend\n\n"
        return str_function

    @staticmethod
    def write_file( data, enum_list, struct_list, function_list ):
        """
        Generate a single julia file. It will contains:
        enums, structs and interfaces of all C functions
        """

        modulefile = open( data['filename'], "w" )

        header_str = wrap_julia.write_header( data )
        modulefile.write( header_str )

        # enums
        if (enum_list and len(enum_list) > 0):
            for enum in enum_list:
                enum_cpy = gen_enum_copy( enum )
                enum_str = wrap_julia.write_enum( data, enum_cpy )
                modulefile.write( enum_str )

        # derived types
        if (struct_list and len(struct_list) > 0):
            for struct in struct_list:
                struct_str = wrap_julia.write_struct( struct )
                modulefile.write( struct_str )

        # functions
        if (function_list and len(function_list) > 0):
            for function in function_list:
                function_str = wrap_julia.write_function( function )
                modulefile.write( function_str )

        footer_str = wrap_julia.write_footer( data )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write( enum_list, struct_list, function_list ):
        f = wrap_julia.write_file( enums, enum_list, struct_list, None )
        print( "Exported file: " + f )

        f = wrap_julia.write_file( common, None, None, function_list )
        print( "Exported file: " + f )
