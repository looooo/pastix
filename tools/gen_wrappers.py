#!/usr/bin/env python
"""
 @file gen_wrappers.py

 Python and Fortran 90 wrapper generator for some of the solverstack
 libraries, inspired from the PLASMA-OMP fortran generator.

 @copyright 2016-2017 University of Tennessee, US, University of
                      Manchester, UK. All rights reserved.
 @copyright 2017-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.2.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date 2021-04-07

"""
import os
import re
import argparse
import wrappers

description = '''\
Generates Fortran 90 and Python wrappers from the spm and pastix header files.'''

help = '''\
----------------------------------------------------------------------
Example uses:

   $PASTIX_SRC_DIR/gen_wrappers.py

----------------------------------------------------------------------
'''

# ------------------------------------------------------------
# command line options
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=description,
    epilog=help )
parser.add_argument('--prefix',        action='store', help='Prefix for variables in Makefile', default='./')
parser.add_argument('args', nargs='*', action='store', help='Files to process')
opts = parser.parse_args()

filenames = [
    "include/pastix/api.h",
    "include/pastix/order.h",
    "include/pastix.h"
]

def main():

    # common cleaned header files
    preprocessed_list = []

    # source header files
    for filename in filenames:
        # source a header file
        c_header_file = open(filename, 'r').read()

        # clean the string (remove comments, macros, etc.)
        clean_file = wrappers.polish_file(c_header_file)

        # convert the string to a list of strings
        initial_list = clean_file.split("\n")

        # process the list so that each enum, struct or function
        # are just one item
        nice_list = wrappers.preprocess_list(initial_list)

        # compose all files into one big list
        preprocessed_list += nice_list

    # register all enums
    enum_list = wrappers.parse_enums( preprocessed_list )

    # register all structs
    struct_list = wrappers.parse_structs( preprocessed_list )

    # register all individual functions and their signatures
    function_list = wrappers.parse_prototypes( preprocessed_list )

    # print( "------------ ENUM ----------------" )
    # print( enum_list )
    # print( "------------ STRUCT ----------------" )
    # print( struct_list )
    # print( "------------ FUNCTION ----------------" )
    # print( function_list )

    for g in [ wrappers.wrap_fortran, wrappers.wrap_python, wrappers.wrap_julia ]:
        g.write( enum_list, struct_list, function_list )

# execute the program
main()
