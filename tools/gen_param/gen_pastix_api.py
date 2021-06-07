import os
import sys
import re
import string
import time

const_str = ''' * This file is generated automatically. If you want to modify it, modify
 * ${PASTIX_HOME}/docs/pastix_params.yaml and run ${PASTIX_HOME}/tools/gen_param/gen_parm_files.py.
 *
 * @copyright 2004-'''+ time.strftime( "%Y" ) +''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @date  '''+ time.strftime( "%Y-%m-%d" )

pastix_enums_begin = '''/**
 *
 * @file pastix/api.h
 *
 * PaStiX API enums parameters.
 *
''' + const_str + '''
 *
 * @addtogroup pastix_api
 * @{
 *
 **/
#ifndef _pastix_api_h_
#define _pastix_api_h_

BEGIN_C_DECLS

'''

close_bracket = '''/**
 * @}
 */
'''

pastix_enums_end = '''
END_C_DECLS

#endif /* _pastix_api_h_ */

''' + close_bracket

blas_str = '''/**
 *
 * @name Constants compatible with CBLAS & LAPACK & PLASMA
 * @{
 *    The naming and numbering of the following constants is consistent with:
 *
 *       - CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz)
 *       - C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/)
 *       - Plasma (http://icl.cs.utk.edu/plasma/index.html)
 *
 */
'''

def get_param_max_str_size( params ) :
    """
    Get the max size of all names, brief description and default value
    (Output purpose).

    @in  params : The params dictionnary
    @out max name, brief and default value size
    """
    maxName   = max( list( map(lambda x : len(x),            params.keys())   ) )

    maxBrief  = max( list( map(lambda x : len(x['BRIEF']),   params.values()) ) )

    maxDefVal = max( list( map(lambda x : len(x['DEFAULT']), params.values()) ) )

    return maxName, maxBrief, maxDefVal

def gen_param( params, prefix ) :
    """
    Generate a string that corresponds to the declaration of the prefix

    @in  params : The params dictionnary
    @in  prefix : IPARM/DPARM
    @out the enum declaration string for the pastix_[i/d]parm_t
    """
    ptype = ""
    if prefix == "IPARM" :
        ptype = "Integer"
    else :
        ptype = "Float"

    result = "/**\n"\
             " * @brief " + ptype + " parameters\n"\
             " */\n"\
             "typedef enum pastix_"+ prefix.lower() +"_e {\n"

    maxName, maxBrief, maxDefVal = get_param_max_str_size( params )

    for field in params :
        param = params[field]
        # Get the fields
        name    = field + ","
        brief   = param['BRIEF']
        default = param['DEFAULT']
        inout   = " " + param['INOUT']

        # Align them correctly
        name    += " " * ( maxName   - len(name) + 1 )
        brief   += " " * ( maxBrief  - len(brief)    )
        default += " " * ( maxDefVal - len(default)  )
        inout   += " " * ( 3         - len(inout) )

        # Add them in the result
        result += "    " + name +" /**< "+ brief +" Default: "+ default + inout +" */\n"

    result += "    "+ prefix +"_SIZE\n"\
              "} pastix_"+ prefix.lower() +"_t;\n\n"
    return result

#
# Wrappers for iparm/dparm
#
def genIparmDeclaration( params ) :
    return gen_param( params, "IPARM" )

def genDparmDeclaration( params ) :
    return gen_param( params, "DPARM" )

def gen_enum_documentation(enum):
    result = "/**\n"\
          " * @brief " + enum['DOC']['BRIEF'] + "\n"

    if 'DETAILS' in enum['DOC'] :
        result += " *\n"

        details = enum['DOC']['DETAILS'].split("\n")
        details = details[:len(details) - 1]

        for line in details:
            if line == '' :
                result += " *\n"
                continue
            result += " * " + line + "\n"
    result += " */\n"

    return result

def gen_enum_declaration( name, enum ) :
    """
    Generate a string that corresponds to the declaration of an enum

    @in  enum : The YAML structure that corresponds to the enum to declare
    @out the enum declaration string for the pastix_[enum]_t
    """
    # Write documentation
    result = gen_enum_documentation( enum )

    # End comment, begin enum declaration
    result += "typedef enum pastix_"+ name + "_e {\n"\

    values = enum['VALUES']
    last   = list(enum['VALUES'].keys())[-1]

    nameMaxSize  = max(list( map(lambda x : len(x['NAME']), values.values()) ))
    briefMaxSize = 0
    valueMaxSize = 0
    if 'BRIEF' in values[0] :
        briefMaxSize = max(list( map(lambda x : len(x['BRIEF']), values.values() ) ))
    if 'VALUE' in values[0] :
        valueMaxSize = max(list( map(lambda x : len(x['VALUE']), values.values() ) ))

    for value in values.values() :
        # Name
        result += "    " + value['NAME']

        # If a value is given, put it in the file
        if 'VALUE' in value :
            result += " " * (nameMaxSize - len(value['NAME']) )
            result += " = "+ value['VALUE']

        if value != enum['VALUES'][ last ] :
            result += ","
        else :
            if 'BRIEF' in value :
                result += " "

        # If a brief comment is given, put it in the file
        if 'BRIEF' in value :
            if 'VALUE' in value :
                result += " " * (valueMaxSize - len(value['VALUE']) ) + " /**< "
            else :
                result += " " * (nameMaxSize  - len(value['NAME']) ) + " /**< "
            result += value['BRIEF'] + " " * (briefMaxSize - len(value['BRIEF']) ) + " */"
        result += "\n"

    # End of the enum declaration
    result += "} pastix_" + name + "_t;\n\n"
    return result

def gen_coeftype( enum ) :
    result = gen_enum_documentation( enum )
    result += "#define pastix_coeftype_t spm_coeftype_t\n"

    maxsize = max( list( map( lambda x : len(x["NAME"]) + 1, enum["VALUES"].values() ) ) )

    for value in enum["VALUES"].values() :
        result += "#define "+ value["NAME"] + ( maxsize - len(value["NAME"]) ) * " " + value["VALUE"] +"\n"
    result += close_bracket + "\n"
    return result

def genPastixEnums( enums ) :
    """
    Generate a string that corresponds to the content of pastix_enums.h, right
    after the declaration of pastix_iparm_t and pastix_dparm_t.

    @in  enums : The enums dictionnary.
    @out The string that declare the enums in api.h
    """
    headerFile  = ""

    for name in enums :
        if name == "coeftype" :
            headerFile += gen_coeftype( enums[name] )
            continue

        # Begin BLAS-like enums
        if name == "layout" :
            headerFile += blas_str + "\n"

        headerFile += gen_enum_declaration( name, enums[name] )

    # Close BLAS-like enums bracket
    headerFile += close_bracket
    return headerFile
