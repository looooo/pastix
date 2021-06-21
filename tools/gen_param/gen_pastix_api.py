import os
import sys
import re
import string
import time

const_str = ''' * This file is generated automatically. If you want to modify it, modify
 * ${PASTIX_HOME}/tools/gen_param/pastix_[iparm/dparm/enums].py and run
 * ${PASTIX_HOME}/tools/gen_param/gen_parm_files.py ${PASTIX_HOME}.
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
    maxName   = max( list( map(lambda x : len(x['name']),    params) ) )

    maxBrief  = max( list( map(lambda x : len(x['brief']),   params) ) )

    maxDefVal = max( list( map(lambda x : len(x['default']), params) ) )

    return maxName, maxBrief, maxDefVal

def gen_param( params ) :
    """
    Generate a string that corresponds to the declaration of the prefix

    @in  params : The params dictionnary
    @in  prefix : IPARM/DPARM
    @out the enum declaration string for the pastix_[i/d]parm_t
    """
    maxName, maxBrief, maxDefVal = get_param_max_str_size( params )
    result = ""
    for param in params :
        # Get the fields
        name    = param['name'].upper() + ","
        brief   = param['brief']
        default = param['default']
        inout   = param['access']

        # Align them correctly
        name    += " " * ( maxName   - len(name) + 1 )
        brief   += " " * ( maxBrief  - len(brief) + 3 )
        default += " " * ( maxDefVal - len(default) + 3  )
        inout   += " " * ( 3         - len(inout) )

        # Add them in the result
        result += "    " + name +" /**< "+ brief +" Default: "+ default + inout +" */\n"
    return result

#
# Wrappers for iparm/dparm
#
def genIparmDeclaration( iparms ) :
    result = "/**\n"\
             " * @brief Integer parameters\n"\
             " */\n"\
             "typedef enum pastix_iparm_e {\n"
    for group in iparms :
        if group["brief"] != "None" :
            result += "    /* "+ group["brief"] +" */\n"
        result += gen_param( group["subgroup"] )
        result += "\n"
    result += "    IPARM_SIZE\n"\
              "} pastix_iparm_t;\n\n"
    return result

def genDparmDeclaration( dparms ) :
    result = "/**\n"\
             " * @brief Float parameters\n"\
             " */\n"\
             "typedef enum pastix_dparm_e {\n"
    result += gen_param( dparms )
    result += "    DPARM_SIZE\n"\
              "} pastix_dparm_t;\n\n"
    return result

def gen_enum_documentation(enum):
    result = "/**\n"\
          " * @brief " + enum['doc']['brief'] + "\n"

    if 'details' in enum['doc'] :
        details = enum['doc']['details'].split("\n")
        details = details[:len(details) - 1]

        for line in details:
            if line == '' :
                result += " *\n"
                continue
            result += " * " + line + "\n"
    result += " */\n"

    return result

def gen_enum_declaration( enum ) :
    """
    Generate a string that corresponds to the declaration of an enum

    @in  enum : The YAML structure that corresponds to the enum to declare
    @out the enum declaration string for the pastix_[enum]_t
    """
    name = enum['name']
    # Write documentation
    result = gen_enum_documentation( enum )

    # End comment, begin enum declaration
    result += "typedef enum pastix_"+ name + "_e {\n"\

    values = enum['values']
    last   = list(enum['values'].keys())[-1]

    nameMaxSize  = max(list( map(lambda x : len(x['name']), values.values()) ))
    briefMaxSize = 0
    valueMaxSize = 0
    if 'brief' in values[0] :
        briefMaxSize = max(list( map(lambda x : len(x['brief']), values.values() ) ))
    if 'value' in values[0] :
        valueMaxSize = max(list( map(lambda x : len(x['value']), values.values() ) ))

    for value in values.values() :
        # Name
        result += "    " + value['name']

        # If a value is given, put it in the file
        if 'value' in value :
            result += " " * (nameMaxSize - len(value['name']) )
            result += " = "+ value['value']

        if value != enum['values'][ last ] :
            result += ","
        else :
            if 'brief' in value :
                result += " "

        # If a brief comment is given, put it in the file
        if 'brief' in value :
            if 'value' in value :
                result += " " * (valueMaxSize - len(value['value']) ) + " /**< "
            else :
                result += " " * (nameMaxSize  - len(value['name']) ) + " /**< "
            result += value['brief'] + " " * (briefMaxSize - len(value['brief']) ) + " */"
        result += "\n"

    # End of the enum declaration
    result += "} pastix_" + name + "_t;\n\n"
    return result

def gen_coeftype( enum ) :
    result = gen_enum_documentation( enum )
    result += "#define pastix_coeftype_t spm_coeftype_t\n"

    maxsize = max( list( map( lambda x : len(x["name"]) + 1, enum["values"].values() ) ) )

    for value in enum["values"].values() :
        result += "#define "+ value["name"] + ( maxsize - len(value["name"]) ) * " " + value["value"] +"\n"
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

    for enum in enums :
        name = enum["name"]
        if name == "coeftype" :
            headerFile += gen_coeftype( enum )
            continue

        # Begin BLAS-like enums
        if name == "layout" :
            headerFile += blas_str + "\n"

        headerFile += gen_enum_declaration( enum )

    # Close BLAS-like enums bracket
    headerFile += close_bracket
    return headerFile
