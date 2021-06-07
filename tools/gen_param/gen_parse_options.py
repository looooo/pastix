import os
import sys
import re
import string
import yaml
import time

const_str = ''' * This file is generated automatically. If you want to modify it, modify
 * ${PASTIX_HOME}/docs/pastix_params.yaml and run ${PASTIX_HOME}/tools/gen_param/gen_parm_files.py.
 *
 * @copyright 2004-'''+ time.strftime( "%Y" ) +''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @date '''+ time.strftime( "%Y-%m-%d" )

parse_options_header_begin = '''/**
 *
 * @file parse_options.h
 *
'''+ const_str +'''
 *
 */
#ifndef _parse_options_h_
#define _parse_options_h_

BEGIN_C_DECLS

pastix_iparm_t parse_iparm( const char *iparm  );
pastix_dparm_t parse_dparm( const char *dparm  );
int            parse_enums( const char *string );

'''

parse_options_header_end = '''
END_C_DECLS

#endif /* _parse_options_h_ */
'''

parse_options_begin = '''/**
 *
 * @file parse_options.c
 *
'''+ const_str +'''
 *
 */
#include "common.h"
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>

'''

parse_iparm_doc ='''/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse iparm keywords to return its associated index in the iparm array
 *
 * This function converts the string only for input parameters, output
 * parameters are not handled.
 *
 *******************************************************************************
 *
 * @param[in] iparm
 *          The iparm string to convert to enum.
 *
 *******************************************************************************
 *
 * @retval The index of the iparm in the array.
 * @retval -1 if the string is not an iparm parameter.
 *
 *******************************************************************************/
'''

parse_dparm_doc = '''/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse dparm keywords to return its associated index in the dparm array
 *
 * This function converts the string only for input parameters, output
 * parameters are not handled.
 *
 *******************************************************************************
 *
 * @param[in] dparm
 *          The dparm string to convert to enum.
 *
 *******************************************************************************
 *
 * @retval The index of the dparm in the array.
 * @retval -1 if the string is not a dparm parameter.
 *
 *******************************************************************************/
'''

parse_enums_doc = '''/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse enum values for iparm array, and return the enum value
 * associated to it.
 *
 *******************************************************************************
 *
 * @param[in] string
 *          The enum name to convert to its value
 *
 *******************************************************************************
 *
 * @retval The value if the enum associated to the string
 * @retval -1 if the string is not an enum in the pastix API.
 *
 *******************************************************************************/
'''

def gen_parse_param( params, parm ) :
    """
    Generate the parse_iparm/parse_dparm routines

    @in  params : The param dictionnary structure.
    @in  parm   : iparm/dparm
    @out The parse_iparm/parse_dparm routines string
    """
    result  = "pastix_"+ parm +"_t\n"\
              "parse_"+ parm +"( const char *"+ parm +" )\n"\
              "{\n"

    maxsize  = max( list( map(lambda x : len(x), params.keys() ) ) )
    for name in params :
        param = params[name]
        if param['INOUT'] != 'IN' :
            continue

        spaces = " " * ( maxsize - len(name) )

        result += "    if(0 == strcasecmp(\"" + name.lower() +"\","+ spaces +" "+ parm +")) { return "+ name +"; }\n"

    result += "\n    return -1;\n"\
              "}\n\n"

    return result

def gen_parse_enums( iparms, enums ) :
    """
    Generate the parse_enums routine

    @in  iparms : The IPARM dictionnary.
    @in  enums  : The ENUMS dictionnary.
    @out The parse_enums routine string
    """
    previous = ""
    result  = "int\n"\
              "parse_enums( const char *string )\n"\
              "{\n"

    for name in iparms :
        iparm = iparms[name]
        if iparm['INOUT'] != 'IN' :
            continue

        # We have an enum IPARM
        if 'ENUM' not in iparm :
            continue
        enumname = iparm['ENUM']

        # Avoid IPARM that may point to the same enums
        if previous == enumname :
            continue
        previous = enumname

        values      = enums[enumname]['VALUES']
        currMaxSize = max( list( map(lambda x : len(x['NAME']), values.values()) ) )
        for field in values.values() :
            name = field["NAME"]
            spaces = " " * ( currMaxSize - len(name) )
            result += "    if(0 == strcasecmp(\"" + name.lower() +"\","+ spaces +" string)) { return "+ name +"; }\n"

        result +="\n"

    # If the value is directly given without string
    result += "    /* If the value is directly given without string */\n"\
              "    {\n"\
              "        int value;\n"\
              "        if ( sscanf( string, \"%d\", &value ) != 1 ) {\n"\
              "            return -1;\n"\
              "        }\n"\
              "        else {\n"\
              "            return value;\n"\
              "        }\n"\
              "    }\n"\
              "}\n\n"

    return result

def gen_enum_get_str_declaration( enum, inHeader ) :
    ftype = "const char*"
    name  = "pastix_"+ enum +"_getstr( pastix_"+  enum +"_t value )"

    if inHeader :
        result = ftype + " " + name + ";\n"
    else :
        result = ftype + "\n" + name + "\n"

    return result

def gen_enum_get_str( name, enum ) :
    """
    Generate a string that corresponds to the declaration of the enum_getstr routine

    @in  enum : The YAML structure that corresponds to the current enum.
    @out the enum declaration string for the pastix_[enum]_getstr routine
    """
    work = enum
    # Factotype has redondant values for its enum:
    # We have to work with this exception
    if name == "factotype" :
        work = { 'VALUES' : {} }
        for i in range(5, 10):
            work['VALUES'][i] = enum['VALUES'][i]

    result = gen_enum_get_str_declaration( name, False) + "{\n"\
             "    switch( value ) {\n"

    for field in work['VALUES'].values() :
        valname = field.get('NAME')
        result += "    case "+ valname + ":\n"
        result += "        return \""+ valname + "\";\n"

    result += "    default :\n        return \"Bad " + name + " given\";\n    }\n}\n\n"
    return result

def genParseOptC( iparms, dparms, enums ) :
    result  = parse_options_begin

    result += parse_iparm_doc
    result += gen_parse_param( iparms, "iparm" )

    result += parse_dparm_doc
    result += gen_parse_param( dparms, "dparm" )

    result += parse_enums_doc
    result += gen_parse_enums( iparms, enums )

    return result

def genParseOptGetstr( enums ) :
    cFile      = ""
    headerFile = parse_options_header_begin

    for name in enums :
        headerFile += gen_enum_get_str_declaration( name, True )
        cFile      += gen_enum_get_str( name, enums[name] )

    headerFile += parse_options_header_end
    return headerFile, cFile