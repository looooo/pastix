"""
 @file gen_parse_options.py

 Declaration of the parse_options.[hc] files that are generated automatically.

 @copyright 2021-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.2.1
 @author Tony Delarue
 @date 2021-04-07

"""
import generation_utils as gu

parse_options_header_begin = gu.headerFileBegin( 'parse_options' ) + '''pastix_iparm_t parse_iparm( const char *iparm  );
pastix_dparm_t parse_dparm( const char *dparm  );
int            parse_enums( const char *string );

'''

parse_options_header_end = "\n" + gu.headerFileEnd( "parse_options" )

parse_options_begin = gu.contentFileBegin( "parse_options" ) +'''#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>

'''

idparm_doc = '''/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse {idparm} keywords to return its associated index in the {idparm} array
 *
 * This function converts the string only for input parameters, output
 * parameters are not handled.
 *
 *******************************************************************************
 *
 * @param[in] {idparm}
 *          The {idparm} string to convert to enum.
 *
 *******************************************************************************
 *
 * @retval The index of the {idparm} in the array.
 * @retval -1 if the string is not {article} {idparm} parameter.
 *
 *******************************************************************************/
'''
parse_iparm_doc = idparm_doc.format( idparm="iparm", article="an" )
parse_dparm_doc = idparm_doc.format( idparm="dparm", article="a" )

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

parse_idparm = gu.declaration.format( ftype='pastix_{idparm}_t', fname='parse_{idparm}', atype='const char', aname='*{idparm}', bracket='{bracket}' )

def gen_parse_iparm( iparms ) :
    """
    Generate the parse_iparm/parse_dparm routines

    @in  dparms : The param dictionnary structure.
    @in  parm   : iparm/dparm
    @out The parse_iparm/parse_dparm routines string
    """
    maxName   = 0
    for group in iparms :
        tmpname = max( list( map(lambda x : len(x['name']), group["subgroup"] ) ) )
        if tmpname > maxName :
            maxName = tmpname

    result  = parse_idparm.format( idparm="iparm", bracket='{' )
    newline = False
    for group in iparms :
        for iparm in group["subgroup"] :
            name = iparm["name"]
            if iparm['access'] != 'IN' :
                continue
            newline = True
            spaces = " " * ( maxName - len(name) )

            result += "    if(0 == strcasecmp(\"" + name +"\","+ spaces +" iparm)) { return "+ name.upper() +"; }\n"

        if newline :
            result +="\n"
            newline = False

    result += "    return -1;\n"\
              "}\n\n"

    return result

def gen_parse_dparm( dparms ) :
    """
    Generate the parse_iparm/parse_dparm routines

    @in  dparms : The param dictionnary structure.
    @in  parm   : iparm/dparm
    @out The parse_iparm/parse_dparm routines string
    """
    result = parse_idparm.format( idparm="dparm", bracket='{' )

    maxsize  = max( list( map(lambda x : len(x['name']), dparms ) ) )
    for dparm in dparms :
        name = dparm["name"]
        if dparm['access'] != 'IN' :
            continue

        spaces = " " * ( maxsize - len(name) )

        result += "    if(0 == strcasecmp(\"" + name +"\","+ spaces +" dparm)) { return "+ name.upper() +"; }\n"

    result += "\n    return -1;\n"\
              "}\n\n"

    return result

def find_enums_name( name, enums ):
    for enum in enums :
        if name == enum["name"] :
            return enum
    return -1

def gen_parse_enums( iparms, enums ) :
    """
    Generate the parse_enums routine

    @in  iparms : The IPARM dictionnary.
    @in  enums  : The ENUMS dictionnary.
    @out The parse_enums routine string
    """
    previous = ""
    result  = gu.declaration.format( ftype='int', fname='parse_enums', atype='const char', aname='*string', bracket='{' )

    for groups in iparms :
        for iparm in groups["subgroup"] :
            if iparm['access'] != 'IN' :
                continue

            # We have an enum IPARM
            if 'enum' not in iparm :
                continue
            enumname = iparm['enum']

            # Avoid IPARM that may point to the same enums
            if previous == enumname :
                continue
            previous = enumname

            values =  find_enums_name( enumname, enums )['values']
            currMaxSize = max( list( map(lambda x : len(x['name']), values) ) )
            for value in values :
                name = value["name"]
                spaces = " " * ( currMaxSize - len(name) )
                result += "    if(0 == strcasecmp(\"" + name.lower() +"\","+ spaces +" string)) { return "+ name +"; }\n"

            result +="\n"

    # If the value is directly given without string
    result += '''    /* If the value is directly given without string */
    {
        int value;
        if ( sscanf( string, "%d", &value ) != 1 ) {
            return -1;
        }
        else {
            return value;
        }
    }
}

'''
    return result

def gen_enum_get_str( enum ) :
    """
    Generate a string that corresponds to the declaration of the enum_getstr routine

    @in  enum : The current enum dictionnary.
    @out the enum declaration string for the pastix_[enum]_getstr routine.
    """
    work = enum
    name = enum['name']
    # Factotype has redondant values for its enum:
    # We have to work with this exception
    if name == "factotype" :
        work = { 'values' : [] }
        for i in range(5, 10):
            work['values'].append( enum['values'][i] )

    result  = gu.enumdecl.format( ftype="const char*", name=name, routine="getstr", bracket="{" )
    result += "    switch( value ) {\n"

    for value in work['values'] :
        valname = value['name']
        result += "    case "+ valname + ":\n"
        result += "        return \""+ valname + "\";\n"

    result += "    default :\n        return \"Bad " + name + " given\";\n    }\n}\n\n"
    return result

def genParseOptC( iparms, dparms, enums ) :
    """
    Generate a string that corresponds to parse_options.c

    @in  iparms : The array containing iparm groups.
    @in  dparms : The array containing dparms.
    @in  enums  : The array containing enums.
    @out the enum declaration string for the pastix_[enum]_getstr routine
    """
    result  = parse_options_begin

    result += parse_iparm_doc
    result += gen_parse_iparm( iparms )

    result += parse_dparm_doc
    result += gen_parse_dparm( dparms )

    result += parse_enums_doc
    result += gen_parse_enums( iparms, enums )

    return result

def genParseOptH( enums ) :
    """
    Write pastix_options.h completely.
    Declare all pastix_[enum]_getstr in the .c file.

    @in  enums  : The array containing enums.
    @out headerFile : The complete parse_options.h string.
    @out cFile      : The string for all the pastix_[enum]_getstr declarations.
    """
    cFile      = ""
    headerFile = parse_options_header_begin

    for enum in enums :
        headerFile += gu.enumprot.format( ftype="const char*", name=enum["name"], routine="getstr", bracket="{\n" )
        cFile      += gen_enum_get_str( enum )

    headerFile += "\nvoid pastix_param2csv( const pastix_data_t *pastix_data, FILE *csv );\n"
    headerFile += gu.prototype.format( ftype="int", fname='iparm_check_values',
                                       atype='const pastix_int_t', aname='*iparm' )
    headerFile += gu.prototype.format( ftype="int", fname='dparm_check_values',
                                       atype='const double      ', aname='*dparm' )

    headerFile  = headerFile[:len(headerFile) - 1] +  parse_options_header_end
    return headerFile, cFile
