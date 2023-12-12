"""
 @file gen_check_options.py

 Declaration of the check_options.c file that is generated automatically.

 @copyright 2021-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.2
 @author Tony Delarue
 @date 2023-07-21

"""
import generation_utils as gu

check_value_header_begin = gu.headerFileBegin( "check_options" )

check_value_header_end = gu.headerFileEnd( "check_options" )

check_value_content_begin = gu.contentFileBegin( "check_options" ) + "\n"

idparm_doc = '''/**
 *******************************************************************************
 *
 * @brief Check the values of the {idparm}.
 *
 *******************************************************************************
 *
 * @param[in] {idparm}
 *          The {idparm} options array.
 *
 *******************************************************************************
 *
 * @return The amount of incorrect values in the {idparm} array.
 *
 *******************************************************************************/
'''

def gen_enum_check_value( enum ) :
    """
    Generate a string that corresponds to the declaration of the enum_check_value routine

    @in  enum : The current enum dictionnary.
    @out the enum declaration string for the pastix_[enum]_check_value routine.
    """
    work = enum
    name = enum['name']
    # Factotype has redondant values for its enum:
    # We have to work with this exception
    if name == "factotype" :
        work = { 'values' : [] }
        for i in range(5, 10):
            work['values'].append( enum['values'][i] )

    result  = gu.enumdecl.format( ftype="static inline int", name=name, routine="check_value", bracket="{" )
    result += "    if("
    spaces  = " "
    indent  = " " * len( "    if( " )

    for value in work['values'] :
        valname = value['name']
        result += spaces + "(value == "+ valname + ") ||\n"
        spaces  = indent

    result  = result[:len(result)-4]
    result += ''' ) {
        return 0;
    }
    return 1;
}

'''
    return result

def iparmRangeCheckValue( iparm ) :
    result  = gu.declaration.format( ftype="static inline int", fname=iparm["name"]+'_check_value',
                                     atype='pastix_int_t', aname='iparm', bracket="{" )
    result += '    /* TODO : Check range iparm['+ iparm['name'].upper() +'] */\n'
    result += '''    (void)iparm;
    return 0;
}

'''
    return result

def iparmEnumCheckValue( iparm ) :
    result  = gu.declaration.format( ftype="static inline int", fname=iparm["name"]+'_check_value',
                                     atype='pastix_int_t', aname='iparm', bracket="{" )
    result += '    int rc;\n'
    result += '    rc = pastix_'+ iparm['enum'] +'_check_value( iparm );\n'
    result += r'''    if ( rc == 1 ) {
        fprintf(stderr, "'''+ iparm["name"].upper() +r''': The value is incorrect\n");
    }
    return rc;
}

'''
    return result

def dparmRangeCheckValue( dparm ) :
    result  = gu.declaration.format( ftype="static inline int", fname=dparm["name"]+'_check_value',
                                     atype='double', aname='dparm', bracket="{" )
    result += '    /* TODO : Check range dparm['+ dparm['name'].upper() +'] */\n'
    result += '''    (void)dparm;
    return 0;
}

'''
    return result

def genCheckOpt( iparms, dparms, enums ) :
    """
    Write pastix_options.h completely.
    Declare all pastix_[enum]_check_value in the content and header file.

    @in  enums  : The array containing enums.
    @out headerFile : The complete parse_options.h string.
    @out cFile      : The string for all the pastix_[enum]_getstr declarations.
    """
    cFile      = check_value_content_begin
    previous = ""

    enums_str = ""
    iparm_str = ""
    iparm_all = gu.declaration.format( ftype="int", fname='iparm_check_values',
                                       atype='const pastix_int_t', aname='*iparm', bracket="{" )
    iparm_all += "    int error = 0;\n"

    for groups in iparms :
        for iparm in groups["subgroup"] :
            if iparm['access'] != 'IN' :
                continue

            # We need an enum iparm
            if 'enum' in iparm :
                enumname = iparm['enum']

                # Avoid IPARM that may point to the same enums
                if previous == enumname :
                    continue
                previous = enumname

                enum = gu.findEnumFromName( enumname, enums )
                enums_str += gen_enum_check_value( enum )
                iparm_str += iparmEnumCheckValue( iparm )
            else :
                iparm_str += iparmRangeCheckValue( iparm )

            iparm_all += '    error += ' + iparm["name"]+'_check_value( iparm['+ iparm['name'].upper() +'] );\n'

    iparm_all += "    return error;\n}\n\n"

    dparm_all  = gu.declaration.format( ftype="int", fname='dparm_check_values',
                                        atype='const double', aname='*dparm', bracket="{" )
    dparm_all += "    int error = 0;\n"
    dparm_str  = ""
    for dparm in dparms :
        if dparm['access'] != 'IN' :
            continue
        dparm_str += dparmRangeCheckValue( dparm )
        dparm_all += '    error += ' + dparm["name"]+'_check_value( dparm['+ dparm['name'].upper() +'] );\n'

    dparm_all += "    return error;\n}\n"

    cFile += enums_str + iparm_str + dparm_str
    cFile += idparm_doc.format(idparm="iparm") + iparm_all
    cFile += idparm_doc.format(idparm="dparm") + dparm_all

    return cFile
