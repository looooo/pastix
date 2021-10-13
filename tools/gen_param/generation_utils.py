"""
 @file generation_utils.py

 Utils for the automatically generated files.

 @copyright 2021-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.2.1
 @author Tony Delarue
 @date 2021-04-07

"""
import time

const_str = ''' * This file is generated automatically. If you want to modify it, modify
 * ${PASTIX_HOME}/tools/gen_param/pastix_[iparm/dparm/enums].py and run
 * ${PASTIX_HOME}/tools/gen_param/gen_parm_files.py ${PASTIX_HOME}.
 *
 * @copyright 2004-'''+ time.strftime( "%Y" ) +''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @date '''+ time.strftime( "%Y-%m-%d" )

def headerFileBegin( filename ):
    return '''/**
 *
 * @file '''+ filename +'''.h
 *
'''+ const_str +'''
 *
 */
#ifndef _'''+filename+'''_h_
#define _'''+filename+'''_h_

BEGIN_C_DECLS

'''

def headerFileEnd( filename ):
    return '''
END_C_DECLS

#endif /* _'''+filename+'''_h_ */
'''

def contentFileBegin( filename ) :
    return '''/**
 *
 * @file '''+ filename +'''.c
 *
'''+ const_str +'''
 *
 */
#include "common.h"
'''

def findEnumFromName( name, enums ):
    for enum in enums :
        if name == enum["name"] :
            return enum
    return -1

prototype   = '''{ftype} {fname}( {atype} {aname} );
'''
declaration = '''{ftype}
{fname}( {atype} {aname} )
{bracket}
'''

enumprot =   prototype.format( ftype='{ftype}', fname='pastix_{name}_{routine}', atype='pastix_{name}_t', aname='value', bracket='{bracket}' )
enumdecl = declaration.format( ftype='{ftype}', fname='pastix_{name}_{routine}', atype='pastix_{name}_t', aname='value', bracket='{bracket}' )
