import sys
import os
from pastix_iparm import iparm
from pastix_dparm import dparm
from pastix_enums import enums
from gen_api_header import genApiHeaderFile
from gen_parse_options import genParseOptH, genParseOptC

def genPastixEnumsFiles( apiHeader, parseoptH, parseoptC ) :

    apiFile = open( apiHeader , "w" )
    apiFile.write( genApiHeaderFile( iparm, dparm, enums ) )
    apiFile.close()

    header, content = genParseOptH( enums )
    parseopt = open( parseoptC , "w" )
    parseopt.write( genParseOptC(iparm, dparm, enums) )
    parseopt.write( content )
    parseopt.close()

    parseopt = open( parseoptH , "w" )
    parseopt.write( header )
    parseopt.close()

##################################################################################
#                                   Main                                         #
##################################################################################

if len(sys.argv) < 2 :
    print("Usage : python3 gen_parm_files.py ${PASTIX_HOME}")
    exit(1)

pastixHome = sys.argv[1]
apiHeader  = os.path.join( pastixHome, "include", "pastix/api.h" )
parseoptC  = os.path.join( pastixHome, "common",  "parse_options.c" )
parseOptH  = os.path.join( pastixHome, "common",  "parse_options.h" )

genPastixEnumsFiles( apiHeader, parseOptH, parseoptC )
