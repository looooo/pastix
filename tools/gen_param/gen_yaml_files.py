import sys
import os
from pastix_iparm import iparm
from pastix_dparm import dparm
from pastix_enums import enums
from gen_pastix_api import *
from gen_parse_options import *

def genPastixEnumsFiles( apiHeader, parseoptH, pasrseoptC ) :

    apiFile = open( apiHeader , "w" )
    apiFile.write( pastix_enums_begin )
    apiFile.write( genIparmDeclaration(iparm) )
    apiFile.write( genDparmDeclaration(dparm) )
    apiFile.write( genPastixEnums(enums) )
    apiFile.write( pastix_enums_end )
    apiFile.close()

    header, content = genParseOptGetstr( enums )
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

pastixHome    = sys.argv[1]
apiHeader     = os.path.join( pastixHome, "include", "pastix/api.h" )
parseoptC     = os.path.join( pastixHome, "common",  "parse_options.c" )
parseOptH     = os.path.join( pastixHome, "common",  "parse_options.h" )

genPastixEnumsFiles( apiHeader, parseOptH, parseoptC )