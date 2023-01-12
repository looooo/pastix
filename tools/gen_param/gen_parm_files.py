#!/usr/bin/env python
"""
 @file gen_parm_files.py

 Python script to generate the parameters files.

 @copyright 2021-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.0
 @author Tony Delarue
 @date 2021-10-19

"""
import sys
import os
from pastix_iparm import iparm
from pastix_dparm import dparm
from pastix_enums import enums
from gen_api_header import genApiHeaderFile
from gen_parse_options import genParseOptH, genParseOptC
from gen_check_options import genCheckOpt
from gen_parm_2_csv import genParm2csv
from gen_pastix_completion import genCompletion

##################################################################################
#                                   Main                                         #
##################################################################################

if len(sys.argv) < 2 :
    print("Usage : python3 gen_parm_files.py ${PASTIX_HOME}")
    exit(1)

pastixHome = sys.argv[1]
apiHeader  = os.path.join( pastixHome, "include", "pastix/api.h" )
parseOptC  = os.path.join( pastixHome, "common",  "parse_options.c" )
parseOptH  = os.path.join( pastixHome, "common",  "parse_options.h" )
checkOptC  = os.path.join( pastixHome, "common",  "check_options.c" )
checkOptH  = os.path.join( pastixHome, "common",  "check_options.h" )
completion = os.path.join( pastixHome, "tools",   "pastix_completion.sh.in" )

apiFile = open( apiHeader , "w" )
apiFile.write( genApiHeaderFile( iparm, dparm, enums ) )
apiFile.close()

header, content = genParseOptH( enums )
parseopt = open( parseOptC , "w" )
parseopt.write( genParseOptC(iparm, dparm, enums) )
parseopt.write( content )
parseopt.write( genParm2csv(iparm, dparm) )
parseopt.close()

parseopt = open( parseOptH , "w" )
parseopt.write( header )
parseopt.close()

checkopt = open( checkOptC , "w" )
checkopt.write( genCheckOpt( iparm, dparm, enums ) )
checkopt.close()

compl = open( completion, "w" )
compl.write( genCompletion( iparm, dparm, enums ) )
compl.close()
