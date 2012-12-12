#! /usr/bin/env python
# -*- coding: utf-8 -*-

# -----------------------------------------
# PLASMA installer
# University of Tennessee Knoxville
# July 4, 2009
# ----------------------------------------

import sys
import os

compilers = { "gcc": { "CC"        : "gcc",
                       "FC"        : "gfortran",
                       "LOADER"    : "$(FC)",
                       "ARCH"      : "ar",
                       "ARCHFLAGS" : "cr",
                       "RANLIB"    : "/usr/bin/ranlib",
                       "CFLAGS"    : "-O2 -DADD_",
                       "FFLAGS"    : "-O2 ",
                       "LDFLAGS"   : ""},
              "intel":{"CC"        : "icc",
                       "FC"        : "ifort",
                       "LOADER"    : "$(FC)",
                       "ARCH"      : "ar",
                       "ARCHFLAGS" : "cr",
                       "RANLIB"    : "/usr/bin/ranlib",
                       "CFLAGS"    : "-O2 -DADD_ -diag-disable vec",
                       "FFLAGS"    : "-O2 -fltconsistency -fp_port",
                       "LDFLAGS"   : "-O2 -nofor_main"},
              "xlc": { "CC"        : "xlc",
                       "FC"        : "xlf",
                       "LOADER"    : "$(FC)",
                       "ARCH"      : "ar",
                       "ARCHFLAGS" : "cr",
                       "RANLIB"    : "ranlib",
                       "CFLAGS"    : "-O2 -qstrict -qthreaded -DNOCHANGE",
                       "FFLAGS"    : "-O2 -qstrict -qthreaded",
                       "LDFLAGS"   : ""},
              "pgi": { "CC"        : "pgcc",
                       "FC"        : "pgfortran",
                       "LOADER"    : "$(FC)",
                       "ARCH"      : "ar",
                       "ARCHFLAGS" : "cr",
                       "RANLIB"    : "ranlib",
                       "CFLAGS"    : "-O2 -DADD_",
                       "FFLAGS"    : "-O2 ",
                       "LDFLAGS"   : "-O2 -Mnomain"}
              }

machines = { #"anaka"      : { "gcc" : { "goto"    : { "dir"    : "/home/fike/buildbot/plasma/slaves/anaka/anaka_gcc_goto",
             #                                         "BLAS"   : "-L/home/mfaverge/opt/ia64/gcc/goto-anaka -lgoto",
             #                                         "CBLAS"  : "",
             #                                         "LAPACK" : "-ltmg",
             #                                         "CLAPACK": "/home/mfaverge/opt/ia64/gcc" }, 
             #                           "goto2"   : { "dir"    : "/home/fike/buildbot/plasma/slaves/anaka/anaka_gcc_goto2",
             #                                         "BLAS"   : "-L/home/mfaverge/opt/ia64/gcc/goto2-anaka -lgoto2",
             #                                         "CBLAS"  : "",
             #                                         "LAPACK" : "-ltmg",
             #                                         "CLAPACK": "/home/mfaverge/opt/ia64/gcc" },
             #                           "refblas" : { "dir"    : "/home/fike/buildbot/plasma/slaves/anaka/anaka_gcc_refblas",
             #                                         "BLAS"   : "-lrefblas",
             #                                         "CBLAS"  : "-L/home/mfaverge/opt/ia64/gcc/lib -lcblas",
             #                                         "LAPACK" : "-ltmg -llapack",
             #                                         "CLAPACK": "/home/mfaverge/opt/ia64/gcc" }
             #                           },
             #                 },
             "battlecat0" : { "gcc" : { "atlas"   : { "BLAS"   : "/iclscratch3/homes/fike/x86/atlas_gcc/lib/libf77blas.a /iclscratch3/homes/fike/x86/atlas_gcc/lib/libatlas.a",
                                                      "CBLAS"  : "-lcblas ",
                                                      "LAPACK" : "-ltmg /iclscratch3/homes/fike/x86/atlas_gcc/lib/liblapack.a",
                                                      "CLAPACK": "/home/mfaverge/opt/x86/gcc" },
                                        #"atlas"   : { "BLAS"   : "-L/iclscratch3/homes/fike/x86/atlas_gcc/lib -lf77blas -latlas",
                                        #              "CBLAS"  : "-lcblas ",
                                        #              "LAPACK" : "-ltmg -L/iclscratch3/homes/fike/x86/atlas_gcc/lib -llapack",
                                        #              "CLAPACK": "/home/mfaverge/opt/x86/gcc" },
                                        "goto2"   : { "BLAS"   : "-L/iclscratch3/homes/fike/x86/goto2_gcc -lgoto2",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "-ltmg ",
                                                      "CLAPACK": "/home/mfaverge/opt/x86/gcc" },
                                        "refblas" : { "BLAS"   : "-lrefblas",
                                                      "CBLAS"  : "-L/home/mfaverge/opt/battlecat0/lib -lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/battlecat0" }
                                        },
                              },
             "beaker"     : { "gcc" : { "refblas" : { "BLAS"   : "-L/silk/homes/buildbot/plasma/slaves/beaker -lrefblas",
                                                      "CBLAS"  : "",   #Cblas missing
                                                      "LAPACK" : "",   #Tmg Missing
                                                      "CLAPACK": ""  }
                                        },
                              },
             "bluegrass"  : { "gcc" : { "atlas"   : { "dir"    : "/home/fike/buildbot/plasma/slaves/bluegrass/bluegrass_gcc_atlas",
                                                      "BLAS"   : "-L/home/fike/ATLAS/my_build_dir/lib -lf77blas -latlas",
                                                      "CBLAS"  : "-lcblas ",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/ppc/gcc" }, 
                                        "goto2"   : { "dir"    : "/home/fike/buildbot/plasma/slaves/bluegrass/bluegrass_gcc_goto2",
                                                      "BLAS"   : "-L/home/fike/lib/gotoblas2_gcc -lgoto2",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/ppc/gcc" },
                                        "refblas" : { "dir"    : "/home/fike/buildbot/plasma/slaves/bluegrass/bluegrass_gcc_refblas",
                                                      "BLAS"   : "-L/home/faverge/opt/ppc/gcc/lib -lrefblas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/ppc/gcc" }
                                        },
                              },
             "brutus"     : { "gcc" : { "acml"    : { "BLAS"   : "/home/mfaverge/opt/x86_64/gcc/acml/gfortran64/lib/libacml.a",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" }, 
                                        "goto2"   : { "BLAS"   : "-L/iclscratch3/homes/fike/x86_64/gotoblas2_gcc -lgoto2",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" },
                                        "refblas" : { "BLAS"   : "-L/home/mfaverge/opt/x86_64/gcc/lib -lrefblas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" }, 
                                        },
                              },
             "ig"         : { "gcc" : { "refblas" : { "BLAS"   : "-L/iclscratch3/homes/fike/x86_64/BLAS -lrefblas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" },
                                        },
                              "pgi" : { "refblas" : { "BLAS"   : "-L/iclscratch3/homes/fike/x86_64/refblas_pgi -lrefblas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/pgi" }
                                        },
                              },
             "jump"       : { "gcc" : { "goto"    : { "dir"    : "/home5/jzam11/jzam1115/public/jump_plasma/jump-plasma-xlc",
                                                      "BLAS"   : "-L/usr/local/lib/ -lblas_Goto.a",
                                                      "CBLAS"  : "",   #Cblas missing
                                                      "LAPACK" : "",   #Lapack Missing / Tmg Missing
                                                      "CLAPACK": "" } 
                                        },
                              },
             "torc10"     : { "gcc" : { "acml"    : { "BLAS"   : "/home/mfaverge/opt/x86_64/gcc/acml/gfortran64/lib/libacml.a",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" }, 
                                        "refblas" : { "BLAS"   : "-L/silk/homes/buildbot/plasma/slaves/torc10 -lrefblas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" },
                                        },
                              },
             "zoot"       : { "gcc" : { "atlas"   : { "BLAS"   : "-L/mnt/scratch/yarkhan/atlas-gcc -lf77blas -latlas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" },
#                                         "goto"    : { "BLAS"   : "-L/mnt/scratch/sw/gotoblas-gcc -lgoto",
#                                                       "CBLAS"  : "",  #Cblas Missing
#                                                       "LAPACK" : "",  #Lapack Missing / Tmg Missing
#                                                       "CLAPACK": ""/home/mfaverge/opt/x86_64/gcc },
                                        "goto2"   : { "BLAS"   : "-L/home/mfaverge/opt/x86_64/gcc/goto-zoot -lgoto2",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" },
                                        "mkl"     : { "BLAS"   : "-L/mnt/scratch/sw/intel/2011.6.233/mkl/lib/em64t -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "",   
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc"  },
                                        "refblas" : { "BLAS"   : "-L/home/mfaverge/opt/x86_64/gcc/lib -lrefblas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/gcc" } 
                                        },
                              "intel":{ "acml"    : { "BLAS"   : "/home/mfaverge/opt/x86_64/intel/acml/ifort64/lib/libacml.a",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/intel" }, 
                                        "atlas"   : { "BLAS"   : "-L/iclscratch3/homes/fike/x86_64/atlas_icc/lib -lf77blas -latlas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/intel" },
                                        "goto"    : { "BLAS"   : "-L/mnt/scratch/sw/gotoblas-intel -lgoto",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/intel" }, 
                                        "goto2"   : { "BLAS"   : "-L/home/mfaverge/opt/x86_64/intel/goto-zoot -lgoto2",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "-ltmg",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/intel" },
                                        "mkl"     : { "BLAS"   : "-L/mnt/scratch/sw/intel/2011.6.233/mkl/lib/em64t -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group",
                                                      "CBLAS"  : "",
                                                      "LAPACK" : "",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/intel"  },
                                        "refblas" : { "BLAS"   : "-L/home/mfaverge/opt/x86_64/intel -lrefblas",
                                                      "CBLAS"  : "-lcblas",
                                                      "LAPACK" : "-ltmg -llapack",
                                                      "CLAPACK": "/home/mfaverge/opt/x86_64/intel" }
                                        },
                              
                              },
             }


def write_makeinc(arch, cc, lib, branch):
    """ Writes the make.inc file """
    
    compiler = compilers[cc]
    conf     = machines[arch][cc][lib]
 
    name = 'make.inc.'+cc+'.'+lib
    dir  = conf.get("dir", '/silk/homes/buildbot/plasma/slaves/'+arch+'/'+arch+'_'+cc+'_'+lib)
    if not (branch == ""):
        name = name+'.'+branch
        dir  = dir+'_'+branch

    name = name+'.'+arch
    dir  = dir+'/build'

    print 'Writing ',name,'...'
    sys.stdout.flush()
    makeinc ="""
#/////////////////// P /// L /// A /// S /// M /// A //////////////////
#/// PLASMA is a software package provided by Univ. of Tennessee,  ///
#/// Univ. of California Berkeley and Univ. of Colorado Denver     ///
#//////////// M /// A /// K /// E /// . /// I /// N /// C /////////////

#///////////// U /// S /// E /// R ////// P /// A /// R /// T //////////

MAKE = make -j 4
PLASMA_DIR  = """+dir+"""

CC          = """+compiler["CC"]+"""
FC          = """+compiler["FC"]+"""
LOADER      = """+compiler["LOADER"]+"""

ARCH        = """+compiler["ARCH"]+"""
ARCHFLAGS   = """+compiler["ARCHFLAGS"]+"""
RANLIB      = """+compiler["RANLIB"]+"""

CFLAGS      = """+compiler["CFLAGS"]+"""
FFLAGS      = """+compiler["FFLAGS"]+"""
LDFLAGS     = """+compiler["LDFLAGS"]+"""

# Blas Library
LIBBLAS     = """+conf["BLAS"]+"""
# CBlas library
LIBCBLAS    = """+conf["CBLAS"]+"""
# lapack and tmg library (lapack is included in acml)
LIBLAPACK   = """+conf["LAPACK"]+"""
INCCLAPACK  = -I"""+conf["CLAPACK"]+"""/include
LIBCLAPACK  = -L"""+conf["CLAPACK"]+"""/lib -llapacke

"""
    writefile(name, makeinc)

    print 'done.'


def writefile(fname, fill):
    """ writes the file fname with content fill """
    fp = open(fname,'w')
    fp.write(fill)
    fp.close()
    

def main(argv):
    
    for arch in machines :
        for cc in machines[arch]:
            for lib in machines[arch][cc] :
                write_makeinc(arch, cc, lib, "")
                branches = machines[arch][cc][lib].get("branches", ())
                for i in branches :
                    print i
                    write_makeinc(arch, cc, lib, i)

    return 0

if "__main__" == __name__:
    sys.exit(main(sys.argv))
