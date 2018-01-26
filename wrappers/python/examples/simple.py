#!/usr/bin/env python3
"""
 @file simple.py

 PaStiX simple python example

 @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-05-04

"""
import pypastix as pastix
import scipy.sparse as sps
import numpy as np

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Load a sparse matrix from RSA driver
spmA = pastix.spm( None, driver=pastix.driver.Laplacian, filename="10:10:10:2.:1." )
#spmA = pastix.spm( None, driver=driver.RSA, filename="$PASTIX_DIR/test/matrix/oilpan.rsa" )
spmA.printInfo()

# Scale A for low-rank: A / ||A||_f
norm = spmA.norm()
spmA.scale( 1. / norm )

# Generate b and x0 vector such that A * x0 = b
x0, b = spmA.genRHS( pastix.rhstype.RndX )

# Initialize parameters to default values
iparm, dparm = pastix.initParam()

# Startup PaStiX
pastix_data = pastix.init( iparm, dparm )

# Change some parameters
iparm[pastix.iparm.factorization] = pastix.factotype.LU

# Perform analyze
pastix.task_analyze( pastix_data, spmA )

# Perform numerical factorization
pastix.task_numfact( pastix_data, spmA )

# Perform solve
x = b.copy()
pastix.task_solve( pastix_data, x )

# Refine the solution
pastix.task_refine( pastix_data, b, x )

# Check solution
spmA.checkAxb( x0, b, x )

pastix.finalize( pastix_data, iparm, dparm )
