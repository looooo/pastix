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

# Set the problem
n = 9
A = sps.spdiags([np.ones(n)*i for i in [4, -1, -1, -1, -1]],
                [0, 1, 3, -1, -3], n, n)
x0 = np.arange(n)
b = A.dot(x0)

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))



# Convert the scipy sparse matrix to spm storage format
spmA = pastix.spm( A )
spmA.printInfo()

# Initialize parameters to default values
iparm, dparm = pastix.initParam()

# Startup PaStiX
pastix_data = pastix.init( iparm, dparm )

factotype = pastix.factotype.LLT
iparm[pastix.iparm.factorization] = factotype

# Perform analyze
pastix.analyze( pastix_data, spmA )

# Perform numerical factorization
pastix.numfact( pastix_data, spmA )

# Perform solve
x = b.copy()
pastix.solve( pastix_data, spmA, x)

# Refine the solution
pastix.refine(pastix_data, b, x)

# Check solution
spmA.checkAxb( x0, b, x )

pastix.finalize( pastix_data, iparm, dparm )
