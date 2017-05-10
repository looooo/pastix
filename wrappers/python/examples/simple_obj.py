"""
 @file simple_obj.py

 PaStiX simple python example with the Solver object

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

# Factorize
solver = pastix.Solver(A, verbose=True)

# Solve
x = solver.solve(b, x_check=x0)
solver.finalize()
