"""
 @file schur-obj.py

 PaStiX Schur python example with an object oriented programing solution.

 @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-05-04

 This example shows how to use pastix solver to solve a system with a
 Schur complement.

"""
import pypastix as pastix
import scipy.sparse as sps
import scipy.linalg as la
import numpy as np

# Set matrix A
n = 9
A = sps.spdiags([np.ones(n)*i for i in [4, -1, -1, -1, -1]],
                [0, 1, 3, -1, -3], n, n)

x0 = np.arange(n).reshape(n,1)
# Construct b as b = A * x_0
b  = A.dot(x0)

tmp = np.eye(2).dot(np.ones(2))  # Hack to make sure that the mkl is loaded
solver = pastix.solver()

solver.schur(A, [2, 3])
S = solver.S
f = solver.b2f(b)
y = la.solve(S, f)
x = solver.y2x(y, b, x0=x0, check=True)
