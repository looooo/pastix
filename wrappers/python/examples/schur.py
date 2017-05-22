#!/usr/bin/env python3
"""
 @file schur.py

 PaStiX schur python example

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
import scipy.linalg as la
import numpy as np

# Set matrix A
n = 9
A = sps.spdiags([np.ones(n)*i for i in [4, -1, -1, -1, -1]],
                [0, 1, 3, -1, -3], n, n)

x0 = np.arange(n).reshape(n,1)
# Construct b as b = A * x_0
b = A.dot(x0)
x = b.copy()

tmp = np.eye(2).dot(np.ones(2))  # Hack to make sure that the mkl is loaded

# Convert the scipy sparse matrix to spm storage format
spmA = pastix.spm( A )
spmA.printInfo()

# Initialize parameters to default values
iparm, dparm = pastix.initParam()

# Startup PaStiX
pastix_data = pastix.init( iparm, dparm )

factotype = pastix.factotype.LLT
iparm[pastix.iparm.factorization] = factotype

# Initialize the Schur list
nschur = 2
schurlist = np.array( [2, 3] )
pastix.setSchurUnknownList ( pastix_data, schurlist )

# Perform analyze
pastix.task_analyze( pastix_data, spmA )

# Perform numerical factorization
pastix.task_numfact( pastix_data, spmA )

# Get the Schur complement
S = np.array( np.zeros( (nschur, nschur) ), order='F', dtype=A.dtype )
pastix.getSchur( pastix_data, S )

# 1- Apply P to b
pastix.subtask_applyorder( pastix_data, pastix.direction.Forward, x )

# 2- Forward solve on the non Schur complement part of the system
if factotype == pastix.factotype.LU:
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Lower, pastix.trans.NoTrans, pastix.diag.Unit, x )
else:
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Lower, pastix.trans.NoTrans, pastix.diag.NonUnit, x )

# 3- Solve the Schur complement part
x[-nschur:] = la.solve(S, x[-nschur:])

# 4- Backward solve on the non Schur complement part of the system
if factotype == pastix.factotype.LU:
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Upper, pastix.trans.NoTrans, pastix.diag.NonUnit, x )
else:
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Lower, pastix.trans.ConjTrans, pastix.diag.NonUnit, x )

# 5- Apply P^t to x
pastix.subtask_applyorder( pastix_data, pastix.direction.Backward, x )

# 6- Refine the solution
pastix.task_refine(pastix_data, b, x)

# Check solution
spmA.checkAxb( x0, b, x )

pastix.finalize( pastix_data, iparm, dparm )
