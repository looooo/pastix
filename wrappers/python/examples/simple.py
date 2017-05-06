"""
 @file simple.py

 PaStiX simple python example

 @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @date 2017-05-04

"""
#from pypastix import pastix
from pastix_api import *
from pyspm import spm
from pypastix import pastix
import scipy.sparse as spp
import numpy as np

flttype = pastix_coeftype.PastixDouble
mtxtype = pastix_mtxtype.PastixGeneral

# Get corresponding numpy type for arithmetic and integers array
nptype = 'f8'
npinttype = 'i8'

iscomplex=0
if flttype == pastix_coeftype.PastixComplex32 or flttype == pastix_coeftype.PastixComplex64:
    iscomplex = 1

# Set matrix A
n    = 9
nrhs = 1
row  = np.array([0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8], dtype=npinttype)
col  = np.array([0, 1, 3, 0, 1, 2, 4, 1, 2, 5, 0, 3, 4, 6, 1, 3, 4, 5, 7, 2, 4, 5, 8, 3, 6, 7, 4, 6, 7, 8, 5, 7, 8], dtype=npinttype)
data = np.array([4.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, -1.0, 4.0,
                 -1.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, 4.0], dtype=nptype)

if iscomplex:
    data.imag += 1

# Construct the matrix in SciPY sparse matrix format
A = spp.coo_matrix((data, (row, col)), shape=(n, n))

# Construct an initial solution
x0 = np.zeros((n, nrhs,), dtype=nptype)
for i in range(nrhs):
    k = i * nrhs
    if iscomplex:
        x0[:,i].real = np.arange(k+1.0, k+n+1.0)
        x0[:,i].imag = np.arange(-n-k, -k)
    else:
        x0[:,i] = np.arange(k+1.0, k+n+1.0, dtype=nptype)

# Construct b as b = A * x_0
b = np.matmul(A.todense(), x0)
x = b.copy()

# Convert the scipy sparse matrix to spm storage format
spmA = spm();
spmA.fromspp( A, mtxtype );
spmA.printInfo()

# Initialize parameters to default values
iparm = np.array( np.zeros( pastix_iparm.iparm_size ), dtype=pastix_np_int )
dparm = np.array( np.zeros( pastix_dparm.dparm_size ), dtype='float64' )
pastix.initParam( iparm, dparm )

# Startup PaStiX
pastix_data = pastix.init( iparm, dparm )

# Perform analyze
pastix.analyze( pastix_data, spmA )

# Perform numerical factorization
pastix.numfact( pastix_data, spmA )

# Perform solve
pastix.solve( pastix_data, spmA, nrhs, x, n )

# Perform refinement
pastix.refine( pastix_data, nrhs, b, n, x, n )

# Check solution
spmA.checkAxb( x0, b, x )

pastix.finalize( pastix_data, iparm, dparm )
