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
from pypastix import pastix
from pyspm import spm
import scipy.sparse as spp
import numpy as np

flttype = pastix.coeftype.PastixDouble
mtxtype = pastix.mtxtype.PastixGeneral

# Get corresponding numpy type
nptype = pastix.get_numpy_type( flttype )

iscomplex=0
if flttype == pastix.coeftype.PastixComplex32 or flttype == pastix.coeftype.PastixComplex64:
    iscomplex = 1

# Set matrix A
n    = 9
nrhs = 1
row  = np.array([0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8])
col  = np.array([0, 1, 3, 0, 1, 2, 4, 1, 2, 5, 0, 3, 4, 6, 1, 3, 4, 5, 7, 2, 4, 5, 8, 3, 6, 7, 4, 6, 7, 8, 5, 7, 8])
data = np.array([4.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, -1.0, 4.0,
                 -1.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, 4.0, -1.0, -1.0, -1.0, 4.0], dtype=nptype)

if iscomplex:
    data.imag += 1

# Construct the matrix in SciPY sparse matrix format
A = spp.coo_matrix((data, (row, col)), shape=(n, n))

# Construct an initial solution
x0 = np.zeros((n, NRHS,), dtype=nptype)
for i in range(NRHS):
    k = i * NRHS
    if iscomplex:
        x0[:,i].real = np.arange(k+1.0, k+n+1.0)
        x0[:,i].imag = np.arange(-n-k, -k)
    else:
        x0[:,i] = np.arange(k+1.0, k+n+1.0, dtype=nptype)

# Construct b as b = A * x_0
b = np.matmul(A.todense(), x0)
x = b

# Initialize parameters to default values
pastix.initParam( iparm, dparm )

# Convert the scipy sparse matrix to spm storage format
spmA = spm.convertSpp( flttype, mtxtype, A );

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
spm.checkAxb( nrhs, spmA, x0, n, b, n, x, n )

#spm.exit( spmA )
pastix.finalize( pastix_data )
