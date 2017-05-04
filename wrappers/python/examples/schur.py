from mpi4py import MPI
from pymumps import Mumps
import scipy.sparse as spp
import scipy.linalg as sl
import numpy as np

# Initialize mumps wrapper
driver_mumps = Mumps('D')
# Get corresponding numpy type
nptype = driver_mumps.nptype

# MPI stuff
comm = MPI.COMM_WORLD
myid = comm.Get_rank()
nprocs = comm.Get_size()

# Set fortran communicator
id = driver_mumps.id
id.comm_fortran = comm.py2f()

# Set par to 1 for sequential version
if nprocs == 1:
    id.par = 1

# Init mumps
driver_mumps.initialize()

# Set matrix A
n = 3
nnz = 6
row = np.zeros((nnz,), dtype='i')
col = np.zeros((nnz,), dtype='i')
data = np.zeros((nnz,),dtype=nptype)
cnt_ = [0] # To use closure in the function

def set_ijv(i, j, v):
    cnt = cnt_[0]
    row[cnt] = i
    col[cnt] = j
    data[cnt] = v
    cnt_[0] = cnt_[0] + 1

set_ijv(0,0, 1.0)
set_ijv(0,1, 2.0)
set_ijv(0,2, 3.0)
set_ijv(1,0, 1.0)
set_ijv(2,0, -2.0)
set_ijv(1,2, 10.0)

A = spp.coo_matrix((data, (row, col)), shape=(n, n))

print (A.todense())

# Set rhs
#rhs = np.arange(n, dtype=nptype)

if myid == 0:
    driver_mumps.set_A(A)
    #driver_mumps.set_RHS(rhs)

# Set icntl
driver_mumps.ICNTL[1] = 6
driver_mumps.ICNTL[2] = 6
driver_mumps.ICNTL[3] = 6
driver_mumps.ICNTL[4] = 6

# Give indices to keep in the schur complement
driver_mumps.set_schur_listvar([1,2]) # Python numbering

# Call mumps
driver_mumps.drive(1) # Analysis
driver_mumps.drive(2) # Facto

if myid == 0:
    schur = driver_mumps.get_schur()
    print("Found:")
    print(schur)

# Finalize mumps
driver_mumps.finalize()

if myid == 0:
    print("done")
