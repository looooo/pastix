# An example of PaStiX's Python interface: resolution of the Helmholtz equation

This example is about the resolution of the Helmholtz equation with [PaStiX](https://gitlab.inria.fr/solverstack/pastix) a sparse direct solver.
The Helmholtz equation:
$$ -\Delta u - k^{2}u=f \text{ on } \Omega $$

This example will show you how to create the matrix, solve the equation with PaStiX, and visualize the solution with [Paraview](https://www.paraview.org/).

## Configuration


Make sure that PaStiX and Paraview are installed on your computer. If this is not the case, the tutorial on how to build and install PaStiX is avalaible on the page [Build and Install PaStiX with CMake](https://solverstack.gitlabpages.inria.fr/pastix/md_docs_doxygen_chapters_Pastix_Install.html). The tutorial to install Paraview is on the [Paraview website](https://www.paraview.org/).

Once PaStiX is installed, you need to set the path to your PaStiX directory:

- On Mac :

```shell
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/pastix-directory/lib
export PYTHONPATH=/pastix-directory/lib/python
```

- On Linux:

```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pastix-directory/lib
export PYTHONPATH=/pastix-directory/lib/python
```

The following libraries are needed for the Helmholtz example:

```python
import spm
import pypastix as pastix
from scipy.sparse import csc_matrix
import numpy as np

from paraview.simple import *

import utilities
import random as random
```

To start with the python interface of PaStiX, you can use the example [simple.py](https://gitlab.inria.fr/solverstack/pastix/-/blob/master/wrappers/python/examples/simple.py). The *spm* and *pypastix* libraries are available on the [PaStiX gitlab](https://gitlab.inria.fr/solverstack/pastix/-/tree/master).

## Creation of the matrix: discretization of the equation

First the matrix's size *n* and number of non-zeros elements *nnz* have to be determined. The Helmholtz matrix is a 3D matrix with: *nrow* the number of rows, *ncol* the number of columns and *nprof* the number of tubes in the 3D matrix. These three dimensions are computed with the wave number *k*, the discretization *h* and the length *L* of the studied object.

```python
k = 15
kh = 0.625
h = kh / k
L = 1.0

nrow = int( L / h ) + 1
ncol = nrow
nprof = nrow

n = nrow * ncol * nprof
nnz = 7 * n - 6 * (ncol) ** 2
```

Then the matrix is created in one of the three storage types of the *spm* library: *COO*, *CSC* and *CSR*. The *CSC* format is the one used in this example but the *COO* and *CSR* formats can be used as well. The list *coltab* corresponds to the index of the start of each column in the *rowtab* and *valtab*. The *rowtab* corresponds to the list of the row indexes in the columns order and the *valtab* corresponds to the list of thevalues in the same order.

The Dirichlet boundary conditions are used to fill the matrix. These conditions are easy to understand and coherent with the modelization of real-life objects with a 3D matrix.

Dirichlet boundary conditions:
$$u(x)=0 \text{ on } \partial \Omega$$

```python
valtab = np.zeros(nnz)
rowtab = np.zeros(nnz)
coltab = np.zeros(n + 1)

index = 0
for i in range(n) :
    coltab[ i ] = index

    if i-ncol ** 2 >= 0 : #in front
    	valtab[ index ] = - 1
        rowtab[ index ] = i - ncol ** 2
        index += 1

    if i % ncol ** 2 >= ncol : #top
        valtab[ index ] = - 1
        rowtab[ index ] = i - ncol
        index += 1

    if i % ncol != 0 : #left
        valtab[ index ] = - 1
        rowtab[ index ] = i - 1
        index += 1

    valtab[ index ] = 6 - kh ** 2 #diag
    rowtab[ index ] = i
    index += 1

    if (i + 1) % ncol != 0 : #right
        valtab[ index ] = - 1
        rowtab[ index ] = i + 1
        index += 1

    if i % ncol ** 2 < ncol ** 2 - ncol : #bottom
        valtab[ index ] = - 1
        rowtab[ index ] = i + ncol
        index += 1

    if (i + ncol ** 2) < n : #back
        valtab[ index ] = - 1
        rowtab[ index ] = i + ncol ** 2
        index += 1

coltab[ n ] = nnz

# Loads the sparse matrix from HB driver
cscA = csc_matrix((valtab, rowtab, coltab), shape = (n, n))/ h ** 2

spmA = spm.spmatrix( cscA )
spmA.printInfo()
```

## Resolution of the Helmholtz equation with PaStiX

The following step is optional, the goal is to make sure that the Intel MKL library is loaded prior to the PaStiX library.

```python
# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))
```

You can normalize the matrix with PaStiX. It is possible to choose the norm used, the default one is Frobenius.

```python
# Scales A for low-rank: A / ||A||_f
norm = spmA.norm()
#pastix_norm = pastix.normtype.Inf
spmA.scale( 1. / norm )
```

You can make your own right-hand-side *b* if you want, but you can also generate it randomly. In this example *b* is generated randomly in order to test if the program returns the correct value. The multiple-right-hand-side is avalaible in PaStiX, it is determined by the variable *nrhs*. Here we use a single-right-hand-side (*nrhs=1*).

```python
# Generates the right-hand-sides b and x (from A * x = b)
nrhs = 1
b = spmA.genRHS( spm.rhstype.RndB, nrhs, getx=False )
x = b.copy()
```

The parameters of PaStiX are initialized with the following command. It is possible to change these parameters, for example: the factorization type or the scheduler. This step is optional, the full list of parameters is available in the documentation: [iparm](https://solverstack.gitlabpages.inria.fr/pastix/group__pastix__api.html#ga6ec72bed3a6ae6acc0291aaa0c94da1c) and [dparm](https://solverstack.gitlabpages.inria.fr/pastix/group__pastix__api.html#ga4c5b7d7a7cd48f69cc27ae29da6b275a).


```python
# Initializes the parameters to their default values
iparm, dparm = pastix.initParam()
# Starts-up PaStiX
pastix_data = pastix.init( iparm, dparm )

# Changes the scheduler
iparm[pastix.iparm.scheduler] = pastix.scheduler.Sequential
# Changes the factorization type
iparm[pastix.iparm.factorization] = pastix.factotype.LU
```

Everything is set up, in order to solve the equation, you need to call the four tasks of PaStiX:

```python
# Performs the analyze
pastix.task_analyze( pastix_data, spmA )

# Performs the numerical factorization
pastix.task_numfact( pastix_data, spmA )

# Performs the solve
pastix.task_solve( pastix_data, spmA, b )

# Perfoms the refinement on the solution x0 = b and returns the final solution in x (this task is optional)
pastix.task_refine( pastix_data, spmA, b, x )
```

## Visualization of the solution

In order to visualize the solution *x* given by PaStiX, you need to create a CSV file.

First, the solution has to be converted into a mesh, with the *ArrayToMesh3D* function:

```python
def ArrayToMesh3D(u, nrow, ncol, nprof):
  U = []
  for z in range(nprof):
    U_prof = []
    for i in range(nrow):
      U_line = []
      for j in range(ncol):
        U_line.append(u[z * ncol ** 2 + i * ncol + j])
      U_prof.append(U_line)
    U.append(U_prof)
  return np.array(U)

X = ArrayToMesh3D( x, nrow, ncol, nprof)
```

Now the corresponding lists and the CSV file can be created. You need to put "x coord, y coord, z coord" at the beginning of the file, so Paraview can read the file correctly.

```python
row=[]
for k in range(ncol):
    row = row + [k] * ncol
row = row * ncol
col = ([k for k in range(ncol)] * nrow) * nprof
prof = []
for k in range(ncol):
    prof = prof + [k] * ncol**2

val = []
for k in range(nprof):
    for i in range(ncol):
        for j in range(nrow):
            val.append(X[i][j][k])

L = np.zeros( n )
for i in range( n ) :
    L[i] = val[i]

with open("test.csv.0", "w") as fich:
    fich.write("x coord, y coord, z coord, scalar\n")
    for k in range(len(row)):
        l = row[k], col[k], prof[k], val[k]
        fich.write(f"{row[k]}, {col[k]}, {prof[k]}, {L[k]}\n")

```

For this last step, you need to open the Python Shell of Paraview (if the shell isn't open in Paraview, you have to search "View" in the tool bar and then select "Python Shell"). Paraview has an option called "Trace Python", which transcribes the actions of the interface into a Python program. The following program is the result of this trace.

```python
import paraview
from paraview.simple import *
paraview.compatibility.major = 5
paraview.compatibility.minor = 10

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
testcsv = CSVReader(registrationName='test.csv', FileName=['"PATH TO COMPLETE"'])

UpdatePipeline(time=0.0, proxy=testcsv)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(registrationName='TableToPoints1', Input=testcsv)
tableToPoints1.XColumn = ' scalar'
tableToPoints1.YColumn = ' scalar'
tableToPoints1.ZColumn = ' scalar'

# properties modified on tableToPoints1
tableToPoints1.XColumn = 'x coord'
tableToPoints1.YColumn = ' y coord'
tableToPoints1.ZColumn = ' z coord'

UpdatePipeline(time=0.0, proxy=tableToPoints1)

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(registrationName='Delaunay3D1', Input=tableToPoints1)

UpdatePipeline(time=0.0, proxy=delaunay3D1)
RenderAllViews()
```

To visualize the solution, you have to split the window, select "Render View" and click on the eye on the left of the screen. At the top of the toolbar at the top, change the option "Solid Color" into "Scalar".

You should obtain the following image corresponding to the solution of the Helmholtz equation in 3D with Dirichlet boundary conditions:![Modelization of the Helmholtz solution with Paraview](https://notes.inria.fr/uploads/upload_86c2cad54f4574fa64fe73605d07d907.png)
