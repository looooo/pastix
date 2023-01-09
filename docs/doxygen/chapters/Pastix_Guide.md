
## Basic testings

Be sure that you set up the correct environment and allocate some
resources. Note that we need only one node for those experiments, so you
can exit the previous allocation on multiple nodes, and do a new one
with a single node.

```bash
    source env-pastix.sh
    qsub -I -l select=1:ncpus=16,walltime=01:00:00
```

To get the collection of examples provided by PaStiX (see
<https://solverstack.gitlabpages.inria.fr/pastix/group__pastix__examples.html>).

To run a simple example on a 3D Laplacian of size \f$ 10 \times 10 \times 10 \f$, on one process:

```bash
    ./simple -lap 10:10:10 -t 16
```

`-t 16` specifies that you want 16 threads.

This examples can all read matrices in Harwell Boeing (`.rsa` or `.rua`),
Matrix Market (`.mtx`), or IJV format (`.ijv`), just be careful to give
the option that specifies the correct format for the file. You can also
use a simple 2D Laplacian matrix generator to the test with the
`-lap size` option. The `-h` option will give you all the options
available in examples.

## Choose ordering library: Scotch or Metis

It is possible to give to the PaStiX solver any pre-computed ordering by
providing the permutation array of the unknowns. However, it is
recommended to use the ordering step from the solver to reduce the
fill-in of the matrix during the factorization and reduce the memory
overhead of the solver. The graph partitioning library can be chosen at
runtime if and only if the solver has been compiled with multiple
choices. If it has been compiled with Scotch and Metis , both will be
available at runtime.

*Remark: with the examples `simple`, or `analyze`, you can switch from
Scotch to METIS with `–ord metis` or `–ord scotch` option.*

## Low-rank compression

PaStiX 6.0 offers low-rank compression to either accelerate the
factorization and solve steps, and/or reduce the memory peak of the
solver. Three different schemes are possible through the parameter
`iparm_compress_when`:

1.  `PastixCompressNever` (0). This is the default behavior with no
    compression at all.

2.  `PastixCompressBegin` (1). This is the behavior that will saves the
    most on the memory peak. However, it will slow down the
    factorization, and this phenomenon will increase when the block
    sizes are larger

3.  `PastixCompressEnd` (2). This is the behavior that will help to save
    time.

Low -rank compression is also strongly dependent from other parameters:

  1. `iparm_compress_method` the will define the family of kernels to use
  for the compression (SVD: `PastixCompressMethodSVD (0)`, or
  Rank-Revealing QR: `PastixCompressMethodRRQR (1)`);
  2. `dparm_compress_tolerance` that defines the tolerance of the
  compression;
  3. `iparm_compress_min_width` that defines the minimal size of required
  to compress a supernode;
  4. and `iparm_compress_min_height` that defines the minimal height required
  for an off-diagonal block to be compressed (when minimal width is
  respected).

Larger are the sizes of the blocks, more memory gain it will induce.
However, with the Minimal Memory (Begin) scenario, it will also increase
the complexity of the update kernels and slow down the computations.
This is not the case for the Just-in-Time (End) scenario.

## Runtime support

PaStiX 6.0 includes support for internal schedulers, but also to use
external runtime supports such as PaRSEC and StarPU to offer access to
accelerators. You can test the runtime support by using the `-s` option
which takes as value: $0$ for sequential scheduling, $1$ for internal
static scheduling, $2$ for PaRSEC, and $3$ for StarPU.
