## How to use PaStiX with GPUs

After you have installed all the components listed in
_Installing PaStiX with GPU support_, you can now launch
PaStiX with GPUs.
The option **-g** allows you to define the number of GPUs
that you want to build.
The options **-s** defines the scheduler:
* **-s 2** uses the PaRSEC scheduler.
* **-s 3** uses the StarPU scheduler.

```sh
$PASTIX_BUILD_DIR/example/simple -g x -s [2|3] --mm your_matrix.mtx
```

Note that GPUs have a significant impact when the number
of computations is large enough. If not, the time to transfer
data to the GPU will be larger than your computation time.

Note that the performance are largely impacted on the first
factorization by the loading of the CUDA library. To better evaluate
the performance on your matrices, we recommend to do multiple runs
with the help of the ` bench_facto` example
(`$PASTIX_SRC_DIR/example/bench_facto.c`) or to do CUDA calls on fake
data before calling the solver.

### Obtaining better performance with PaRSEC

In your home directory, you can create a `~/.parsec/mca-params.conf`
file to better configure PaRSEC.

This file defines all MCA parameters of the PaRSEC library. One per line. We list
here a few parameters that can be helpful to get better performance
with PaStiX:

  * Display the compute capabilities of the discovered devices
```sh
device_show_capabilities = 1
```

  * Display device's statistics to check load balancing and data transfers
```sh
device_show_statistics = 1
```

  * Define the number of CUDA streams. This value must be `>= 3`, as
    the first two streams are reserved to data transfers from/to the
    GPUs. For example, we recommend a value of `8` for Nvidia V100,
    and larger values do not seem to improve the performance.
```sh
device_cuda_max_streams = 8
```

  * Define the number of CUDA events that are pipelined per
  stream. `4` is the common value that we usually recommend.
```sh
device_cuda_max_events_per_stream = 4
```

  * For MPI runs, we **do recommend** to disable the PaRSEC short
    communications by setting their size limit to `0`. Enabling may
    cause communication issues.
```sh
runtime_comm_short_limit = 0
```

Note that you can set all other MCA parameters of PaRSEC through this
file such as the scheduler, some limits, enabling/disabling profiling,
...
These parameters can be found by looking for `mca_param_reg_.*` in
the PaRSEC code.

Then, all these variables will automatically be used by PaStiX when it
uses the PaRSEC scheduler as in:

```sh
$PASTIX_BUILD_DIR/example/bench_facto -g 1 -s 2 --mm your_matrix.mtx
```

Note that this example uses 1 GPU (`-g 1`) and PaRSEC (`-s 2`), and
the total number of threads is let to the topology discovery by hwloc.

### Obtaining better performance with StarPU.

In the same way, you can obtain better performance with StarPU thanks to
environment variables. You can find the complete list in StarPU
[Handbook](https://files.inria.fr/starpu/starpu-1.3.7/html/ExecutionConfigurationThroughEnvironmentVariables.html).

We list here a subset of these variables that may have an impact of
the PaStiX performance:

  * The number of CUDA events that are pipelined per
  stream. `4` is the common value that we usually recommend.
```sh
STARPU_CUDA_PIPELINE=4
```

  * The number of CUDA streams. For example, we recommend a value of
    `8` for Nvidia V100, and larger values do not seem to improve the
    performance. Not that, as opposed to PaRSEC, the data transfer
    streams are not included in this value.
```sh
STARPU_NWORKER_PER_CUDA=8
```

  * One thing that may have an impact on the performance is to
    dedicate on thread per CUDA stream. However if this option is
    enabled, the total number of threads  may be too large, and you
    may require to decrease the number of CPU workers.
    Thus, this may be beneficial on an architecture with very large
    number of cores, but may degrade performance on a small architecture.
```sh
STARPU_CUDA_THREAD_PER_WORKER=[0|1]
```

To use them, you can simply export this variables in your environment or put them
at the beginning of you command line:

```sh
STARPU_NWORKER_PER_CUDA=8 STARPU_CUDA_PIPELINE=4 \
$PASTIX_BUILD_DIR/example/bench_facto -g 1 -s 3 --mm your_matrix.mtx
```

or

```sh
export STARPU_NWORKER_PER_CUDA=8
export STARPU_CUDA_PIPELINE=4
$PASTIX_BUILD_DIR/example/bench_facto -g 1 -s 3 --mm your_matrix.mtx
```

Note that these examples use a single gpu (`-g 1`) with StarPU (`-s 3`).
The total number of threads is let to the discovery of the
topology by hwloc.
