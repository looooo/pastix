## How to use PaStiX with GPUs

After you have installed all the components listed in
_Installing PaStiX with GPU support_, you can now launch
PaStiX with GPUs.
The option **-g** allows you to define the number of GPU
that you want to build.
The options **-s** defines the scheduler:
* **-s 2** uses the PaRSEC scheduler.
* **-s 3** uses the StarPU scheduler.

```sh
$PASTIX_HOME_DIR/build/examples/simple -g x -s [2|3] --mm your_matrix.mtx
```

Note that GPUs have a significative impact when the number
of computations is large enough. If not, the time to transfer
datas to the GPU will be bigger than your computation time.

To have better performances, you should also do a first
factorization, for warming-up the data transfers. You can
find an example of multiple factorizations in
$PASTIX_HOME_DIR/examples/bench_facto.c

### Obtaining better performances with PaRSEC

At the root of your session, create a .parsec directory:

```sh
cd ~
mkdir .parsec/
```

In this directory, create a file named **mca-params.conf** which will
define some variables for the PaRSEC configuration:

```sh
device_show_capabilities = 1
device_show_statistics = 1            # Display how much your devices are used.
debug_verbose = 0                     # Debug mode.
debug_color = 0
device_cuda_max_streams = 1           # Number of streams per cuda device.
                                      # Above 8, no major upgrades are observed.
device_cuda_max_events_per_stream = 4 # Number of tasks pipelined in your device.
runtime_comm_short_limit = 0
```

These variables will the be used if you run PaStiX with the
PaRSEC scheduler:

```sh
$PASTIX_HOME_DIR/build/examples/bench_facto -g 1 -s 2 --mm your_matrix.mtx
```

### Obtaining better performances with StarPU.

In the same way, you can obtain better performances with StarPu thanks to
environment variables. You can find the complete list
[here](https://files.inria.fr/starpu/starpu-1.3.7/html/ExecutionConfigurationThroughEnvironmentVariables.html)

To name some of them:

```sh
STARPU_NWORKER_PER_CUDA=xx         # Number of streams per cuda device.
STARPU_CUDA_THREAD_PER_WORKER=xx   # Number of threads per cuda streams.
STARPU_CUDA_THREAD_PER_DEV=xx      # Number of threads per cuda device.
STARPU_CUDA_PIPELINE=xx            # Number of tasks pipelined in your device.
STARPU_SILENT=1                    # Will make StarPU silent. Good to use if you're
                                   # sure that StarPu didn't find conflicts with your
                                   # configuration.
```

You can then export this variables in your environment or put them
at the beginning of you command line:

```sh
STARPU_NWORKER_PER_CUDA=8 STARPU_CUDA_PIPELINE= 4 \
$PASTIX_HOME_DIR/build/examples/bench_facto -g 1 -s 3 --mm your_matrix.mtx

# OR

export STARPU_NWORKER_PER_CUDA=8
export STARPU_CUDA_PIPELINE= 4
$PASTIX_HOME_DIR/build/examples/bench_facto -g 1 -s 3 --mm your_matrix.mtx
```