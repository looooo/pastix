#!/usr/bin/env bash
#
#  @file run.sh
#
#  @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.4.0
#  @author Tony Delarue
#  @author Florent Pruvost
#  @date 2024-07-05
#
# This script build the correct environment to benchmark PaStiX.
#
set -x

# Unset the binding environment of the CI for this specific case
unset STARPU_MPI_NOBIND
unset STARPU_WORKERS_NOBIND

# Set it to current directory if not run through the CI
CI_PROJECT_DIR=${CI_PROJECT_DIR:-.}

# to avoid a lock during fetching pastix branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# save guix commits
#guix describe --format=json > guix.json
guix time-machine -C ./tools/bench/guix-channels.scm -- describe --format=json > guix.json

# define env var depending on the node type
if [ "$NODE" = "bora" ]
then
  export SLURM_CONSTRAINTS="bora,omnipath"
  export PASTIX_BUILD_OPTIONS="-DPASTIX_WITH_MPI=ON -DCMAKE_BUILD_TYPE=Release -DPASTIX_WITH_STARPU=ON -DPASTIX_WITH_PARSEC=ON"
  export MPI_OPTIONS="--map-by ppr:1:node:pe=36"
  GUIX_MANIFEST="./tools/bench/guix-manifest-openmpi.scm"
elif [ "$NODE" = "sirocco" ]
then
  export SLURM_CONSTRAINTS="sirocco,omnipath,v100"
  export PASTIX_BUILD_OPTIONS="-DPASTIX_WITH_MPI=ON -DPASTIX_WITH_CUDA=ON -DCMAKE_BUILD_TYPE=Release -DPASTIX_WITH_STARPU=ON -DPASTIX_WITH_PARSEC=ON"
  export MPI_OPTIONS="--map-by ppr:1:node:pe=40"
  export LD_PRELOAD="/usr/lib64/libcuda.so"
  GUIX_MANIFEST="./tools/bench/guix-manifest-openmpi-cuda.scm"
else
  echo "$0: Please set the NODE environnement variable to bora or sirocco."
  exit -1
fi

export STARPU_HOSTNAME=$NODE
export PARSEC_HOSTNAME=$NODE

# Create envrionment and submit jobs
exec guix time-machine -C ./tools/bench/guix-channels.scm -- shell --pure \
       --preserve=PLATFORM \
       --preserve=NODE \
       --preserve=LD_PRELOAD \
       --preserve=^CI \
       --preserve=proxy$ \
       --preserve=^SLURM \
       --preserve=^JUBE \
       --preserve=^MPI \
       --preserve=^STARPU \
       --preserve=^PARSEC \
       --preserve=^PASTIX \
       -m $GUIX_MANIFEST \
       -- /bin/bash --norc ./tools/bench/plafrim/slurm.sh
err=$?

# clean tmp
rm -rf /tmp/guix-$$

# exit with error code from the guix command
exit $err
