#!/usr/bin/env bash
#
#  @file run.sh
#
#  @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.2.1
#  @author Tony Delarue
#  @author Florent Pruvost
#  @date 2021-07-19
#
# This script build the correct environment to benchmark PaStiX.
#

echo "######################### PaStiX environment benchmarks #########################"

# Unset the binding environment of the CI for this specific case
unset STARPU_MPI_NOBIND
unset STARPU_WORKERS_NOBIND

# to avoid a lock during fetching pastix branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# save guix commits
guix describe --format=json > guix.json

GUIX_ENV="pastix"
export MPI_OPTIONS=""
# define env var depending on the node type
if [ "$NODE" = "bora" ]
then
  export SLURM_CONSTRAINTS="bora,omnipath"
  export PASTIX_BUILD_OPTIONS="-DPASTIX_WITH_MPI=ON -DCMAKE_BUILD_TYPE=Release -DPASTIX_WITH_STARPU=ON -DPASTIX_WITH_PARSEC=ON"
  export MPI_OPTIONS="--map-by ppr:1:node:pe=36"
elif [ "$NODE" = "sirocco" ]
then
  export SLURM_CONSTRAINTS="sirocco,omnipath,v100"
  export PASTIX_BUILD_OPTIONS="-DPASTIX_WITH_MPI=OFF -DPASTIX_WITH_CUDA=ON -DCMAKE_BUILD_TYPE=Release -DPASTIX_WITH_STARPU=ON -DPASTIX_WITH_PARSEC=ON"
  export LD_PRELOAD="/usr/lib64/libcuda.so /usr/lib64/libnvidia-fatbinaryloader.so.440.33.01"
  GUIX_ENV="pastix-cuda" # Doesn't exist for the moment. Create it in guix-hpc/inria/hiepacs.scm
else
  echo "$0: Please set the NODE environnement variable to bora or sirocco."
  exit -1
fi

export STARPU_HOSTNAME=$NODE
export PARSEC_HOSTNAME=$NODE

if [ "$MPI" = "openmpi" ]
then
  GUIX_ENV_MPI=""
  GUIX_ADHOC_MPI="openssh openmpi"
fi

GUIX_ADHOC="coreutils gawk grep jube perl python python-click python-gitpython python-elasticsearch python-certifi sed slurm@19 mkl libltdl zlib"

GUIX_RULE="$GUIX_ENV $GUIX_ENV_MPI --ad-hoc $GUIX_ADHOC $GUIX_ADHOC_MPI"

# Create envrironment and submit jobs
exec guix environment --pure \
                      --preserve=PLATFORM \
                      --preserve=NODE \
                      --preserve=LD_PRELOAD \
                      --preserve=^CI \
                      --preserve=^SLURM \
                      --preserve=^JUBE \
                      --preserve=^MPI \
                      --preserve=^STARPU \
                      --preserve=^PARSEC \
                      --preserve=^PASTIX \
                      $GUIX_RULE \
                      -- /bin/bash --norc ./tools/bench/plafrim/slurm.sh

echo "####################### End PaStiX environment benchmarks #######################"

# Clean tmp
rm -rf /tmp/guix-$$
