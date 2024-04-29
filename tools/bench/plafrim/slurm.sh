#!/usr/bin/env bash
#
#  @file slurm.sh
#
#  @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.2
#  @author Tony Delarue
#  @author Florent Pruvost
#  @date 2023-07-21
#
# This script launch the slurm jobs to benchmark PaStiX.
#

# Check the environment
echo $PLATFORM
echo $NODE
env |grep ^CI
env |grep ^SLURM
env |grep ^JUBE
env |grep ^MPI
env |grep ^STARPU
env |grep ^PASTIX

set -x

# Parameters of the Slurm jobs
TIME=01:30:00
PART=routage
NP=$SLURM_NP
CORES=$SLURM_CORES
CONS=$SLURM_CONSTRAINTS
EXCL=

# Submit jobs
export JOB_NAME=pastix\-$NODE\-$MPI\-$NP
sbatch --wait --job-name="$JOB_NAME" --output="$JOB_NAME.out" --error="$JOB_NAME.err" \
        -N $NP --time=$TIME --partition=$PART --constraint=$CONS --exclude=$EXCL --exclusive \
        -n $NP -c $CORES ./tools/bench/pastix_guix.sh
err=$?

cat pastix\-$NODE\-$MPI\-$NP.err
cat pastix\-$NODE\-$MPI\-$NP.out

# exit with error code from the sbatch command
exit $err
