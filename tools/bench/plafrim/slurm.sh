#!/usr/bin/env bash
#
#  @file slurm.sh
#
#  @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.0
#  @author Tony Delarue
#  @author Florent Pruvost
#  @date 2021-08-25
#
# This script launch the slurm jobs to benchmark PaStiX.
#

echo "######################### PaStiX job benchmarks #########################"

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

function wait_completion {
    # Wait for completion of jobs
    echo "JOB_LIST $JOB_LIST"
    while [ "$NJOB" -gt 0 ]
    do
        for JOB in $JOB_LIST
        do
            IS_JOB_IN_QUEUE=`squeue |grep "$JOB"`
            if [[ -z "$IS_JOB_IN_QUEUE" ]]
            then
                NJOB=$[NJOB-1]
                JOB_LIST=`echo $JOB_LIST | sed "s#$JOB##"`
                echo "JOB $JOB finished"
            else
                echo "$IS_JOB_IN_QUEUE"
            fi
        done
        sleep 30
    done
}

# Parameters of the Slurm jobs
TIME=04:00:00
PART=routage
NP=$SLURM_NP
CORES=$SLURM_CORES
CONS=$SLURM_CONSTRAINTS
EXCL=

# Submit jobs
NJOB=0
JOB_ID=`JOB_NAME=pastix\-$NODE\-$MPI\-$NP && sbatch --job-name="$JOB_NAME" --output="$JOB_NAME.out" --error="$JOB_NAME.err"\
        -N $NP --time=$TIME --partition=$PART --constraint=$CONS --exclude=$EXCL --exclusive \
        -n $NP -c $CORES $CI_PROJECT_DIR/tools/bench/pastix_guix.sh | sed "s#Submitted batch job ##"`
if [[ -n "$JOB_ID" ]]
then
    JOB_LIST="$JOB_LIST $JOB_ID"
    NJOB=$[NJOB+1]
fi

# Wait for completion of jobs
wait_completion

echo "####################### End PaStiX job benchmarks #######################"

exit 0
