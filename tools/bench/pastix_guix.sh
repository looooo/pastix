#!/usr/bin/env bash
#
#  @file pastix_guix.sh
#
#  @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.0
#  @author Tony Delarue
#  @author Florent Pruvost
#  @date 2021-08-25
#
# This script compiles, runs and analyzes PaStiX experiments
# for benchmarking purpose.
#

# Configure and build PaStiX
mkdir -p $CI_PROJECT_DIR/build-$NODE-$MPI
cd $CI_PROJECT_DIR/build-$NODE-$MPI
rm CMake* -rf
cmake $PASTIX_BUILD_OPTIONS ..
make -j20
export PASTIX_BUILD=$PWD

# clean old benchmarks
cd $CI_PROJECT_DIR/tools/bench/$PLATFORM/results
jube remove --force --id $JUBE_ID

# Execute jube benchmarks
cd $CI_PROJECT_DIR/tools/bench/$PLATFORM/
jube run pastix.xml --tag $JUBE_RUN --include-path parameters/$NODE --id $JUBE_ID

# jube analysis
jube analyse results --id $JUBE_ID

# jube report
jube result results --id $JUBE_ID > pastix_$JUBE_ID.csv

# send results to the elasticsearch server
cp $CI_PROJECT_DIR/guix.json .
python3 $CI_PROJECT_DIR/tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "pastix" -h $NODE -m $MPI pastix_$JUBE_ID.csv
