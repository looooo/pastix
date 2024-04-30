#!/usr/bin/env bash
#
#  @file pastix_guix.sh
#
#  @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.2
#  @author Tony Delarue
#  @author Florent Pruvost
#  @date 2023-07-21
#
# This script compiles, runs and analyzes PaStiX experiments
# for benchmarking purpose.
#
set -ex

# Configure and build PaStiX
if [ -d build-$NODE-$MPI ]; then
  rm build-$NODE-$MPI -r
fi
cmake -B build-$NODE-$MPI $PASTIX_BUILD_OPTIONS
cmake --build build-$NODE-$MPI -j20 --verbose
export PASTIX_BUILD=$PWD/build-$NODE-$MPI

# clean old benchmarks
if [ -d tools/bench/$PLATFORM/results ]; then
  rm tools/bench/$PLATFORM/results -r
fi

# Execute jube benchmarks
jube run tools/bench/$PLATFORM/pastix.xml --tag $JUBE_RUN --include-path tools/bench/$PLATFORM/parameters/$NODE --id $JUBE_ID
#jube run tools/bench/$PLATFORM/pastix-test.xml --tag $JUBE_RUN --include-path tools/bench/$PLATFORM/parameters/$NODE --id $JUBE_ID

# jube analysis
jube analyse tools/bench/$PLATFORM/results --id $JUBE_ID

# jube report
jube result tools/bench/$PLATFORM/results --id $JUBE_ID > pastix_$JUBE_ID.csv
cat pastix_$JUBE_ID.csv

# send results to the elasticsearch server
#ls guix.json
python3 tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "pastix" pastix_$JUBE_ID.csv
#python3 tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "pastix-test" pastix_$JUBE_ID.csv
