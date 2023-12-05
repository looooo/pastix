#!/usr/bin/env bash
###
#
#  @file test.sh
#  @copyright 2023-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.1
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2023-11-29
#
###

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

cd build
source ${CI_PROJECT_DIR}/.gitlab/env.sh ${VERSION} || fatal
source ${CI_PROJECT_DIR}/install-${VERSION}/bin/pastix_env.sh || fatal
CTESTCOMMAND=`echo "ctest --output-on-failure --no-compress-output $TESTS_RESTRICTION -T Test --output-junit ../${LOGNAME}-junit.xml"`
$CTESTCOMMAND || fatal
if [[ "$SYSTEM" == "linux" ]]; then
  # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
  # so that we can only make the coverage report on the linux runner with gcc
  cd ..
  lcov --directory build --capture --output-file ${LOGNAME}.lcov
fi
