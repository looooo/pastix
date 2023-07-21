#!/usr/bin/env bash
###
#
#  @file coverage.sh
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.0
#  @author Mathieu Faverge
#  @date 2023-03-22
#
###
#
# Group all the coverage files together and display the summary
# to produce the final pastix.lcov file.
#
###
export INPUT_FILES=""
for name in $( ls -1 pastix-*.lcov );
do
  lcov --remove $name '*spm/*' '*parsec/*trf_sp*' -q --output-file /tmp/$name
  mv /tmp/$name .;
  export INPUT_FILES="$INPUT_FILES -a $name";
done
lcov $INPUT_FILES -o pastix.lcov
lcov --summary pastix.lcov
lcov_cobertura pastix.lcov --output pastix_coverage.xml
gcovr --xml-pretty --exclude-unreachable-branches --print-summary -o coverage.xml --root ${CI_PROJECT_DIR}
