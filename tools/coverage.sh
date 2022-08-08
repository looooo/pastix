#!/usr/bin/env bash
###
#
#  @file coverage.sh
#  @copyright 2013-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.1.0
#  @author Mathieu Faverge
#  @date 2022-02-08
#
###
#
# Group all the coverage files together and display the summary
# to produce the final chameleon.lcov file.
#
###
export INPUT_FILES=""
for name in $( ls -1 pastix-*.lcov );
do
  lcov --remove $name '*spm/*' '*parsec/*trf_sp*'
     -q --output-file /tmp/$name
  mv /tmp/$name .;
  export INPUT_FILES="$INPUT_FILES -a $name";
done
lcov $INPUT_FILES -o pastix.lcov
lcov --summary pastix.lcov
