#!/usr/bin/env bash
###
#
#  @file coverage.sh
#  @copyright 2013-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.2
#  @author Mathieu Faverge
#  @date 2023-12-07
#
###
#
# Group all the coverage files together and display the summary
# to produce the final ${PROJECT_NAME}.lcov file.
#
###
if [ $# -lt 1 ]
then
    echo "Usage $0 project_name"
fi

PROJECT_NAME=$1

INPUT_FILES=""
for name in $( ls -1 ${PROJECT_NAME}-*.lcov );
do
    lcov --remove $name '*spm/*' '*parsec/*trf_sp*' -q --output-file /tmp/$name
    mv /tmp/$name .;
    INPUT_FILES="$INPUT_FILES -a $name";
done
lcov $INPUT_FILES -o ${PROJECT_NAME}.lcov
lcov --summary ${PROJECT_NAME}.lcov
lcov_cobertura ${PROJECT_NAME}.lcov --output ${PROJECT_NAME}-coverage.xml
