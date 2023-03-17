#!/usr/bin/env bash
###
#
#  @file analysis.sh
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.0
#  @author Mathieu Faverge
#  @date 2022-08-09
#
###
PROJECT=pastix

# Performs an analysis of the project source code:
# - we consider to be in the project's source code root
# - we consider having coverage files $PROJECT_*.cov in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

TOOLSDIR=$(dirname $0)
$TOOLSDIR/find_sources.sh $BUILDDIR

# Generate coverage xml output
$TOOLSDIR/coverage.sh

# to get it displayed and captured by gitlab to expose the badge on the main page
lcov --summary $PROJECT.lcov | tee $PROJECT-gcov.log

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"
export UNDEFINITIONS="$UNDEFINITIONS -UPARSEC_PROF_DRY_BODY -UPARSEC_PROF_DRY_RUN -UPARSEC_PROF_TRACE -UPARSEC_PROF_GRAPHER -UPARSEC_SIM -DPINS_ENABLE -UBUILD_PARSEC"
export UNDEFINITIONS="$UNDEFINITIONS -UPARSEC_DEBUG_NOISIER -UPARSEC_DEBUG_PARANOID -UPARSEC_DEBUG_HISTORY -UPARSEC_C_PARSEC_HAVE_VISIBILITY"
export UNDEFINITIONS="$UNDEFINITIONS -UBUILDING_STARPU -UPASTIX_WITH_STARPU_DIST"
export UNDEFINITIONS="$UNDEFINITIONS -UNAPA_SOPALIN"

# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS} --file-list=./filelist-c.txt 2> ${PROJECT}_cppcheck.xml

# Set the default for the project key
SONARQUBE_PROJECTKEY=${SONARQUBE_PROJECTKEY:-topal:$CI_PROJECT_NAMESPACE:$PROJECT}

ls $BUILDDIR/*.json

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectKey=$SONARQUBE_PROJECTKEY
sonar.projectDescription=Parallel Sparse direct Solver
sonar.projectVersion=6.3.0

sonar.scm.disabled=false
sonar.scm.provider=git
sonar.scm.exclusions.disabled=true

sonar.sources=$BUILDDIR, bcsc, blend, common, example, graph, include, kernels, order, refinement, sopalin, symbol, test
sonar.inclusions=`cat filelist.txt | grep -v spm | xargs echo | sed 's/ /, /g'`
sonar.sourceEncoding=UTF-8
sonar.cxx.file.suffixes=.h,.c
sonar.cxx.errorRecoveryEnabled=true
sonar.cxx.gcc.encoding=UTF-8
sonar.cxx.gcc.regex=(?<file>.*):(?<line>[0-9]+):[0-9]+:\\\x20warning:\\\x20(?<message>.*)\\\x20\\\[(?<id>.*)\\\]
sonar.cxx.gcc.reportPaths=$PROJECT_build*.log
sonar.cxx.xunit.reportPaths=*.junit
sonar.cxx.cobertura.reportPaths=*.cov
sonar.cxx.cppcheck.reportPaths=${PROJECT}_cppcheck.xml
sonar.cxx.clangsa.reportPaths=$BUILDDIR/analyzer_reports/*/*.plist
sonar.cxx.jsonCompilationDatabase=$BUILDDIR/compile_commands-seq.json, $BUILDDIR/compile_commands-mpi.json
EOF
echo "====== sonar-project.properties ============"
cat sonar-project.properties
echo "============================================"

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log
