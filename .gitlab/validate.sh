#!/usr/bin/env bash
###
#
#  @file validate.sh
#  @copyright 2023-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2023-09-22
#
###

# Check some metrics on sonarqube (https://sonarqube.inria.fr/sonarqube/)
# and depending on the value return 0 (success) or 1 (failure).

if [ $# -gt 0 ]; then
  METRIC=$1
fi
METRIC=${METRIC:-BUG}

if [[ -z $CI_MERGE_REQUEST_IID || -z $CI_PROJECT_NAMESPACE || -z $CI_PROJECT_NAME ]]; then
  echo "One of the variables CI_MERGE_REQUEST_IID, CI_PROJECT_NAMESPACE,
  CI_PROJECT_NAME is empty. This script must be used during a gitlab merge
  request only -> Failure."
  exit 1
fi

if [[ -z $SONARQUBE_LOGIN ]]; then
  echo "SONARQUBE_LOGIN is empty, please give a valid sonarqube user's token,
  with permissions set on the project -> Failure."
  exit 1
fi

if [[ $METRIC == "BUG" ]]; then
  BUG=`curl -u $TOKEN: -X GET "https://sonarqube.inria.fr/sonarqube/api/measures/component?component=${CI_PROJECT_NAMESPACE}%3A${CI_PROJECT_NAME}&pullRequest=${CI_MERGE_REQUEST_IID}&metricKeys=new_bugs" |jq '.component.measures[0].period.value' | sed -e "s#\"##g"`
  if [[ $BUG > 0 ]]; then
    echo "%{BUG} new bugs detected by Sonarqube -> Failure."
    exit 1
  else
    echo "No new bugs detected by Sonarqube -> Success."
    exit 0
  fi
elif [[ $METRIC == "COVERAGE" ]]; then
  COV=`curl -u $TOKEN: -X GET "https://sonarqube.inria.fr/sonarqube/api/measures/component?component=${CI_PROJECT_NAMESPACE}%3A${CI_PROJECT_NAME}&pullRequest=${CI_MERGE_REQUEST_IID}&metricKeys=new_coverage" |jq '.component.measures[0].period.value' | sed -e "s#\"##g"`
  if [[ $COV == "null" || -z $COV  ]]; then
    echo "Coverage is empty, certainly that there are no lines of new code (considered during the analysis) to compare -> Success."
  else
    if [[ $COV < 80 ]]; then
      echo "Coverage on new lines is ${COV}%, which is < 80% -> Failure."
      exit 1
    else
      echo "Coverage on new lines is ${COV}%, which is >= 80% -> Success."
      exit 0
    fi
  fi
fi
