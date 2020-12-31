#!/bin/bash

PATHS_TO_RUN=( "." "mapping_to_reference" "denovo_assembly" "diff_expr")

NUMBER_OF_CORES=2

ARGUMENTS="-p -j ${NUMBER_OF_CORES} --dryrun"

########################################################################

BASE_PATH=$(pwd)

for path in "${PATHS_TO_RUN[@]}"
do
  cd ${path} || exit
  echo "In path: ${path}"
  snakemake ${ARGUMENTS}
  cd ${BASE_PATH} || exit
done
