#!/bin/bash
# a script to estimate pathwise completeness of a custom butyrate biosynthesis module in Lachnospiraceae genomes 
#   using input annotation sets from each annotation tool
# designed to work with the output of the Snakemake workflow
# creates a matrix of the results wherein each row describes a genome and each column describes an annotation tool
# USAGE: bash generate_butyrate_completeness_matrix.sh FOLDER_OF_ANNOTATION_DATA
# the only argument to this script should be the folder where your annotation data is stored, organized by annotation tool

# exit upon error
set -e

# INPUT VARIABLES
if [ -z "$1" ]; then
  echo "ERROR. Must provide input folder with annotation data as the first argument to this script."
  exit 1
fi
FUNCTIONS_OUTPUT_DIR=$1

# INPUT DATA
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # location of the current script
USER_MODULES_FOLDER="${SCRIPT_DIR}/../../resources/Butanoate_module/"
GENOMES_FILE="${SCRIPT_DIR}/../../resources/lachno.tsv"

# create user-defined modules database with the butyrate biosynthesis module
# modules file must be stored at this relative path: `${USER_MODULES_FOLDER}/modules/BUTANOATE``
anvi-setup-user-modules -u $USER_MODULES_FOLDER

# convert functions output from each tool to enzymes-txt input files for anvi-estimate-metabolism
# relies on the data being stored at these relative paths: `${FUNCTIONS_OUTPUT_DIR}/{TOOL}/functions/default/{GENOME_ACCESSION}.tsv`
# also relies on the input file $GENOMES_FILE containing the accessions of the Lachnospiraceae genomes
mkdir -p ${FUNCTIONS_OUTPUT_DIR}/enzymes-txt-files

# make enzyme-txt input files for each Lachnospiraceae genome (from each annotation tool)
for tool in anvio kofamscan microbeannotator; do \
  echo "Generating enzymes-txt files from $tool annotations for Lachnospiraceae genomes...."; \
  mkdir -p ${FUNCTIONS_OUTPUT_DIR}/enzymes-txt-files/$tool
  mkdir -p ${FUNCTIONS_OUTPUT_DIR}/enzymes-txt-files/$tool/default
  while read acc; do \
    python ${SCRIPT_DIR}/functions_output_to_enzymes_txt.py ${FUNCTIONS_OUTPUT_DIR}/$tool/functions/default/${acc}.tsv \
      ${FUNCTIONS_OUTPUT_DIR}/enzymes-txt-files/$tool/default/${acc}_enzymes.txt \
      $tool; \
  done < <(tail -n+2 $GENOMES_FILE | cut -f 1 ); \
done

# estimate completeness of the butyrate biosynthesis pathways
for tool in anvio kofamscan microbeannotator; do \
  echo "Estimating butyrate metabolism from $tool annotations for Lachnospiraceae genomes...."; \
  mkdir -p ${FUNCTIONS_OUTPUT_DIR}/BUTANOATE_MODULE_OUTPUT
  mkdir -p ${FUNCTIONS_OUTPUT_DIR}/BUTANOATE_MODULE_OUTPUT/$tool
  mkdir -p ${FUNCTIONS_OUTPUT_DIR}/BUTANOATE_MODULE_OUTPUT/$tool/default
  while read acc; do \
    anvi-estimate-metabolism --enzymes-txt ${FUNCTIONS_OUTPUT_DIR}/enzymes-txt-files/$tool/default/${acc}_enzymes.txt \
      -u $USER_MODULES_FOLDER \
    --only-user-modules \
    -O ${FUNCTIONS_OUTPUT_DIR}/BUTANOATE_MODULE_OUTPUT/$tool/default/${acc} \
    --include-zeros; \
  done < <(tail -n+2 $GENOMES_FILE | cut -f 1 ); \
done

# finally, combine the results into one matrix
python gen_butyrate_matrix.py $GENOMES_FILE ${FUNCTIONS_OUTPUT_DIR}/BUTANOATE_MODULE_OUTPUT/
