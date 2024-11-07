#!/bin/bash
# a driver script to run find_and_compare_operons_in_modules.py on all genomes
# designed to work with the output of the Snakemake workflow
# USAGE: bash find_operon_examples.sh FOLDER_OF_ANNOTATION_DATA
# the only argument to this script should be the folder where your annotation data is stored, organized by annotation tool

# exit upon error
set -e

# INPUT VARIABLES
if [ -z "$1" ]; then
  echo "ERROR. Must provide input folder with annotation data as the first argument to this script."
  exit 1
fi
FUNCTIONS_OUTPUT_DIR=$1
OPERONS_IN_GENOMES="${FUNCTIONS_OUTPUT_DIR}/operons_in_genomes.txt"
GENOMES_LIST_FILE="${FUNCTIONS_OUTPUT_DIR}/genome_list.txt"  # used as input to KO_annotation_counts.py later

TOOL_DIR_1="anvio"
TOOL_DIR_2="kofamscan"
TOOL_DIR_3="microbeannotator"
OUTPUT_TABLE_1="${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_1}-vs-${TOOL_DIR_2}-operon-examples.txt"
OUTPUT_TABLE_2="${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_1}-vs-${TOOL_DIR_3}-operon-examples.txt"

# prep for output (remove previously-existing output)
rm -f $OUTPUT_TABLE_1 $OUTPUT_TABLE_2 $GENOMES_LIST_FILE
echo -e "genome\ttool\tmodules_in_operon" > $OPERONS_IN_GENOMES # this output will be appended to by `find_and_compare_operons_in_modules.py``

for dir in ${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_1}/metabolism/default/*/; do \
    dir=${dir%*/}      # remove the trailing "/"
    genome="${dir##*/}"
    echo $genome >> $GENOMES_LIST_FILE
    echo "Processing $genome"
    IN_1="${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_1}/metabolism/default/${genome}/anvio_estimate_metabolism_modules.txt"
    IN_2="${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_2}/metabolism/default/${genome}/anvio_estimate_metabolism_modules.txt"
    OUT="${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_1}/metabolism/default/${genome}/anvio_vs_kofamscan_operons.txt"
    python find_and_compare_operons_in_modules.py $IN_1 $IN_2 $OUT $genome $OPERONS_IN_GENOMES
    # combine with the overall table
    if [ ! -f $OUTPUT_TABLE_1 ]; then
        cat $OUT > $OUTPUT_TABLE_1
    else
        cat $OUTPUT_TABLE_1 <(tail -n+2 $OUT) > new && mv new $OUTPUT_TABLE_1
    fi

    IN_3="${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_3}/metabolism/default/${genome}/anvio_estimate_metabolism_modules.txt"
    OUT="${FUNCTIONS_OUTPUT_DIR}/${TOOL_DIR_1}/metabolism/default/${genome}/anvio_vs_microbeannotator_operons.txt"
    python find_and_compare_operons_in_modules.py $IN_1 $IN_3 $OUT $genome $OPERONS_IN_GENOMES
     # combine with the overall table
    if [ ! -f $OUTPUT_TABLE_2 ]; then
        cat $OUT > $OUTPUT_TABLE_2
    else
        cat $OUTPUT_TABLE_2 <(tail -n+2 $OUT) > new && mv new $OUTPUT_TABLE_2
    fi
done

echo "Differentially-annotated operon examples are printed to $OUTPUT_TABLE_1 and $OUTPUT_TABLE_2"
echo "Modules with operon-like structure in each genome are described in $OPERONS_IN_GENOMES"
echo "List of genome acessions is printed to $GENOMES_LIST_FILE"