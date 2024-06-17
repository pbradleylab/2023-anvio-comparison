#!/usr/bin/python
# a script for converting functions output from multiple tools to enzymes-txt input for `anvi-estimate-metabolism`
# input: a functions output file (from one tool)
# output: an enzymes-txt file corresponding to the input functions file
# usage: python functions_output_to_enzymes_txt.py INPUT_FILE [OUTPUT_FILE] [TOOL_TYPE]

import os
import sys
import pandas as pd

# COMMAND LINE PARAMETERS
if len(sys.argv) > 1: 
    input_functions_file = sys.argv[1]
else:
    print("No input file specified")
    sys.exit(1)

if len(sys.argv) > 2: 
    output_file = sys.argv[2]
else:
    output_file = input_functions_file + "-enzymes.txt"

if len(sys.argv) > 3:
    tool = sys.argv[3]
else:
    tool = "anvio"

# METHODS to convert functions output from each tool type into an enzymes-txt table
# each method accepts a file path to the functions output as input and returns a pandas DataFrame
# for anvi'o and kofamscan, you can also choose whether to keep all hits or just the best one per gene call
def convert_anvio(file, keep_best_hit_per_gene=True):
    df = pd.read_csv(file, sep="\t")
    df = df[df.source == 'KOfam']
    if keep_best_hit_per_gene:
        # keep hit with lowest e-value per gene
        df = df.sort_values('e_value').drop_duplicates('gene_callers_id', keep='first')
    enzymes_df = df[['gene_callers_id', 'accession', 'source']]
    enzymes_df.columns = ['gene_id', 'enzyme_accession', 'source']
    return enzymes_df

def convert_kofamscan(file, keep_best_hit_per_gene=True):
    df = pd.read_csv(file, sep="\t", comment="#", \
                     names = ['passes_threshold','gene','KO','thrshld','score','evalue','definition'])
    df = df[df.passes_threshold == "*"]
    if keep_best_hit_per_gene:
        # keep hit with lowest e-value per gene
        df = df.sort_values('evalue').drop_duplicates('gene', keep='first')
    enzymes_df = df[['gene','KO']]
    enzymes_df["source"] = "KOfam"
    enzymes_df.columns = ['gene_id', 'enzyme_accession', 'source']
    return enzymes_df

def convert_microbeannotator(file):
    df = pd.read_csv(file, sep="\t")
    enzymes_df = df[df.ko_number.notna()][['query_id','ko_number']]
    enzymes_df["source"] = "KOfam"
    enzymes_df.columns = ['gene_id', 'enzyme_accession', 'source']
    return enzymes_df

print(f"Input File: {input_functions_file}")
enzymes_df = None
if tool == 'anvio':
    enzymes_df = convert_anvio(input_functions_file)
elif tool == 'kofamscan':
    enzymes_df = convert_kofamscan(input_functions_file)
elif tool == "microbeannotator":
    enzymes_df = convert_microbeannotator(input_functions_file)
else:
    print(f"ERROR. No function defined for processing tool '{tool}'.")
    sys.exit(1)

print(f"Output File: {output_file}")
enzymes_df.to_csv(output_file, sep="\t", index=False)