#!/usr/bin/python
# a script for comparing modules output from `anvi-estimate-metabolism` (for one genome) to find
#  pathways encoded in operons that are more complete when using anvi'o annotations than when using kofamscan annotations
# input: two modules mode output files for one genome, one from anvi'o and one from kofamscan
# output: a table of pathways in operons that are more complete using anvi'o
# usage: python find_and_compare_operons_in_modules.py ANVIO_MODULES_FILE KOFAMSCAN_MODULES_FiLE [OUTPUT_FILE]

import os
import sys
import pandas as pd

# COMMAND LINE PARAMETERS
if len(sys.argv) > 2:
    input_1 = sys.argv[1]
    input_2 = sys.argv[2]
else:
    print("Not enough input files. Please provide modules mode files from, i.e. (1) anvi'o and (2) kofamscan.")
if len(sys.argv) > 3:
    output_file = sys.argv[3]
else:
    output_file = os.path.dirname(anvio_input) + "anvio_vs_kofamscan_operons.txt"

## INTERNAL PARAMETERS
MAX_GENE_ID_DIFFERENCE_FOR_OPERON = 5   # how distant are subsequent gene calls allowed to be to consider them in an 'operon'?
FRACTION_SMALL_DIFFERENCES_FOR_OPERON = 0.5  # required proportion of gene calls that are close enough to each other to predict an 'operon'

## GUESS THE TOOL NAMES FROM THE INPUT PATHS (for operons in genomes table output)
split_input_path_1 = input_1.split('/')
TOOL_1 = split_input_path_1[split_input_path_1.index('metabolism') - 1]
split_input_path_2 = input_2.split('/')
TOOL_2 = split_input_path_2[split_input_path_2.index('metabolism') - 1]

## FUNCTIONS
def is_operon(gene_call_list):
    """Given a list of gene calls in the module, returns True if we predict it is encoded in an operon.
    
    The input gene_call_list should be a string with comma-separated integer values.
    """

    # we consider a module to be encoded in an operon if *most* of its component genes are located close to each other on the chromosome.
    # we can use the gene caller IDs to roughly identify genes that are close to each other (their IDs will be similar in value)
    # if we order the gene caller IDs, we can compute differences between subsequent gene IDs. If most of those differences (let's say over 50%)
    # are less than a certain value (let's say ~5 to allow for missing annotations), then we predict that it is encoded in an operon
    gcs = [int(x) for x in gene_call_list.split(',')]
    gcs = sorted(gcs)
    if len(gcs) < 2:
        return False
    diffs = [gcs[i+1] - gcs[i] for i in range(0, len(gcs) - 1)]
    small_enough_diffs = [d for d in diffs if d <= MAX_GENE_ID_DIFFERENCE_FOR_OPERON]

    if len(small_enough_diffs) / len(diffs) > FRACTION_SMALL_DIFFERENCES_FOR_OPERON:
        return True
    else:
        return False


a_df = pd.read_csv(input_1, sep="\t", index_col=0)
k_df = pd.read_csv(input_2, sep="\t", index_col=0)

## identify modules that are more complete with anvi'o annotations than with kofamscan's
mods_in_anvio = set(a_df[a_df.pathwise_module_completeness > 0].index.tolist())
mods_in_ks = set(k_df[k_df.pathwise_module_completeness > 0].index.tolist())

mods_to_consider = set([])
for mod in mods_in_anvio:
    if (mod in mods_in_ks and a_df.loc[mod, "pathwise_module_completeness"] > k_df.loc[mod, "pathwise_module_completeness"]) or \
        mod not in mods_in_ks:
            mods_to_consider.add(mod)

## identify subset of modules in potential operons
operon_mods = [mod for mod in mods_to_consider if is_operon(a_df.loc[mod, "gene_caller_ids_in_module"])]
cols_to_include = ['genome_name', 'pathwise_module_completeness', 'gene_caller_ids_in_module', 'enzyme_hits_in_module']

a_operons = a_df[a_df.index.isin(operon_mods)][cols_to_include]
a_operons.rename(columns={'pathwise_module_completeness': 'anvio_completeness', 'gene_caller_ids_in_module': 'anvio_gene_caller_ids',
                          'enzyme_hits_in_module': 'anvio_enzyme_hits'}, inplace=True)

k_operons = k_df[k_df.index.isin(operon_mods)][cols_to_include[1:]]
k_operons.rename(columns={'pathwise_module_completeness': 'kofamscan_completeness', 'gene_caller_ids_in_module': 'kofamscan_gene_caller_ids',
                          'enzyme_hits_in_module': 'kofamscan_enzyme_hits'}, inplace=True)

op_table = a_operons.join(k_operons, on='module')
print(f"Output File: {output_file}")
op_table.to_csv(output_file, sep="\t")