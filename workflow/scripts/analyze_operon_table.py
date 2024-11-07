## script for analyzing the operon table
## filters the operon example table to get some good examples for the supplementary figure, and computes some numbers about them

import pandas as pd

import anvio
import argparse
from anvio import kegg

MIN_GENOMES_WITH_MODULE_AS_OPERON = 30 # we want to see the module encoded in operon-like structure in MORE than this many genomes

# INPUT FILES (OTHER USERS MAY NEED TO CHANGE PATHS HERE)
mod_db_path = "KEGG_12_15_2023/MODULES.db"                     # PATH TO MODULES DATABASE THAT WE USED FOR OUR ANALYSIS
a_vs_k_operon_example_file = "anvio-vs-kofamscan-operon-examples.txt" # output of find_operon_examples.sh
a_vs_m_operon_example_file = "anvio-vs-microbeannotator-operon-examples.txt"
operon_table_file = "operons_in_genomes.txt"                   # output of find_operon_examples.sh

# OUTPUT FILES
a_vs_k_sorted_out = "anvio-vs-kofamscan-operons-sorted.txt"
a_vs_k_annotation_differences_out = "example_modules_KO_annotation_differences.txt"

# take the original output of find_operon_examples.sh and 
# sort it by (1) magnitude of completeness difference between the two tools and (2) operon length
df = pd.read_csv(a_vs_k_operon_example_file, sep="\t")
df['diff'] = df.anvio_completeness - df.kofamscan_completeness
df['len'] = df.anvio_gene_caller_ids.str.split(',').str.len()
df.sort_values(['diff', 'len'], ascending=[False, False], inplace=True)
df.to_csv(a_vs_k_sorted_out, sep="\t", index=False)

args = argparse.Namespace()
db = kegg.ModulesDatabase(mod_db_path, args)

# count the total number of modules that are encoded in operon-like structure
ops = pd.read_csv(operon_table_file, sep="\t")
modules_x_operons_dict = {}
for mod in db.get_all_modules_as_list():
    modules_x_operons_dict[mod] = {}
    for tool in ['anvio', 'kofamscan', 'microbeannotator']:
        modules_x_operons_dict[mod][f"num_genomes_it_is_operon_in_{tool}"] = ops[(ops.modules_in_operon.str.contains(mod)) & (ops.tool == tool)].shape[0]
mod_x_ops_df = pd.DataFrame().from_dict(modules_x_operons_dict).T
operon_like_mods = mod_x_ops_df[(mod_x_ops_df.num_genomes_it_is_operon_in_anvio > MIN_GENOMES_WITH_MODULE_AS_OPERON) | \
             (mod_x_ops_df.num_genomes_it_is_operon_in_kofamscan > MIN_GENOMES_WITH_MODULE_AS_OPERON) | \
             (mod_x_ops_df.num_genomes_it_is_operon_in_microbeannotator > MIN_GENOMES_WITH_MODULE_AS_OPERON)].index.to_list()
num_operon_like_mods = len(operon_like_mods)
print(f"Number of operon-like modules (in at least {MIN_GENOMES_WITH_MODULE_AS_OPERON + 1} genomes): {num_operon_like_mods}")

df2 = pd.read_csv(a_vs_m_operon_example_file, sep="\t")
diff_annotated_by_anvio = set([m for m in df.module.to_list() if m in operon_like_mods]).union(set([m for m in df2.module.to_list() if m in operon_like_mods]))
num_diff_annotated_by_anvio = len(diff_annotated_by_anvio)
print(f"The following {num_diff_annotated_by_anvio} modules are more complete using anvi'o annotations:\n{','.join(diff_annotated_by_anvio)}")

df = pd.read_csv(a_vs_k_sorted_out, sep="\t") # this file was generated above in this script
# how many of these modules appear to be encoded in an operon in MORE than N genomes?
df.module.value_counts()[df.module.value_counts() > MIN_GENOMES_WITH_MODULE_AS_OPERON]
good_module_examples = df.module.value_counts()[df.module.value_counts() > MIN_GENOMES_WITH_MODULE_AS_OPERON].index.to_list()

example_mod_table = pd.DataFrame()
for mod in good_module_examples:
    num_genomes_diff_completeness_for_mod = df[df.module == mod].shape[0]
    kos_in_mod = db.get_kos_in_module(mod)
    mod_dict = {k: {'module': mod, 'num_genomes_with_differential_module_completeness': num_genomes_diff_completeness_for_mod} for k in kos_in_mod}
    for k in kos_in_mod:
        num_genomes_anvio_finds_ko = df[df.module == mod].anvio_enzyme_hits.str.contains(k).value_counts()[True] if True in df[df.module == mod].anvio_enzyme_hits.str.contains(k).value_counts() else 0
        num_genomes_kofamscan_finds_ko = df[df.module == mod].kofamscan_enzyme_hits.str.contains(k).value_counts()[True] if True in df[df.module == mod].kofamscan_enzyme_hits.str.contains(k).value_counts() else 0

        mod_dict[k]['num_genomes_anvio_finds_ko'] = num_genomes_anvio_finds_ko
        mod_dict[k]['num_genomes_kofamscan_finds_ko'] = num_genomes_kofamscan_finds_ko

    mod_df = pd.DataFrame.from_dict(mod_dict).T
    example_mod_table = pd.concat([example_mod_table, mod_df])

example_mod_table.index.rename('KO', inplace=True)
example_mod_table.to_csv(a_vs_k_annotation_differences_out, sep="\t")