#!/usr/bin/python
# a script to get annotation counts for all KOs in example modules

import os
import pandas as pd

import anvio
import argparse
from anvio import kegg

# INPUT FILES (OTHER USERS MAY NEED TO CHANGE PATHS HERE)
mod_db_path = "KEGG_12_15_2023/MODULES.db" # PATH TO MODULES DATABASE THAT WE USED FOR OUR ANALYSIS
WORKFLOW_OUTPUT_PATH = "annotation"                 # PATH TO SNAKEMAKE WORKFLOW OUTPUT DIRECTORY
ks_functions_folder = os.path.join(WORKFLOW_OUTPUT_PATH, "kofamscan/functions/default/")
an_functions_folder = os.path.join(WORKFLOW_OUTPUT_PATH, "anvio/functions/default/")
ma_functions_folder = os.path.join(WORKFLOW_OUTPUT_PATH, "microbeannotator/functions/default/")
genomes_file = os.path.join(WORKFLOW_OUTPUT_PATH, "genome_list.txt") # this is output from find_operon_examples.sh

OUTPUT = "example_modules_KO_annotation_counts.txt"

def convert_anvio(file, keep_best_hit_per_gene=False):
    df = pd.read_csv(file, sep="\t")
    df = df[df.source == 'KOfam']
    if keep_best_hit_per_gene:
        # keep hit with lowest e-value per gene
        df = df.sort_values('e_value').drop_duplicates('gene_callers_id', keep='first')
    enzymes_df = df[['gene_callers_id', 'accession']]
    enzymes_df.columns = ['gene', 'KO']
    return enzymes_df

def convert_kofamscan(file, keep_best_hit_per_gene=False):
    df = pd.read_csv(file, sep="\t", comment="#", \
                     names = ['passes_threshold','gene','KO','thrshld','score','evalue','definition'])
    df = df[df.passes_threshold == "*"]
    if keep_best_hit_per_gene:
        # keep hit with lowest e-value per gene
        df = df.sort_values('evalue').drop_duplicates('gene', keep='first')
    enzymes_df = df[['gene','KO']]
    return enzymes_df

def convert_microbeannotator(file):
    df = pd.read_csv(file, sep="\t")
    enzymes_df = df[df.ko_number.notna()][['query_id','ko_number']]
    enzymes_df.columns = ['gene', 'KO']
    return enzymes_df

args = argparse.Namespace()
db = kegg.ModulesDatabase(mod_db_path, args)


genomes_list = [g.strip() for g in open(genomes_file, 'r').readlines()]

module_example_list = ["M00023", "M00026", "M00149", "M00157", "M00924", "M00925"]
ko_count_dict = {}
for mod in module_example_list:
    kos_in_mod = db.get_kos_in_module(mod)
    for k in kos_in_mod:
        ko_count_dict[k] = {'module': mod, 'total_genes_annotated': 0,
                            'all_tools': 0,
                            'anvio_and_kofamscan': 0,
                            'anvio_and_microbeannotator': 0,
                            'kofamscan_and_microbeannotator': 0,
                            'only_anvio': 0,
                            'only_microbeannotator': 0,
                            'only_kofamscan': 0,
                            }

for g in genomes_list:
    ks_file = ks_functions_folder + g + ".tsv"
    an_file = an_functions_folder + g + ".tsv"
    ma_file = ma_functions_folder + f"/{g}/annotation_results/{g}.faa.annot"
    ks_df = convert_kofamscan(ks_file)
    an_df = convert_anvio(an_file)
    ma_df = convert_microbeannotator(ma_file)

    for k in ko_count_dict:
        ks_gcids = set(ks_df[ks_df.KO == k]['gene'].to_list())
        an_gcids = set(an_df[an_df.KO == k]['gene'].to_list())
        ma_gcids = set(ma_df[ma_df.KO == k]['gene'].to_list())

        total = len(an_gcids.union(ks_gcids).union(ma_gcids))
        all_count = len(an_gcids.intersection(ks_gcids).intersection(ma_gcids))
        an_ks_count = len((an_gcids.intersection(ks_gcids)).difference(ma_gcids))
        an_ma_count = len((an_gcids.intersection(ma_gcids)).difference(ks_gcids))
        ma_ks_count = len((ma_gcids.intersection(ks_gcids)).difference(an_gcids))
        only_an_count = len((an_gcids.difference(ks_gcids)).difference(ma_gcids))
        only_ma_count = len((ma_gcids.difference(ks_gcids)).difference(an_gcids))
        only_ks_count = len((ks_gcids.difference(an_gcids)).difference(ma_gcids))
        
        ko_count_dict[k]['total_genes_annotated'] += total 
        ko_count_dict[k]['all_tools'] += all_count
        ko_count_dict[k]['anvio_and_kofamscan'] += an_ks_count
        ko_count_dict[k]['anvio_and_microbeannotator'] += an_ma_count
        ko_count_dict[k]['kofamscan_and_microbeannotator'] += ma_ks_count
        ko_count_dict[k]['only_anvio'] += only_an_count
        ko_count_dict[k]['only_microbeannotator'] += only_ma_count
        ko_count_dict[k]['only_kofamscan'] += only_ks_count
    print(f"Finished processing {g}")

count_table = pd.DataFrame().from_dict(ko_count_dict).T 
count_table.index.rename('KO', inplace=True)
count_table.to_csv(OUTPUT, sep="\t")

