#!/usr/bin/python
# a script to combine various metabolism output files into one matrix
# used by the driver script generate_butyrate_completeness_matrix.sh 

import os;
import pandas as pd

genome_info_file = "Lachno_genomes.txt"
genome_info = pd.read_csv(genome_info_file, sep="\t")
acc_list = genome_info.accession_short.to_list()
tool_list = ["anvio", "kofamscan", "MicrobeAnnotator"]

completeness_matrix = pd.DataFrame(index = acc_list, columns = tool_list)
for acc in acc_list:
  for tool in tool_list:
    infile = os.path.join("Archive/BUTANOATE_MODULE_OUTPUT", tool, "default", acc + "_modules.txt")
    df = pd.read_csv(infile, sep="\t")
    completeness_matrix.loc[acc, tool] = df[df.module == "BUTANOATE"].loc[0,"pathwise_module_completeness"]

completeness_matrix.index.name = "genome"
completeness_matrix.to_csv("butanoate_completeness_matrix.txt", sep="\t")

# rename rows by species name rather than GTDB accession
acc_to_species = {acc: genome_info.loc[acc, "species"].replace("s__", "") + f"  ({acc})" for acc in genome_info.index}
completeness_matrix.rename(acc_to_species, inplace=True)
completeness_matrix.to_csv("butanoate_completeness_matrix_labeled.txt", sep="\t")