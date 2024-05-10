#!/usr/bin/python
# a script to combine various metabolism output files into one matrix
# used by the driver script generate_butyrate_completeness_matrix.sh 
# usage: python gen_butyrate_matrix.py GENOMES_FILE MODULE_OUTPUT_FOLDER

import os
import sys
import argparse
import pandas as pd

# COMMAND LINE PARAMETERS
if len(sys.argv) < 3:
    print("USAGE ERROR: not enough command-line arguments\n" +
          "python gen_butyrate_matrix.py GENOMES_FILE MODULE_OUTPUT_FOLDER")
    sys.exit(1)

genome_info_file = sys.argv[1]
module_output_folder = sys.argv[2]

# CONSTANTS
tool_list = ["anvio", "kofamscan", "MicrobeAnnotator"]
module_list = ["BUTANOATE", "SUBSPEC", "BUCASYNOP"]

genome_info = pd.read_csv(genome_info_file, sep="\t", index_col=0)
acc_list = genome_info.index.to_list()

for mod in module_list:
  completeness_matrix = pd.DataFrame(index = acc_list, columns = tool_list)
  for acc in acc_list:
    for tool in tool_list:
      infile = os.path.join(module_output_folder, tool, "default", acc + "_modules.txt")
      df = pd.read_csv(infile, sep="\t", index_col = 0)
      completeness_matrix.loc[acc, tool] = df.loc[mod,"pathwise_module_completeness"]

  completeness_matrix.index.name = "genome"
  completeness_matrix.to_csv(f"{mod}_completeness_matrix.txt", sep="\t")