import pandas as pd
import sys


original_data = pd.read_csv(sys.argv[1])
target_individuals = 5000

scaling_factor = (target_individuals / sum(original_data["n_complete"]))
original_data["downsampled_numbers"] = [round(scaling_factor * count) for count in original_data["n_complete"]]
original_data[["family","n_complete","downsampled_numbers"]].to_csv(sys.argv[2], index=False, header=None)
