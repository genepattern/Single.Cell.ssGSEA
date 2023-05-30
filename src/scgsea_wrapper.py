import subprocess
import argparse
import os

import humanfriendly
from timeit import default_timer as timer

start = timer()

parser = argparse.ArgumentParser()

# ~~~~Module Required Scripts~~~~~ #
parser.add_argument("--rscript", 
                    type=str, 
                    help="R script",
                    default='False')

parser.add_argument("--pyscript", 
                    type=str, 
                    help="Python script",
                    default='False')

parser.add_argument("--helper_functions",
                    type=str, 
                    help="Python helper functions",
                    default='False')

# ~~~~Module Required Arguments~~~~~ #
parser.add_argument("--input_file",
                    type=str,
                    help="Input file",
                    default='False')

parser.add_argument("--gene_set_database_file",
                    type=str,
                    help="gene set database file(s)",
                    default='False')

parser.add_argument("--output_file_name",
                    type=str,
                    help="filename to use for output files",
                    default='False')

parser.add_argument("--n_threads",
                    type=str,
                    help="job CPU count",
                    default = None)

parser.add_argument("--cluster_data_label",
                    type=str,
                    help="Metadata label to use for aggregating cells",
                    default="seurat_clusters")

# Optional parameter
parser.add_argument("--chip_file",
                    type=str,
                    help="chip file",
                    default = None)

args = parser.parse_args()

print('==========================================================')
print("Proprocessing the RDS data...")
print('==========================================================')

r_command = ['Rscript', args.rscript, "--input_file", args.input_file, "--cluster_data_label", args.cluster_data_label]
r_command_str = " ".join(r_command)
print(f"Running command: {r_command_str}")
subprocess.run(r_command)
print("Finished Preprocessing!\n")

print('==========================================================')
print("Performing scGSEA...")
print('==========================================================')

python_command = ['python3', args.pyscript, "--gene_set_database_file", args.gene_set_database_file, \
                 "--output_file_name", args.output_file_name, "--n_threads", args.n_threads]
if args.chip_file is not None:
    python_command.extend(["--chip_file", args.chip_file])
python_command_str = " ".join(python_command)
print(f"Running command: {python_command_str}") 
subprocess.run(python_command)
print("scGSEA Complete!\n")

end = timer()
print("We are done! Wall time elapsed:", humanfriendly.format_timespan(end- start))
