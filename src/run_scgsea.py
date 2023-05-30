#!/usr/bin/env python3

from scgsea_helper import *
import argparse
import pandas as pd
import os

import humanfriendly
from timeit import default_timer as timer

parser = argparse.ArgumentParser()
# ~~~~Module Required Arguments~~~~~ #
parser.add_argument("--gene_set_database_file",
                    type=str,
                    help="Gene set file",
                    default='False')

parser.add_argument("--output_file_name",
                    type=str,
                    help="Output file name",
                    default='False')

parser.add_argument("--chip_file",
                    type=str,
                    help="Chip file")

parser.add_argument("--n_threads",
                    type=int,
                    help="number of CPUs")

args = parser.parse_args()

print('==========================================================')
print("Running scGSEA for")
print(args.gene_set_database_file)

print("Now getting work done.")
print('==========================================================\n')

# Open the input file
print("About to read the metacell expression")
if os.path.exists("cluster_expression.csv"):
  cluster_exp = pd.read_csv("cluster_expression.csv", index_col = 0)
else:
  print("cluster_expression.csv not available")

# Load the chip file and convert to gene symbol
if args.chip_file:
  print("Loading CHIP file to convert to Gene Symbol")
  print('==========================================================')
  chip = read_chip(args.chip_file)
  cluster_exp = convert_to_gene_symbol(chip, cluster_exp)
  print("Loaded CHIP file!\n")
else:
  print("Chip file not found -- Using the original gene identifiers without conversion")

# Load the gene set database files
print(f"Loading {args.gene_set_database_file} to convert to Gene Symbol")
print('==========================================================')
# if args.gene_set_database_file.endswith(".txt"):
gs, gs_desc = read_gmts(args.gene_set_database_file)
print("Loaded gene set file!\n")

print("Number of gene sets loaded for scGSEA: {}".format(gs.shape[0]))

print("Running scGSEA...")
print('==========================================================')
start = timer()
scGSEA_scores = run_ssgsea_parallel(
    cluster_exp,
    gs,
    n_job = args.n_threads,
    file_path = None
)
end = timer()

orig_name = scGSEA_scores.columns[0]
if orig_name.startswith("RNA."):
    column_count = len(scGSEA_scores.columns)
    new_cols = ['cluster' + str(i) for i in range(1, column_count + 1)]
    scGSEA_scores.columns = new_cols

scGSEA_scores.to_csv(args.output_file_name + ".csv", sep="\t", mode = 'w')

write_gct(scGSEA_scores, args.output_file_name, gs_desc)

print("We are done!")

print(f"scGSEA Runtime using {args.n_threads} CPUs: ", humanfriendly.format_timespan(end - start))
