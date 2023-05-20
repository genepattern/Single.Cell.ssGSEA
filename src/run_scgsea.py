#!/usr/bin/env python3

from scgsea_helper import *
import argparse
import pandas as pd
import os

import humanfriendly
from timeit import default_timer as timer
beginning_of_time = timer()

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
                    help="Chip file",
                    default='False')

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
print("Loading CHIP file to convert to Gene Symbol")
print('==========================================================')
chip = read_chip(args.chip_file)
cluster_exp = convert_to_gene_symbol(chip, cluster_exp)
print("Loaded CHIP file!\n")

# Load the gene set database files
print(f"Loading {args.gene_set_database_file} to convert to Gene Symbol")
print('==========================================================')
# if args.gene_set_database_file.endswith(".txt"):
gs = read_gmts(args.gene_set_database_file)
print("Loaded gene set file!\n")

print("Running scGSEA...")
scGSEA_scores = single_sample_gseas(
    cluster_exp,
    gs,
    plot=True,
    title=None,
    gene_score_name=None,
    annotation_text_font_size=16,
    annotation_text_width=88,
    annotation_text_yshift=64,
    html_file_path="{}_plot".format(args.output_file_name),
    plotly_html_file_path=None)

write_gct(scGSEA_scores, args.output_file_name)

end_of_time = timer()
print("We are done! Wall time elapsed:", humanfriendly.format_timespan(end_of_time - beginning_of_time))
