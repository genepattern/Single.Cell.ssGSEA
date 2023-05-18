#!/usr/bin/env python3

from scgsea_helper import *
import argparse
import humanfriendly
import pandas as pd
from timeit import default_timer as timer
beginning_of_time = timer()

parser = argparse.ArgumentParser()
# ~~~~Module Required Arguments~~~~~ #
parser.add_argument("--input_file",
                    type=str,
                    help="Input file",
                    default='False')

parser.add_argument("--gene_set_database_file",
                    type=str,
                    help="gene set",
                    default='False')

parser.add_argument("--output_file_name",
                    type=str,
                    help="scGSEA scores",
                    default='False')

parser.add_argument("--chip_file",
                    type=str,
                    help="chip file",
                    default='False')

args = parser.parse_args()
# if args.verbose:
#     print("Ah! The old verbosaroo")

print("~~~~~~~~~~~~~~~~~~~~~~")
print("Using arguments:")
print(args)
print("Now getting work done.")
print("~~~~~~~~~~~~~~~~~~~~~~")

# Open the input file
# print("About to read the metacell expression")
# if os.path.exists("metacell_expression.csv"):
#   metacell_exp = pd.read_csv("metacell_expression.csv", index_col = 0)
# print(metacell_exp)

# # Open the output file
# out_filename = args.output_filename

# Load the metacell expression data
metacell_exp = pd.read_csv(args.input_file, index_col = 0)

# Load the chip file and convert to gene symbol
chip = read_chip(args.chip_file)
metacell_exp = convert_to_gene_symbol(chip, metacell_exp)
print("converted to gene symbol")

# Load the gene set database files
gs = pd.read_csv(args.gene_set_database_file, index_col = 0, header = None)

print("loaded both metacell and gene set")

scGSEA_scores = single_sample_gseas(
    metacell_exp,
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
