import argparse
import json
import sys

from sc_ssGSEA import read_gmts, run_ssgsea_parallel, Expression

print(sys.argv)
parser = argparse.ArgumentParser()

parser.add_argument(
	"--input_file",
	type=str,
	help="Input file",
	default='False'
)

parser.add_argument(
	"--gene_set_database_file",
	type=str,
	help="gene set database file(s)",
	default='False'
)

parser.add_argument(
	"--output_file_name",
	type=str,
	help="filename to use for output files",
	default='False'
)

parser.add_argument(
	"--n_threads",
	type=str,
	help="job CPU count",
	default = 1
)

parser.add_argument(
	"--cluster_data_label",
	type=str,
	help="Metadata label to use for aggregating cells",
	default="seurat_clusters"
)

parser.add_argument(
	"--chip_file",
	type=str,
	help="chip file",
	default = None,
	required = False
)

args = parser.parse_args()

## Load and parse expression + cell labels

expr = Expression.get_expression_object(
	args.input_file,
	args.cluster_data_label
)

expr.load()

## Load gene sets

gs, gs_desc = read_gmts(args.gene_set_database_file)

## Run single cell ssGSEA
sc_ssGSEA_scores = run_ssgsea_parallel(
	expr.metacells,
	gs,
	n_job = args.n_threads,
	file_path = None
)

## Save results

sc_ssGSEA_scores.to_csv(
	args.output_file_name+".tsv",
	sep = '\t'
)











