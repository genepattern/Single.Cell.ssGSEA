import subprocess
import argparse
import os

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
                    help="Python script",
                    default='False')

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

# For when aggregating cells not using the annotation from Seurat.Clustering module
# parser.add_argument("--cluster_data",
#                     type=str,
#                     help="Name of the metadata to be used for aggregating cells",
#                     default="seurat_clusters")

args = parser.parse_args()

print('==========================================================')
print("Proprocessing the RDS data...")
print('==========================================================')
# r_command = ['Rscript', args.rscript, "--input_file", args.input_file, "--cluster_data", args.cluster_data]
r_command = ['Rscript', args.rscript, "--input_file", args.input_file]
subprocess.run(r_command)
print("Finished Preprocessing!\n")

print('==========================================================')
print("Performing scGSEA...")
print('==========================================================')
python_command = ['python3', args.pyscript, "--gene_set_database_file", args.gene_set_database_file, \
                 "--output_file_name", args.output_file_name, "--chip_file", args.chip_file]
subprocess.run(python_command) 
print("scGSEA Complete!\n")
