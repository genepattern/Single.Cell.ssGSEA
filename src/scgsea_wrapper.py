import subprocess
import argparse
import os

# parser = argparse.ArgumentParser()
# # ~~~~Module Required Arguments~~~~~ #
parser.add_argument("--rscript", 
                    type=str, 
                    help="R script",
                    default='False')

parser.add_argument("--pyscript", 
                    type=str, 
                    help="Python script",
                    default='False')

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

print("Proprocessing the RDS data")
print("===========================")
r_command = ['Rscript', args.rscript, "--input_file", args.input_file]
subprocess.run(r_command)

# python_command = ['python3', args.pyscript, "--gene_set_database_file", args.gene_set_database_file, \
#                  "--output_file_name", args.output_file_name, "--chip_file", args.chip_file]
# subprocess.run(python_command) 
