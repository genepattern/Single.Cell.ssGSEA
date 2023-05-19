import subprocess
import argparse
import os

current_directory = os.getcwd()
print(current_directory)

file_list = os.listdir()
for file_name in file_list:
  print(file_name)

if os.path.exists('normalize_and_aggregate.R'):
  print("found it")
if os.path.exists('scgsea.py'):
  print("found it")
print("NOT FOUND")


# parser = argparse.ArgumentParser()
# # ~~~~Module Required Arguments~~~~~ #
# parser.add_argument("--input_file",
#                     type=str,
#                     help="Input file",
#                     default='False')

# parser.add_argument("--gene_set_database_file",
#                     type=str,
#                     help="gene set",
#                     default='False')

# parser.add_argument("--output_file_name",
#                     type=str,
#                     help="scGSEA scores",
#                     default='False')

# parser.add_argument("--chip_file",
#                     type=str,
#                     help="chip file",
#                     default='False')

# args = parser.parse_args()

# subprocess.run(['Rscript', ]
