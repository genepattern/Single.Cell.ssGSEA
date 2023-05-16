#!/usr/bin/env python3

from ExampleModule_functions import *
import argparse
import humanfriendly
import pandas as pd
from timeit import default_timer as timer
beginning_of_time = timer()

# parser = argparse.ArgumentParser()
# ~~~~Module Required Arguments~~~~~ #
# parser.add_argument("--gene_set_database_files",
#                     type=str,
#                     help="What message to add (if any)",
#                     default='False')

# args = parser.parse_args()
# if args.verbose:
#     print("Ah! The old verbosaroo")

print("~~~~~~~~~~~~~~~~~~~~~~")
print("Using arguments:")
print(args)
print("Now getting work done.")
print("~~~~~~~~~~~~~~~~~~~~~~")

# Open the input file
print("About to read the metacell expression")
if os.path.exists("metacell_expression.csv"):
  metacell_exp = pd.read_csv("metacell_expression.csv", index_col = 0)
print(metacell_exp)

# # Open the output file
# out_filename = args.output_filename

# # Load the single cell transcriptomic dataset
# data = sc.read_h5ad('../data/17_liver.h5ad')

# # Check for cluster annotation metadata
# if 'cluster' not in data.obs.keys():
#   raise ValueError("cluster annotation metadata not found.")

# # Normalize the dataset using CLR across cells
# norm_data = do_normalize(data, axis = 0, method = "clr")

# # Subset the data by celltype
# for ct in norm_data.obs['celltype'].unique():
#   norm_data_by_celltype['{}'.format(ct)] = norm_data[norm_data.obs['cell_type'] == ct]
  
# # Aggregate the cells by metacell
# exp_norm_by_celltype_aggregated = {}
# for key in exp_norm_by_celltype.keys():
#     exp_norm_by_celltype_aggregated[key] = do_aggregate(exp_norm_by_celltype[key], method = 'mean', aggregate_all = True)

# # Perform ssGSEA
# start_all = time.time()

# for comb in keys:
#     ct = comb[:comb.find('_')].replace(" ", "_")
#     comb_w_underscore = comb.replace(" ", "_")

# #### ---------- GENERATE PERMUTATIONS ---------- ####
# #     rank = "counts" in comb
# #     print("Generating Permutation for {} {}".format(comb, "-- Rank Normalize" if rank else ""))
# #     perm_cell = permute_ties_parallel(norm_exp_aggregated[comb], 50, rank, n_job = 8)
# #     perm_cell = [s for subset in perm_cell for s in subset]

    
# # Compute the score using KS
# temp = run_ssgsea_parallel_with_permutation_permlist(
#   perm_cell, gmt, n_perm = 50, stat = 'ks', n_job = 32)
#   temp.to_csv("../results/{}_{}_50perms_computational_genesets_20expgenes.csv".format(comb_w_underscore.replace("/", "_"), metric), sep='\t', header=True, index=True)
        
#         print("--------------------------------------------------------------------------------")        
#         print("\t Completed ssGSEA for {} using {} - 50 permutations".format(comb, metric))
        
# end_all = time.time()

# print("--------------------------------------------------------------------------------")        
# print("\t Total Time spent: {} ".format(end_all - start_all))


# with open(out_filename, 'w') as g:
#   for line in f.readlines():
#     g.write(line)
#     g.write('\n')  # This could be an argument, wether or not to add a newline.
# f.close()

# end_of_time = timer()
# print("We are done! Wall time elapsed:", humanfriendly.format_timespan(end_of_time - beginning_of_time))
