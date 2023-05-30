################################################################################
#Parse Parameters
################################################################################
print('==========================================================')
print("Loading library: optparse")
library("optparse")

# Parse input arguments
parser = OptionParser()
# parameter types: 'character', 'integer', 'logical', 'double', or 'complex'
# ====================================
# Paramter for the RDS input file
parser <- add_option(parser, c("--input_file"), help = "RDS file to load.")
parser <- add_option(parser, c("--cluster_data_label"), help = "Metadata label to use for aggregating cells")

print('==========================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('==========================================================')

cat('Loading Seurat...')
suppressMessages(library(Seurat))
print('...done loading libraries!')

################################################################################
#Begin Running the functions
################################################################################

print('==========================================================')
cat("About to read the Seurat object named")
print(args$input_file)
SeuratObj = readRDS(args$input_file)
cat(args$input_file)
print("has been read to memory!")

# Normalize the raw data
print("Normalizing the data...")
print('==========================================================')
SeuratObj <- NormalizeData(SeuratObj, normalization.method = "CLR", margin = 2, verbose = TRUE)
print("Normalizing complete!")

# Group by metacell 
print("Aggregating cells by clusters...")
print('==========================================================')
aggregatedObj <- AverageExpression(SeuratObj, group.by = args$cluster_data_label, use.raw = FALSE)
print("Aggregating cells by clusters complete!")

print("Saving metacell average expression profile...")
print('==========================================================')
write.csv(aggregatedObj, file = "cluster_expression.csv", row.names = TRUE)
print("cluster_expression.csv file saved!")
