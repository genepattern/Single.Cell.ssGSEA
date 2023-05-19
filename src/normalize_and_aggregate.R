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
# Paramter for the input file
parser <- add_option(parser, c("--input_file"), help = "RDS file to load.")

# ====================================
# parameter for the output file name
# parser <- add_option(parser, c("--gene_set_file"), help = "GMT file to load.")
# ====================================

# ====================================
# parameter for the output file name
# parser <- add_option(parser, c("--output_file_name"),type='character',default='scGSEA_scores', help = "Basename of the file to be saved.")
# ====================================

print('==========================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('==========================================================')
# Setting up the PDF file for the plots
# pdf(file=paste(args$output_file_name,'.pdf',sep=''))

# Processing how to merge plots

cat('Loading Seurat...')
suppressMessages(library(Seurat))
print('...done loading libraries!')

################################################################################
#Begin Running the functions
################################################################################

print('==========================================================')
# Call the setupR function

cat("About to read the Seurat object named")
print(args$input_file)
SeuratObj = readRDS(args$input_file)
cat(args$input_file)
print("has been read to memory!")
library(Seurat)

# Normalize the raw data
print("Normalizing the data")
print('==========================================================')
SeuratObj <- NormalizeData(SeuratObj, normalization.method = "CLR", margin = 2, verbose = TRUE)

# Group by metacell 
print("Aggregating the cells by clusters")
print('==========================================================')
aggregatedObj <- AverageExpression(SeuratObj, group.by = "seurat_clusters", use.raw = FALSE, verbose = TRUE)

print("Saving metacell average expression profile.")
print('==========================================================')
write.csv(aggregatedObj, file = "metacell_expression.csv", row.names = TRUE)

#dev.off() # Close the PDF file
