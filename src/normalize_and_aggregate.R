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
# Paramter for setupR
parser <- add_option(parser, c("--input_file"), help = "RDS file to load.")

# ====================================
#parameter for saving plot
parser <- add_option(parser, c("--output_file_name"),type='character',default='metacells', help = "Basename of the file to be saved.")
# ====================================

print('==========================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('==========================================================')
# Setting up the PDF file for the plots
pdf(file=paste(args$output_file_name,'.pdf',sep=''))

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

# Perform CLR normalization across cells
SeuratObj <- NormalizeData(SeuratObj, normalization.method = "CLR", margin = 2, verbose = TRUE)

# Aggregate cells by metacell annotation
aggregatedObj <- AggregateExpression(SeuratObj, group.by = "idents", verbose = TRUE)

# Save the metacell expression vector in csv format

dev.off() # Close the PDF file
