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
# PARAMETERS for plotting
parser <- add_option(parser, c("--genes"), type='character', default='MYC, CD8A', help = "List of genes to visualize. If you write multiple genes they must be separated by a comma and a space (for example: MS4A1, CD8A, CD14).")
parser <- add_option(parser, c("--group_plots"), type='character', default='Horizontally', help = "How to group plots that are associated with the same marker.")
# ====================================
#parameter for saving plot
parser <- add_option(parser, c("--output_file_name"),type='character',default='SeuratMarkers', help = "Basename of the file to be saved.")
# ====================================


print('==========================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('==========================================================')
# Setting up the PDF file for the plots
pdf(file=paste(args$output_file_name,'.pdf',sep=''))

# Processing how to merge plots
print('Processing how to merge plots...')
if (args$group_plots=='Horizontally'){
  n_cols=1
}else if (args$group_plots=='Vertically'){
  n_cols=2
}else if (args$group_plots=='No'){
  n_cols=0
}else{
  print('This should not have happened')
  n_cols=0
}
print('... done!')

# print('Loading Seurt, scatter, and dplyr...')
cat('Loading Seurat...')
suppressMessages(library(Seurat))
# suppressMessages(library(scater))
# suppressMessages(suppressWarnings(library(dplyr)))
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

for (gene in unlist(strsplit(args$genes, ', '))){
    tryCatch(
        expr = {
            
            cat('Plotting gene named ')
            cat(gene)
            if (n_cols){
              vin_plot <- VlnPlot(SeuratObj, features = gene)
              elend_plot <- FeaturePlot(SeuratObj, features = gene)
              print(CombinePlots(plots = list(vin_plot, elend_plot),ncol=n_cols))
            }else{
              print(VlnPlot(SeuratObj, features = gene))
              print(FeaturePlot(SeuratObj, features = gene))
            }
            print('... done!')
        },
        error = function(e){
            cat('The gene named ')
            cat(gene)
            cat(' was not present in your dataset.')
            print(e)
            
            par(mar = c(0,0,0,0))
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, paste("The gene named\n",gene,"\nwas not present in your dataset."), 
                 cex = 1.6, col = "black")
            par(mar = c(5, 4, 4, 2) + 0.1)
        }
    )    
}

dev.off() # Close the PDF file
