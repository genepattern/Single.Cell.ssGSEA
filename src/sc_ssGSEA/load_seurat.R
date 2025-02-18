#' Script to be invoked by Expression.SeuratRDS_Expression.load()
#' Saves a sparse matrix, genes, barcodes, and metadata in a Python-compatible
#' format in `/tmp/`


library(Matrix)
library(Seurat)

args <- commandArgs(trailingOnly=TRUE)

rds.path <- args[1]

seurat.obj <- readRDS(rds.path)

## Save expression

Matrix::writeMM(
	seurat.obj@assays[["RNA"]]@counts,
	"/tmp/expr.mtx"
)

## Save genes

genes.file <- file("/tmp/genes.tsv")
writeLines(rownames(seurat.obj), genes.file)
close(genes.file)

## Save barcodes

barcodes.file <- file("/tmp/barcodes.tsv")
writeLines(colnames(seurat.obj), barcodes.file)
close(barcodes.file)

## Save metadata

write.table(
	seurat.obj@meta.data,
	"/tmp/metadata.tsv",
	sep = '\t',
	quote = FALSE,
)