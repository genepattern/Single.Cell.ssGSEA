# scGSEA

This module performs scGSEA (single-cell GSEA) for single cell transcriptomics data. This module is intended to be used subsequent to Seurat.Clustering module and the user can supply the Seurat RDS file from the Seurat.Clustering module.

The scGSEA score can be used to identify the biological mechanism and in downstream visualzations.

## Parameters
- **Seurat RDS file** (required)
- **Gene Set Database file** (required)
  - MSigDB gene set collections are available. Users can supply independent gene set database file in GMT format.
- Chip file
  - Mappings for human, mouse, and rat gene identifiers
- Output file name

## Docker
- This module uses genepattern/scgsea:v0.1

## Contact
- John Jun -- UCSD Mesirov Lab
- johnjun094@cloud.ucsd.edu


