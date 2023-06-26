<!-- remove all comments before releasing -->
<!-- This is the name of the module as it will appear in GenePatter, and its version, for clarity -->
# scGSEA

<!-- A brief text description of the module, usually one sentence in length. -->
**Description**: scGSEA is an extension of ssGSEA designed to improve the assessment of pathway activity in single-cell data by addressing sparsity and ensuring stable enrichment scoring.

<!-- This field is for the author or creator of the module. If the algorithm of the module is from a published paper, this is usually the first or corresponding author from the paper. If the module algorithm is unpublished, this is usually the developer of the module itself. This field can simply be a name of a person or group. -->
**Authors**: John Jun; UCSD - Mesirov Lab, UCSD

<!--This field is used for responding to help requests for the module, and should be an email address or a link to a website with contact information or a help forum. -->
**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help)

<!-- All modules have a version number associated with them (the last number on the LSID) that is used to differentiate between modules of the same name for reproducibility purposes. However, for publicly released software packages that are wrapped as GenePattern modules, sometimes this version number will be different that the version number of the algorithm itself (e.g. TopHat v7 in GenePattern uses version 2.0.8b of the TopHat algorithm). Since this information is often important to the user, the algorithm version field is an optional attribute that can be used to specify this different version number. Remove this field if not applicable -->
<!-- **Algorithm Version**: _OPTIONAL_ and Not applicable for this particular module -->

<!-- Why use this module? What does it do? If this is one of a set of modules, how does this module fit in the set? How does it work? write overview as if you are explaining to a novice. Include any links or images which would serve to clarify -->
## Summary

scGSEA is an extension of ssGSEA tailored for single-cell data analysis. It addresses the challenges of sparsity and unreliable enrichment scoring by employing specialized normalization methods and scoring metrics. By utilizing scGSEA, scientists can explore and interpret pathway activity and functional alterations within heterogeneous populations of cells, thereby advancing our understanding of complex biological systems.

## Discussion
scGSEA (single-cell Gene Set Enrichment Analysis) is an extension of ssGSEA (single-sample Gene Set Enrichment Analysis) specifically designed for analyzing single-cell data. While ssGSEA is commonly used for bulk gene expression data, scGSEA addresses the challenges posed by the sparsity and unreliability of enrichment scoring in single-cell datasets.

ssGSEA is a computational method used to assess the activity of predefined gene sets within individual samples. It quantifies the enrichment of gene sets by calculating an enrichment score for each sample based on the expression levels of genes within the gene sets. This approach allows researchers to determine the relative activity of biological pathways or sets of genes in a given sample, enabling the identification of biologically relevant signatures and functional alterations.

However, when applied to single-cell data, ssGSEA encounters certain challenges. Single-cell datasets are sparse, meaning that only a small fraction of genes are expressed in each individual cell. This sparsity can lead to unreliable enrichment scores and hinder the accurate interpretation of pathway activity. To overcome these limitations, scGSEA introduces modifications to the normalization and scoring metrics used in ssGSEA to ensure stability in the computed gene set enrichment scores.

In scGSEA, normalization methods are adapted to account for the sparsity of single-cell data. Traditional normalization techniques, such as log-transformations or scaling by total read counts, may not be suitable for single-cell data due to the presence of many zeros. Instead, normalization method that has been robust to sparsity is employed. Additionally, the scoring metric in scGSEA is modified to account for the variability and uncertainty associated with single-cell data.

By employing scGSEA, researchers can address various scientific research questions in the context of single-cell data analysis. For example, they can investigate the activity of specific pathways or gene sets across different cell types or conditions within a heterogeneous population of cells. scGSEA can also be used to identify genes or pathways that are differentially regulated between different cell clusters or states, aiding in the discovery of cellular heterogeneity and functional diversity.

Furthermore, scGSEA enables the integration of single-cell data with prior knowledge of biological pathways or gene sets, providing insights into the underlying mechanisms driving cellular processes. By comparing enrichment scores between different samples or conditions, researchers can gain a better understanding of how gene expression patterns relate to specific biological functions, disease phenotypes, or treatment responses.

<!-- appropriate papers should be cited here -->
## References
1. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102

2. Barbie, D. A., Tamayo, P., Boehm, J. S., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462:108-112. http://doi.org/10.1038/nature08460

## Source Links
* [The GenePattern scGSEA module source repository](https://github.com/genepattern/scGSEA/)
<!-- * scGSEA uses the [genepattern/scgsea Docker image](https://hub.docker.com/layers/150060459/genepattern/example-module/2/images/sha256-ae4fffff67672e46b251f954ad226b7ad99403c456c1c19911b6ac82f1a27f2f?context=explore)
* The Dockerfile used to build that image is [here.](https://github.com/genepattern/scgsea/docs/Dockerfile) -->

## Parameters
<!-- short description of the module parameters and their default values, as well as whether they are required -->

| Name | Description <!--short description--> | Default Value | 
---------|--------------|----------------|
| input_file * |  File containing raw counts data to be read in |
| chip_file  | Chip file used for conversion to gene symbols |
| gene_set_database_file *  | Gene set data in GMT format |
| output_file_name * | Basename to use for output file | scGSEA_scores |
| cluster_data_label | Metadata label for cell grouping information | seurat_clusters |
| cluster_data_file | Metadata file for cell grouping information |

*  Required


When supplying the cell grouping information,
* For `h5seurat`, `h5ad`, or `Seurat RDS` files : use the `cluster_data_label` parameter.
* For `10x MEX`, `h5` files : use the `cluster_data_file` parameter.


## Input Files
<!-- longer descriptions of the module input files. Include information about format and/or preprocessing...etc -->

1. `input_file`  
    This is a file containing raw counts gene expression data that is not normalized. The scGSEA module supports multiple input file formats including Seurat RDS, H5seurat, H5ad formats as well as 10x Market Exchange (MEX) and hdf5 (h5) formats. For a Seurat object, `$RNA@counts` slot will be accessed. For an AnnData object, `raw.X` slot will be accessed.
2. `chip_file`  
    This parameter’s drop-down allows you to select CHIP files from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current version of MSigDB. You can also upload your own gene set file(s) in [CHIP](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CHIP:_Chip_file_format_.28.2A.chip.29) format.
3. `gene_set_database_file`
    * This parameter’s drop-down allows you to select gene sets from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current version of MSigDB. You can also upload your own gene set file(s) in [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format.
    * If you want to use files from an earlier version of MSigDB you will need to download them from the archived releases on the [website](https://www.gsea-msigdb.org/gsea/downloads.jsp).
4. `output_file_name`  
    The prefix used for the name of the output GCT and CSV file. The default output prefix is `scGSEA_scores`. The output CSV and GCT files will contain the projection of input dataset onto a space of gene set enrichments scores.
5. `cluster_data_label`  
    The name of the metadata label for cluster information within the input Seurat/AnnData object. This label will be used to access the cell grouping information utilized for aggregating cells. The default value for this parameter is `seurat_clusters`, which is the metadata label for cluster annotations generated upon running Seurat.Clustering module. Use the default value when using the RDS file generated from the [Seurat.Clustering](https://github.com/genepattern/Seurat.Clustering) module. Otherwise, provide the appropriate metadata label for the slot that stores cell grouping information.
6. `cluster_data_file`  
    If your input file is `h5` or `10x MEX` format, a separate cluster data file (tab-delimited .txt file) must be supplied here. The grouping information in this file is used to aggregate cells prior to computing scGSEA scores. Therefore, if you have `h5` or `10x MEX` formatted files and do not have cluster data file, please use [ScanpyUtilities](https://github.com/genepattern/ScanpyUtilities) module or the GenePattern Seurat Suite ([Seurat.QC](https://github.com/genepattern/Seurat.QC) > [Seurat.Preprocess](https://github.com/genepattern/Seurat.Preprocess) > [Seurat.Clustering](https://github.com/genepattern/Seurat.Clustering)) to compute clusters.
   
    <img width="350" alt="Screenshot 2023-06-20 at 10 46 48 AM" src="https://github.com/genepattern/scGSEA/assets/111310290/9d9f624b-8d50-46ba-90ba-4ad3a8279150">

## Output Files
<!-- list and describe any files output by the module -->

1. `<output_file_name>.csv`   
    This is a gene set by cell cluster data consisted of scGSEA scores. 
2. `<output_file_name>.gct`   
    This is a gene set by cell cluster data consisted of scGSEA scores. The [HeatmapViewer module](https://github.com/genepattern/HeatMapViewer) can accept this file as input for generating heatmap visualizations.
3. `cluster_expression.csv`   
    This is a gene by cell cluster data consisted of normalized gene expression level. 
4. `stdout.txt`  
    This is standard output from the script.

## Example Data
<!-- provide links to example data so that users can see what input & output should look like and so that they and we can use it to test -->
Input:  
[input_file](https://datasets.genepattern.org/data/test_data/scGSEA/local.clustered.rds)   
[gene_set_database_file](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/h.all.v2023.1.Hs.symbols.gmt)   
[chip_file](https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations/human/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip)

Output:  
[scgsea_scores.csv](https://github.com/genepattern/scGSEA/blob/develop/data/91737/scGSEA_scores.csv)   
[scgsea_scores.gct](https://github.com/genepattern/scGSEA/blob/develop/data/91737/scGSEA_scores.gct)

## Platform Dependencies
Task Type:
Projection

CPU Type:
any

Operating System:
any

Language:
R 4.1.0, Python 3.8.5

## Version Comments
<!--For each version of a module, provide a short comment about what was changed in the new version of a module. Version comments consist of 3 parts: a date, a version number, and a short description. The date should be the release date of that version of the module, and the version number should match the version of the module for which it corresponds to. The description can be short, but should be informative (e.g. "added support for log transformed data", or "fixed bug with out of memory exception"). When a user views the documentation, all version comments up to and including the current version will be displayed, and act as a short version history for the module. -->

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
| 1 | May 21, 2023 | Preproduction version |
