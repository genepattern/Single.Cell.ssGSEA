# Single Cell ssGSEA

**Description**: Single Cell ssGSEA (sc_ssGSEA) is an extension of ssGSEA tailored for single-cell data analysis. It addresses the challenge of sparsity by employing a normalization method and scoring metric chosen to minimize any variability. By utilizing sc_ssGSEA, scientists can explore and interpret pathway activity and functional alterations within heterogeneous populations of cells.

**Authors**: John Jun, Alexander T. Wenzel; UCSD - Mesirov Lab, UCSD

**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help)

## Discussion
Single Cell ssGSEA is an extension of [ssGSEA](https://cloud.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00270:10.1.0)<sup>[2]</sup> (single-sample Gene Set Enrichment Analysis) specifically designed for analyzing single-cell data. sc_ssGSEA is a computational method used to assess the activation or repression of an *a priori* defined set of genes associated with a particular biological process or pathway in an individual sample (Barbie et al., Nature 2009). Genes are ranked by their absolute expression in the sample and assessed for their over-representation, i.e., enrichment, at the top or bottom of the list.

While ssGSEA was designed for use with bulk gene expression data, sc_ssGSEA addresses the challenge posed by the sparsity which may cause some variability in enrichment scoring due to the large number of genes which are not expressed or lowly expressed in single-cell datasets. To overcome these limitations, sc_ssGSEA aggregates (averages) cells within user-defined groups, such as clusters and introduces modifications to the normalization and scoring metrics used in ssGSEA to ensure stability in the computed gene set enrichment scores independent of the order that genes with the same mRNA abundance (or number of reads) appear in the ranked list. sc_ssGSEA utilizes the centered-log ratio normalization across cell prior to aggregation and uses the Kolmogorov-Smirnov (KS) statistic for the scoring metric. We find that this approach best reduces variability in the enrichment scores resulting from reordering of “tied” genes. 

## Source Links
* [The GenePattern scGSEA module source repository](https://github.com/genepattern/scGSEA/)
<!-- * scGSEA uses the [genepattern/scgsea Docker image](https://hub.docker.com/layers/150060459/genepattern/example-module/2/images/sha256-ae4fffff67672e46b251f954ad226b7ad99403c456c1c19911b6ac82f1a27f2f?context=explore)
* The Dockerfile used to build that image is [here.](https://github.com/genepattern/scgsea/docs/Dockerfile) -->

## Parameters
<!-- short description of the module parameters and their default values, as well as whether they are required -->

<table>
    <thead>
        <tr>
            <th>Parameter Group</th>
            <th>Name</th>
            <th>Description</th>
            <th>Default Value</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td colspan="1" rowspan="4" align="center">Input Files</td>
            <td>input file *</td> 
            <td>File containing raw counts or mRNA abundance estimates</td>
            <td></td>
        </tr>
        <tr>
            <td>gene set database file *</td>
            <td>Gene sets in GMT format</td>
            <td></td>
        </tr>
        <tr>
            <td>output file name *</td>
            <td>Basename to use for output file</td>
            <td><i>scGSEA_scores</i></td>
        </tr>
        <tr>
            <td>chip file</td>
            <td>Chip file used for conversion to gene symbols</td>
            <td></td>
        </tr>
        <tr>
            <td colspan="1" rowspan="2" align="center">Cell Grouping Data</td>
            <td>cluster data label</td> 
            <td>Metadata label for cell grouping (metacell) information; clustering data</td>
            <td><i>seurat_clusters</i></td>
        </tr>
    </tbody>
</table>

*  Required


### Input Files

1. `input file`  
    This is a file containing unnormalized gene expression data in raw read counts or estimated RNA abundance. The scGSEA module supports multiple input file formats including Seurat RDS, H5seurat, H5AD. For a Seurat object, the $RNA@counts slot will be used. For an AnnData object, the raw.X slot will be used. 

   * If you come across the following message in the `stderr.txt` file, please verify that the input file contains unnormalized raw counts data.

    &nbsp;<pre><code>The raw counts matrix was not composed of integer values. This may represent an issue with the processing pipeline. Please be advised...</code></pre>
   * If you have used `kallisto` or `salmon.alevin` for alignment, please disregard the message about the raw counts data not being in integer format; the aforementioned tools generate estimated RNA abundances, which may consist of non-integer count values.

2. `gene set database file`
    * This parameter’s drop-down allows you to select gene sets from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current (2023) version of MSigDB. You can also upload your own gene set file(s) in [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format.

    * If you want to use files from an earlier version of MSigDB you will need to download them from the archived releases on the [website](https://www.gsea-msigdb.org/gsea/downloads.jsp).

3. `chip file`  
    This parameter’s drop-down allows you to select CHIP files from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current version (2023) of MSigDB.

4. `output file name`  
    The prefix used for the name of the output GCT and CSV file. The default output prefix is <i>scGSEA_scores</i>. The output CSV and GCT files will contain a gene set x metacell matrix of enrichments scores.

5. `cluster data label`  
    The name of the metadata label for cell grouping information within the input Seurat/AnnData object. This label will be used to access the cell grouping information utilized for aggregating cells to create metacells. The default value for this parameter is <i>seurat_clusters</i>, which is the metadata label for the slot that stores cell-to-cluster mapping generated by the Genepattern Seurat.Clustering module. Use the default value when using the RDS file generated by the [Seurat.Clustering](https://github.com/genepattern/Seurat.Clustering) module. Otherwise, provide the appropriate metadata label for the slot that stores cell grouping information.

### Job Options
* `Job memory`  
    This is the memory allocated to your scGSEA job. The default memory is 32GB; however, appropriate adjustment is recommended based on the size of your input file. 
* `Job cpuCount`  
    scGSEA supports parallelization of the enrichment scoring process through multi-threading. The default value is 3 and the job will be divided into 3 subprocesses. It is recommended that you increase the **job cpuCount** when computing enrichment scores for a large number of gene sets.

## Output Files
<!-- list and describe any files output by the module -->

1. `<output_file_name>.tsv`   
    This is a gene set x metacell matrix consisted of scGSEA scores. 

## Example Data
<!-- provide links to example data so that users can see what input & output should look like and so that they and we can use it to test -->
Input:  
[input_file](https://datasets.genepattern.org/data/test_data/scGSEA/local.clustered.rds)   
[gene_set_database_file](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/h.all.v2023.1.Hs.symbols.gmt)   
[chip_file](https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations/human/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip)

Output:  
[scgsea_scores.csv](https://github.com/genepattern/scGSEA/blob/develop/data/91737/scGSEA_scores.csv)

## References
1. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102

2. Barbie, D. A., Tamayo, P., Boehm, J. S., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462:108-112. http://doi.org/10.1038/nature08460
