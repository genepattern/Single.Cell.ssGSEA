<!-- remove all comments before releasing -->
<!-- This is the name of the module as it will appear in GenePatter, and its version, for clarity -->
# scGSEA (v1)

<!-- A brief text description of the module, usually one sentence in length. -->
**Description**: Performs single sample GSEA for single cell data

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

2. Reich, M., Liefeld, T., Gould, J., Lerner, J., Tamayo, P., & Mesirov, J. P. (2006). GenePattern 2.0. Nature Genetics, 38(5), 500-501. https://doi.org/10.1038/ng0506-500
<!-- links to your source repository **specific to the release version**, the Docker image used by the module (as specified in your manifest), and (if applicable) the sha link to the Dockerfile used to build your Docker image -->

## Source Links
* [The GenePattern ExampleModule v2 source repository](https://github.com/genepattern/scGSEA/)
* scGSEA uses the [genepattern/scgsea:2 Docker image](https://hub.docker.com/layers/150060459/genepattern/example-module/2/images/sha256-ae4fffff67672e46b251f954ad226b7ad99403c456c1c19911b6ac82f1a27f2f?context=explore)
* The Dockerfile used to build that image is [here.](https://github.com/genepattern/scgsea/docs/Dockerfile)

## Parameters
<!-- short description of the module parameters and their default values, as well as whether they are required -->

| Name | Description <!--short description--> | Default Value |
---------|--------------|----------------
| input_file * |  File to be read in RDS format |
| chip_file  | Chip file used for conversion to gene symbols |
| gene_set_database_file *  | Gene set data in GMT format |
| output_file_name * | The basename to use for output file | scGSEA_scores

\*  required

## Input Files
<!-- longer descriptions of the module input files. Include information about format and/or preprocessing...etc -->

1. filename  
    A long form explanation of the parameter. For example: This is the file which will be read in by the python script and to which text will be added, if add_custom_message is set to true. The parameter expects a text file with a .txt extension (e.g. file.txt)
    
## Output Files
<!-- list and describe any files output by the module -->

1. \<output_filename\>.txt
    The input file plus any text you added, if you chose to add text.
2. stdout.txt
    This is standard output from the Python script. Sometimes helpful for debugging.

## Example Data
<!-- provide links to example data so that users can see what input & output should look like and so that they and we can use it to test -->

Input:  
[data_placeholder.txt](https://github.com/genepattern/ExampleModule/blob/v2/data/data_placeholder.txt)

Output:  
[created_file_ground_truth.txt](https://github.com/genepattern/ExampleModule/blob/v2/gpunit/output/basic_test/created_file_ground_truth.txt)

## Platform Dependencies
Task Type:
Projection

CPU Type:
any

Operating System:
any

Language:
R-3.2

## Version Comments
<!--For each version of a module, provide a short comment about what was changed in the new version of a module. Version comments consist of 3 parts: a date, a version number, and a short description. The date should be the release date of that version of the module, and the version number should match the version of the module for which it corresponds to. The description can be short, but should be informative (e.g. "added support for log transformed data", or "fixed bug with out of memory exception"). When a user views the documentation, all version comments up to and including the current version will be displayed, and act as a short version history for the module. -->

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
| 1 | May 09, 2023 | Initial version |
