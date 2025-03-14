# Single Cell ssGSEA

Single Cell ssGSEA is an extension of the Single Sample Gene Set Enrichment Analysis (ssGSEA) method\[1\] for use with single cell RNA-sequencing (scRNA-seq) data. Because of the sparsity of scRNA-seq data, ssGSEA scores computed in individual cells are subject to uncertainty. The Single Cell ssGSEA approach reduces this uncertainty by:

- Aggregating cells according to user-defined cluster or cell type labels
- Using the weighted Kolmogorov-Smirnov score\[2\] instead of ssGSEA's area under the curve
- Normalizing counts using a centered log-ratio\[3\] prior to aggregation. 

For more details and benchmarking results, see this [preprint](https://www.biorxiv.org/content/10.1101/2024.06.03.597180v2).

# Installation

Single Cell ssGSEA is available as a [Docker image]([https://hub.docker.com/r/atwenzel/sc_ssgsea](https://hub.docker.com/r/genepattern/sc_ssgsea)) and as a [PyPI Python package](https://pypi.org/project/sc-ssGSEA/). This method is also available as a GenePattern module on the [GenePattern Cloud Server](https://cloud.genepattern.org/gp/pages/login.jsf).

To install the Docker image, run `docker pull genepattern/sc_ssgsea`. 

To install the Python package, run `python -m pip install sc-ssGSEA`. 

**Note**: Single Cell ssGSEA accepts multiple input formats (see below). To run Single Cell ssGSEA on a Seurat object saved in an `.rds` file, R and Seurat must be installed in the same environment. The Single Cell ssGSEA Docker image has both installed, but users of the PyPI package will need to install R and Seurat separately. 

# Using Single Cell ssGSEA

## GenePattern module

Create an account on the [GenePattern Cloud Server](https://cloud.genepattern.org/gp/pages/login.jsf) and search for the "Single Cell ssGSEA" module. 

**Parameters**

- `input.file` : A file containing a Seurat object in RDS format, an AnnData object in H5AD format, or a Seurat object in H5Seurat format. The file must end in `.rds`, `.h5ad`, or `.h5seurat` respectively. 

- `gene.sets` : The gene sets to test for enrichment. Choose one or more options from either the Human or Mouse collections listed in the dropdown. 

- `cluster.data.label` : The name of the column in the metadata contained in `input.file` which contains a grouping of cells. 

- `chip.file` : Optionally, provide a `.chip` file to transform genes into a different namespace or convert between human and mouse orthologous genes. 

- `output.file` : The name of the file that contains the results matrix. By default, this is `scores`. Input will have the suffix `.tsv` appended. 

**Please note that H5AD files must contain untransformed counts. Please see [Using H5AD files](#using-h5ad-files) for more information.**

## Python Package

The following code calls Single Cell ssGSEA on an RDS file containing a Seurat object, assuming that `seurat_object.rds` is the file, which contains a metadata column called `seurat_clusters`, and that gene sets are defined in `gene_sets.gmt`.

```python
from sc_ssGSEA import read_gmt, run_ssgsea_parallel, Expression

## Load and parse expression + cell labels

expr = Expression.get_expression_object(
    "seurat_object.rds",
    "seurat_clusters"
)

expr.load()

## Load gene sets

gs, _ = read_gmt("gene_sets.gmt")

## Run single cell ssGSEA

sc_ssGSEA_scores = run_ssgsea_parallel(
    expr.metacells,
    gs,
)
```

The first argument to `Expression.get_expression_object()` can be an RDS file containing a Seurat object, H5AD file containing an AnnData object, or an H5Seurat file containing a Seurat object. Single Cell ssGSEA relies on the file suffix, which must be either `.rds`, `.h5ad`, or `.h5seurat` respectively. 

**Please note that H5AD files must contain untransformed counts. Please see [Using H5AD files](#using-h5ad-files) for more information.**

### Advanced usage: Custom inputs

If you wish to use another file format, create a Python class that inherits from `Expression` and implements the `load()` method. The `load()` method may assume that it has access to the fields `_filepath: str` and `_group_name: str` (the metadata column), and it must populate the fields `_gene_names: List[str]`, `_cell_names: List[str]`, and `_group_labels: List[str]`. `load()` should then create a `scipy.sparse.csr_matrix` containing the expression data, call `Expression._normalize_sparse_matrix(sparse_mat: csr_matrix)` and then populate `_metacells: pandas.DataFrame` using the function`Expression._get_metacells(sparse_mat: csr_matrix)`. 

## Using H5AD files

While Scanpy workflows often involve transformations to the raw counts populating the `AnnData` object, Single Cell ssGSEA requires untransformed integer counts and will expect them in `adata.raw` field. Users should call `adata.raw = adata.copy()` prior to any normalization or other numerical transformations. Single Cell ssGSEA will check that `adata.raw` contains integer counts by sampling some expression values and checking that they are equivalent to the nearest integer value within a small tolerance threshold. 

# References

1. Barbie, D. A. et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature 462, 108–112 (2009).
2. Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences 102, 15545–15550 (2005).

