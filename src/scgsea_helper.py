#!/usr/bin/env python3
#NB - all of these import statements should specify their versions and be executed in a separate script at Docker build time.
# Here is where I'd put my functions, if I had any!
def permute_ties(s):
    new_cell = pd.Series([], dtype = 'float64')
    for uv in sort(s.unique())[::-1]:
        sub_s = s[s == uv]
        
        sub_s_shuffle = sub_s.sample(frac = 1, replace = False)
        
        new_cell = pd.concat([new_cell, sub_s_shuffle])      
    return new_cell

def permute_ties_ntimes(s, n, rank):
    new_cells = []
    for i in range(n):
        temp = permute_ties(s)
        if rank:
            temp = temp.rank(method = "dense")
        new_cells.append(temp)
    return new_cells

def permute_ties_parallel(s, n_perm, rank, n_job = 1):
    n_list = [int(n_perm / n_job) for i in range(n_job)]
    extra = n_perm - sum(n_list)
    if extra > 0:
        n_list[-1] += extra
    
    permutations = list(chain(
        multiprocess(
            # Function to perform multiple ssGSEAs 
            # Perm_cell_ - A subset of n_perm permutations
            # gene_sets - a list of gene sets for which to calculate the enrichment
            # stat - scoring metric
            # The result of each run of single_sample_gseas should be geneset x permutation(n_perm/n_job) matrix
            permute_ties_ntimes,
            ((s, n, rank) for n in n_list), n_job,
        ))
    )
    return permutations

def single_sample_gsea(
    gene_score,
    gene_set_genes,
    statistic="ks",
    plot=True,
    title=None,
    gene_score_name=None,
    annotation_text_font_size=16,
    annotation_text_width=88,
    annotation_text_yshift=64,
    html_file_path=None,
    plotly_html_file_path=None,
):

    gene_score = gene_score.dropna()
    ##############gene_score_sorted = gene_score.sort_values(ascending=False)
    gene_score_sorted = gene_score

    in_ = in1d(gene_score_sorted.index, gene_set_genes.dropna(), assume_unique=True)

    in_sum = in_.sum()

    if in_sum == 0:

        warn("Gene scores did not have any of the gene-set genes.")

        return

    gene_score_sorted_values = gene_score_sorted.values

    gene_score_sorted_values_absolute = absolute(gene_score_sorted_values)

    in_int = in_.astype(int)

    hit = (
        gene_score_sorted_values_absolute * in_int
    ) / gene_score_sorted_values_absolute[in_].sum()

    miss = (1 - in_int) / (in_.size - in_sum)

    y = hit - miss

    cumulative_sums = y.cumsum()

    if statistic not in ("ks", "auc"):

        raise ValueError("Unknown statistic: {}.".format(statistic))

    if statistic == "ks":

        max_ = cumulative_sums.max()

        min_ = cumulative_sums.min()

        if absolute(min_) < absolute(max_):

            score = max_

        else:

            score = min_

    elif statistic == "auc":

        score = cumulative_sums.sum()
        
    return score

def convert_adata_counts_to_gct(data):
    genexcell = pd.DataFrame(data.raw.X.toarray().T, index = data.var_names, columns = data.obs_names)
    return genexcell

def convert_adata_to_gct(data):
    df = data.to_df()
    exp, cell_names, gene_names = df.values, df.index, df.columns
    genexcell = pd.DataFrame(exp.T, index = gene_names, columns = cell_names)
    return genexcell

def CLR_normalize(x):
    # Select the positive values
    pos_x = x[x > 0]

    # Compute the geometric mean of the positive values
    gm_pos_x = np.exp(np.sum(np.log1p(pos_x)) / len(x))

    # Compute the centered log-ratio transform
    clr = np.log1p(x / gm_pos_x)

    return clr

def do_normalize(
    adata,
    axis = 1,
    method = "log"):
    df_copy = adata.copy()
    
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html
    if method == "cpm": # COUNTS PER MILLION -- SCANPY
        sc.pp.normalize_total(df_copy, target_sum = 1e6)
    
    # LogNormalize - As implemented in Seurat package
    # Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
    # This is then natural-log transformed using log1p.
    elif method == "lognorm":
        sc.pp.normalize_total(df_copy, target_sum=1e4)
        sc.pp.log1p(df_copy)
    
    # The expression matrix is log-transformed
    # The geometric mean is computed for each sample
    # Divide each count by the geometric mean
    # Log transform
    # In R: log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
    elif method == "clr":
        df_copy.X = np.apply_along_axis(CLR_normalize, axis, df_copy.X)

    # Normalize count data to relative counts per cell by dividing by the total per cell. 
    # Optionally use a scale factor. Default is 1.
    # Unnecessary as we have the CPM normalization method -- The only difference is the scale!
#     elif method == "rc":
#         sc.pp.normalize_total(df_copy, target_sum=1)
    
    # Log transform of the dataset
    elif method == "log":
        sc.pp.log1p(df_copy)
        
    else:
        raise ValueError("Requested normalization method unavailable. Please enter one of the following - (1)log (2)cpm (3)clr (4)rc (5)sfn (6)rn")
    return df_copy

def do_aggregate(
    gene_x_sample,
    n_cells = 0,
    method = "mean",
    aggregate_all = False):
    
    norm_exp = gene_x_sample.copy()
    if aggregate_all:
        n_cells = len(gene_x_sample.columns)

    if len(norm_exp.columns) < n_cells:
        raise ValueError("Error: number of cells to aggregate should be less than {}".format(len(norm_exp.columns)))
    else:
        norm_exp = norm_exp.loc[:, np.random.choice(norm_exp.columns, size=n_cells, replace=False)]
        if method == "mean":
            norm_exp = norm_exp.mean(axis = 1)
        elif method == "sum":
            norm_exp = norm_exp.sum(axis = 1)
        else:
            raise ValueError("Error: invalid method input - Use mean or sum")
    return norm_exp

# Run multiple ssGSEAs
def single_sample_gseas(
    gene_scores,
    gene_sets,
    perm_index,
    statistic="ks",
    plot=True,
    title=None,
    gene_score_name=None,
    annotation_text_font_size=16,
    annotation_text_width=88,
    annotation_text_yshift=64,
    html_file_path=None,
    plotly_html_file_path=None,
):
    perm_scores = {
        gs_name: []
        for gs_name in gene_sets.index
    }
    
    for i, gene_score in enumerate(gene_scores):
        for gs_name in gene_sets.index:
            perm_scores[gs_name].append(
                single_sample_gsea(
                    gene_score = gene_score,
                    gene_set_genes = gene_sets.loc[gs_name,:].dropna(),                    
                    statistic = statistic
                )
            )
    perm_scores = pd.DataFrame(perm_scores).T
    perm_scores.columns = ["perm {}".format(perm_index+i) for i in range(len(gene_scores))]
    return perm_scores

def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]
        
def run_ssgsea_parallel_with_permutation_permlist(
    perm_cell,
    gene_sets,
    n_perm,
    stat = 'ks',
    n_job = 1,
    file_path = None,
):

    np.random.seed(49)
    
    subarray_size = int(n_perm / n_job)
    
    # Score gene set x permutation
    print("Begin job delegation")
    score__gene_set_x_sample = pd.concat(
        multiprocess(
            # Function to perform multiple ssGSEAs 
            # Perm_cell_ - A subset of n_perm permutations
            # gene_sets - a list of gene sets for which to calculate the enrichment
            # stat - scoring metric
            # The result of each run of single_sample_gseas should be geneset x permutation(n_perm/n_job) matrix
            single_sample_gseas,
            (
                (perm_cell_, gene_sets, i*subarray_size, stat) for i, perm_cell_ in enumerate(list(divide_chunks(perm_cell, subarray_size)))
            ), n_job,
        ), sort = False, axis = 1
    )

    # ---------- Since different permutations are run, the column must be different ---------- #
    ## Assure columns come out in same order they came in
    #score__gene_set_x_sample = score__gene_set_x_sample[gene_x_sample.columns]

    if file_path is not None:
        score__gene_set_x_sample.to_csv(file_path, sep = '\t')

    return score__gene_set_x_sample

def multiprocess(callable_, args, n_job, random_seed=20121020):

    seed(random_seed)

    with Pool(n_job) as process:

        return process.starmap(callable_, args)
