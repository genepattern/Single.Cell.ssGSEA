#!/usr/bin/env python3
#NB - all of these import statements should specify their versions and be executed in a separate script at Docker build time.

import pandas as pd
from numpy import absolute, in1d, sort

# Here is where I'd put my functions, if I had any!
def single_sample_gsea(
    gene_score,
    gene_set_genes,
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
    gene_score_sorted = gene_score.sort_values(ascending=False)

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
    
    # KS scoring
    max_ = cumulative_sums.max()
    min_ = cumulative_sums.min()
    if absolute(min_) < absolute(max_):
        score = max_
    else:
        score = min_

    return score

# Run multiple ssGSEAs
def single_sample_gseas(
    gene_scores,
    gene_sets,
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
    
    for i, metacell in enumerate(gene_scores.columns):
        for gs_name in gene_sets.index:
            perm_scores[gs_name].append(
                single_sample_gsea(
                    gene_score = gene_scores[metacell],
                    gene_set_genes = gene_sets.loc[gs_name,:].dropna()
                )
            )
    perm_scores = pd.DataFrame(perm_scores).T
    perm_scores.columns = [gene_scores.columns]
    return perm_scores

def read_chip(chip):
    chip_df=pd.read_csv(chip, sep='\t', index_col=0, skip_blank_lines=True)
    return chip_df

def convert_to_gene_symbol(chip, exp):
    joined_df = chip.join(exp, how='inner')
    joined_df.reset_index(drop=True, inplace=True)
    annotations = joined_df[["Gene Symbol", "Gene Title"]].drop_duplicates().copy()
    joined_df.drop("Gene Title", axis = 1, inplace = True)
    collapsed_df = joined_df.groupby(["Gene Symbol"]).max()
    return collapsed_df

def write_gct(out_matrix, file_name):
    text_file = open(filename + ".gct", "w")
        text_file.write('#1.2\n')
        text_file.write(str(len(out_matrix)) + "\t" +
                        str(len(out_matrix.columns) - 1) + "\n")
        text_file.close()
        out_matrix.to_csv(filename + ".gct", sep="\t", mode='a')
