#!/usr/bin/env python3
#NB - all of these import statements should specify their versions and be executed in a separate script at Docker build time.

import pandas as pd
import numpy
import warnings
from numpy import absolute, in1d, sort


# ssGSEA code from PheNMF repository
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
    scgsea_scores = {
        gs_name: []
        for gs_name in gene_sets.index
    }
    
    for i, metacell in enumerate(gene_scores.columns):
        for gs_name in gene_sets.index:
            scgsea_scores[gs_name].append(
                single_sample_gsea(
                    gene_score = gene_scores[metacell],
                    gene_set_genes = gene_sets.loc[gs_name,:].dropna()
                )
            )
    scgsea_scores = pd.DataFrame(scgsea_scores).T
    scgsea_scores.columns = [gene_scores.columns]
    return scgsea_scores

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

def write_gct(out_matrix, filename):
    text_file = open(filename + ".gct", "w")
    text_file.write('#1.2\n')
    text_file.write(str(len(out_matrix)) + "\t" +
                        str(len(out_matrix.columns) - 1) + "\n")
    text_file.close()
    
    # Change the column names from RNA.<cluster_number> to cluster<cluster_number>
    column_count = len(out_matrix.columns)
    new_cols = ['cluster' + str(i) for i in range(1, column_count + 1)]
    out_matrix.columns = new_cols
    
    # Save as a gct file
    out_matrix.to_csv(filename + ".gct", sep="\t", mode='a')
    # Save as a csv file
    out_matrix.to_csv(filename + ".csv", sep="\t", mode='w')
    
def read_gmt(gs_db, thres_min=2, thres_max=2000):
    with open(gs_db) as f:
        temp=f.read().splitlines()
    max_Ng=len(temp)
    # temp_size_G will contain size of each gene set
    temp_size_G=list(range(max_Ng))
    for i in range(max_Ng):
        temp_size_G[i]=len(temp[i].split("\t")) - 2
    max_size_G=max(temp_size_G)
    gs=pd.DataFrame(numpy.nan, index=range(max_Ng), columns=range(max_size_G))
    temp_names=list(range(max_Ng))
    temp_desc=list(range(max_Ng))
    gs_count=0
    for i in range(max_Ng):
        gene_set_size=len(temp[i].split("\t")) - 2
        gs_line=temp[i].split("\t")
        gene_set_name=gs_line[0]
        gene_set_desc=gs_line[1]
        gene_set_tags=list(range(gene_set_size))
        for j in range(gene_set_size):
            gene_set_tags[j]=gs_line[j + 2]
        if numpy.logical_and(gene_set_size >= thres_min, gene_set_size <= thres_max):
            temp_size_G[gs_count]=gene_set_size
            gs.iloc[gs_count]=gene_set_tags + \
                list(numpy.full((max_size_G - temp_size_G[gs_count]), numpy.nan))
            temp_names[gs_count]=gene_set_name
            temp_desc[gs_count]=gene_set_desc
            gs_count=gs_count + 1
    Ng=gs_count
    gs_names=list(range(Ng))
    gs_desc=list(range(Ng))
    size_G=list(range(Ng))
    gs_names=temp_names[0:Ng]
    gs_desc=temp_desc[0:Ng]
    size_G=temp_size_G[0:Ng]
    gs.dropna(how='all', inplace=True)
    gs.index=gs_names
    return gs
#    return {'N_gs': Ng, 'gs': gs, 'gs_names': gs_names, 'gs_desc': gs_desc, 'size_G': size_G, 'max_N_gs': max_Ng}

def read_gmts(gs_dbs):
    gs = pd.DataFrame()
    for gs_db in gs_dbs:
        gs = pd.concat([gs, read_gmt(gs_db)], ignore_index=False)
    return gs
