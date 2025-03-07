"""
File: annotate_spatial_clusters.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 5

"""

import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from pathlib import Path
import os
import warnings

def clustering(
        df_struc, df_pvals, 
        workdir, input_gene, screen_name = 'Meta', 
        columns=[f'SUM_LFC3D_neg_psig', f'SUM_LFC3D_pos_psig'], 
        names=['negative', 'positive'], 
        max_distances=20, pthr_cutoff=['p<0.05'], score_type='LFC3D', merge_col='unipos',

        clustering_kwargs = {"n_clusters": None, "metric": "euclidean", "linkage": "single", }, 
        subplots_kwargs={'figsize':(10,7)}, 
): 
    """
    Description
        Calculates number of clusters for a range of clustering radii

    Params
        df_struc: pandas dataframe, required
            DataFrame containing the x y z coordinates of the same length
        df_pvals: pandas dataframe
            from previous step average_split_bin_lfc3d() or average_split_bin_metaaggregation()
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        structureid: str, required
            the name of the AF and uniprot input
        screen_name: str, optional
            the name of the input screen
        max_distances: int, optional
            the upper limit of the range to test for clustering threshold
        pthr: float, optional
            the p-value threshold
        score_type: str, optional
            'LFC' or 'LFC3D'
        clustering_kwargs: dict, optional
            input params into AgglomerativeClustering

    Returns
    arr_d_thr: list of ints
        values for clustering radius
    yvals: list of lists of ints
        2 lists (neg and pos) where the values are number of clusters
    """

    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / f'cluster_{score_type}'):
        os.mkdir(edits_filedir / f'cluster_{score_type}')

    df_hits_clust = df_pvals.copy()
    if len(df_struc.columns) == 3: 
        df_hits_clust[["x_coord", "y_coord", "z_coord"]] = df_struc
    elif len(df_struc.columns) == 4: 
        df_hits_clust = pd.merge(df_hits_clust, df_struc, on=merge_col, how='left') 
        df_hits_clust = df_hits_clust.rename(columns=dict(zip(df_struc.columns, ['x_coord', 'y_coord', 'z_coord'])))
    else: 
        warnings.warn("df_struc must have 3 or 4 columns")

    # CLUSTERING #
    arr_d_thr = [float(i+1) for i in range(max_distances)] # CLUSTERING DISTANCE HYPERPARAM
    column_lists = [[] for _ in columns]
    columns_dict = dict(zip(columns, column_lists)) # FOR EVERY DATA COLUMN TO CLUSTER
    
    for name, arr in columns_dict.items(): 
        # EXTRACT ROWS ABOVE CUTOFF #
        dict_hits = {}
        df_pvals_temp = df_hits_clust.loc[(df_hits_clust[name].isin(pthr_cutoff)), ].reset_index(drop=True)
        # REMOVE ROWS WITHOUT POSITION INFO FOR PDBs #
        df_pvals_temp = df_pvals_temp[~df_pvals_temp[['x_coord', 'y_coord', 'z_coord']].isin(['-']).any(axis=1)]
        dict_hits[merge_col] = list(df_pvals_temp[merge_col])

        # EXTRACT X Y Z OF HITS ABOVE CUTOFF #
        np_META_hits_coord = np.array(df_pvals_temp[['x_coord', 'y_coord', 'z_coord']].copy())
        if np_META_hits_coord.shape[0] < 2: # NO DATA TO CLUSTER ON #
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            return None, None, None

        # FOR RANGE OF RADIUS, RUN CLUSTERING #
        for dist in arr_d_thr: 
            func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=dist)
            clus_lbl = func_clustering.fit(np_META_hits_coord).labels_

            n_c_output = int(max(clus_lbl)+1)
            print(f'Number of clusters of {name} hits: @ d = {dist} {n_c_output}')
            arr.append(n_c_output)

            dict_hits[f"{name}_Clust_{str(int(dist))}A"] = clus_lbl
        df_hits_clust = df_hits_clust.merge(pd.DataFrame(dict_hits), how='left', on=[merge_col])

    df_hits_clust.fillna('-')
    hits_filename = edits_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_Aggr_Hits.tsv"
    df_hits_clust.to_csv(hits_filename, sep='\t', index=False)

    # PLOT #
    yvals = list(columns_dict.values())
    plot_cluster_distance(arr_d_thr, yvals, names, edits_filedir, 
                          input_gene, screen_name, score_type, subplots_kwargs)

    return df_hits_clust, arr_d_thr, yvals

def plot_cluster_distance(
        x, ys, names, edits_filedir, input_gene, screen_name, score_type, 
        subplots_kwargs={'figsize':(10,7)}, 
): 
    dist_dict = {'clust_dist': x}
    for n, y in zip(names, ys): 
        dist_dict[n] = y
    dist_stat = pd.DataFrame(dist_dict)
    clust_filename = edits_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_Aggr_Hits_List.tsv" 
    dist_stat.to_csv(clust_filename, sep = '\t', index=False)

    fig, ax = plt.subplots(**subplots_kwargs)
    for n in names: 
        sns.lineplot(data=dist_stat, x="clust_dist", y=n)

    plt.xlabel('Cluster Radius')
    plt.ylabel('Number of Clusters')
    plt.title(f'Positive vs Negative Clusters {input_gene}')
    plot_filename = edits_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_cluster_distance.png"
    plt.savefig(plot_filename, dpi=300) 


def clustering_distance(
        df_struc, df_pvals, df_pvals_clust, 
        dist, 
        workdir, input_gene, screen_name = 'Meta', score_type='LFC3D',  
        columns=[f'SUM_LFC3D_neg_psig', f'SUM_LFC3D_pos_psig'], 
        names=['negative', 'positive'], 
        pthr_cutoff='p<0.05', merge_col='unipos', 
        clustering_kwargs = {"n_clusters": None, "metric": "euclidean", "linkage": "single", }, 
        subplots_kwargs={'figsize':(15, 12)}, horizontal=False, 
): 
    """
    Description
        Calculates number of clusters for one clustering radius

    Params
        df_struc: pandas dataframe, required
            DataFrame containing the x y z coordinates of the same length
        df_pvals: pandas dataframe
            from previous step average_split_bin_lfc3d() or average_split_bin_metaaggregation()
        thr_distance: int
            the chosen clustering radius
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        structureid: str, required
            the name of the AF and uniprot input
        screen_name: str, optional
            the name of the input screen
        pthr: float, optional
            the p-value threshold
        score_type: str, optional
            'LFC' or 'LFC3D'
        clustering_kwargs: dict, optional
            input params into AgglomerativeClustering

    Returns
    arr_d_thr: list of ints
        values for clustering radius
    yvals: list of lists of ints
        2 lists (neg and pos) where the values are number of clusters
    """
    
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    df_hits_clust = df_pvals.copy()
    if len(df_struc.columns) == 3: 
        df_hits_clust[["x_coord", "y_coord", "z_coord"]] = df_struc
    elif len(df_struc.columns) == 4: 
        df_hits_clust = pd.merge(df_hits_clust, df_struc, on=merge_col, how='left') 
        df_hits_clust = df_hits_clust.rename(columns=dict(zip(df_struc.columns, ['x_coord', 'y_coord', 'z_coord'])))
    else: 
        warnings.warn("df_struc must have 3 or 4 columns")

    columns_dict = dict(zip(names, columns))
    # OPEN CLUSTERING FILE #

    for name, colname in columns_dict.items(): 
        df_pvals_temp = df_hits_clust.loc[(df_pvals[colname].isin(pthr_cutoff)), ].reset_index(drop=True)
        df_pvals_temp = df_pvals_temp[~df_pvals_temp[['x_coord', 'y_coord', 'z_coord']].isin(['-']).any(axis=1)]

        np_META_hits_coord = np.array(df_pvals_temp[["x_coord", "y_coord", "z_coord"]]).copy()
        if np_META_hits_coord.shape[0] < 2: 
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            raise None

        func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=dist)
        clustering = func_clustering.fit(np_META_hits_coord)

        fig = plot_dendrogram(clustering, df_pvals_temp, edits_filedir, name, input_gene, 
                              screen_name, score_type, dist, merge_col, horizontal, subplots_kwargs)

        # CLUSTER INDEX AND LENGTH
        df_pvals_clust = df_pvals_clust.loc[(df_pvals_clust[colname].isin(pthr_cutoff)), ].reset_index(drop=True)
        clust_indices = df_pvals_clust[f'{colname}_Clust_{str(int(dist))}A'].unique()

        txt_filename = edits_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_{name}_Dendogram_{str(int(dist))}A.txt"
        with open(txt_filename, "w") as f:
            for c in clust_indices: 
                c_data = df_pvals_clust.loc[df_pvals_clust[f'{colname}_Clust_{str(int(dist))}A'] == c, ].reset_index(drop=True)
                if len(c_data) == 0: # IN THE CASE OF INCOMPLETE STRUCTURE DATA (PDB) SOME CLUSTERS ARE REMOVED #
                    continue
                all_unipos = c_data[merge_col].tolist()
                f.write(f'{c} : {len(c_data)} : {c_data.at[0, merge_col]} - {c_data.at[len(c_data)-1, merge_col]}\n')
                f.write(f'   {all_unipos}\n')

    return clust_indices

def plot_dendrogram(
        clustering, df_pvals_temp, 
        edits_filedir, name, input_gene, screen_name, 
        score_type, dist, merge_col, horizontal, 
        subplots_kwargs={'figsize':(20,15)}, 
):
    
    fig, ax = plt.subplots(**subplots_kwargs)
    counts = np.zeros(clustering.children_.shape[0]) # CREATE COUNTS OF SAMPLE
    n_samples = len(clustering.labels_)
    
    for i, merge in enumerate(clustering.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [clustering.children_, clustering.distances_, counts]).astype(float)
    xlbl = list(df_pvals_temp[merge_col])

    # PLOT CORRESPONDING DENDROGRAM #
    if horizontal: 
        dendrogram(linkage_matrix, color_threshold=6.0, labels=xlbl, orientation='right')
    else: 
        dendrogram(linkage_matrix, color_threshold=6.0, labels=xlbl, leaf_rotation=90.)
    plt.title(f'{input_gene} {score_type} {name} Clusters')
    dendogram_filename = edits_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_{score_type}_{name}_Dendogram_{str(int(dist))}A.png"
    plt.savefig(dendogram_filename, dpi=300)
    plt.show()
