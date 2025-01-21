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
        workdir, input_gene, structureid, 
        screen_name = '', 
        columns=[f'SUM_LFC3D_neg_psig', f'SUM_LFC3D_pos_psig'], 
        names=['negative', 'positive'], 
        max_distances=20, pthr_cutoff='p<0.05', score_type='LFC3D', 

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
    df_hits_clust[["x_coord", "y_coord", "z_coord"]] = df_struc

    # CLUSTERING #
    arr_d_thr = [float(i+1) for i in range(max_distances)]
    column_lists = [[] for _ in columns]
    hits_dict = dict(zip(columns, column_lists))
    
    for name, arr in hits_dict.items(): 
        # EXTRACT ROWS ABOVE CUTOFF #
        dict_hits = {}
        df_pvals_temp = df_hits_clust.loc[(df_hits_clust[name] == pthr_cutoff), ].reset_index(drop=True)
        dict_hits['unipos'] = list(df_pvals_temp['unipos'])

        # EXTRACT X Y Z OF HITS ABOVE CUTOFF #
        np_META_hits_coord = np.array(df_pvals_temp[['x_coord', 'y_coord', 'z_coord']].copy())
        if np_META_hits_coord.shape[0] < 2: 
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            return None, None

        # FOR RANGE OF RADIUS, RUN CLUSTERING #
        for thr_distance in arr_d_thr: 
            func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=thr_distance)
            clus_lbl = func_clustering.fit(np_META_hits_coord).labels_

            n_c_output = int(max(clus_lbl)+1)
            print(f'Number of clusters of {name} hits: @ d = {thr_distance} {n_c_output}')
            arr.append(n_c_output)

            colname_hits = f"{name}_hit_clust_{str(int(thr_distance))}"
            # df_pvals_temp[colname_hits] = clus_lbl
            dict_hits[colname_hits] = clus_lbl
        df_hits_clust = df_hits_clust.merge(pd.DataFrame(dict_hits), how='left', on=['unipos'])

    df_hits_clust.fillna('-')
    hits_clust_filename = edits_filedir / f"cluster_{score_type}/{structureid}_{screen_name}_Aggr_Hits.tsv"
    df_hits_clust.to_csv(hits_clust_filename, sep='\t', index=False)

    # PLOT #
    clust_dist_filename = edits_filedir / f"cluster_{score_type}/{structureid}_{screen_name}_Aggr_Hits_List.tsv" 
    yvals = list(hits_dict.values())
    plot_cluster_distance(arr_d_thr, yvals, names, clust_dist_filename, edits_filedir, input_gene, screen_name, score_type, 
                          subplots_kwargs)

    return arr_d_thr, yvals

def plot_cluster_distance(
        x, ys, names, out_filename, edits_filedir, input_gene, screen_name, score_type, 
        subplots_kwargs={'figsize':(10,7)}, 
): 
    dist_dict = {'clust_dist': x}
    for n, y in zip(names, ys): 
        dist_dict[n] = y
    dist_stat = pd.DataFrame(dist_dict)
    fig, ax = plt.subplots(**subplots_kwargs)

    for n in names: 
        sns.lineplot(data=dist_stat, x="clust_dist", y=n)

    plt.xlabel('Cluster Radius')
    plt.ylabel('Number of Clusters')
    plt.title(f'Positive vs Negative Clusters {input_gene}')
    plt.savefig(edits_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_cluster_distance.png") 

    dist_stat.to_csv(out_filename, sep = '\t', index=False)


def clustering_distance(
        df_struc, df_pvals, 
        thr_distance, 
        workdir, input_gene, input_uniprot, structureid, hits_clust_filename, 
        screen_name='', 
        columns=[f'SUM_LFC3D_neg_psig', f'SUM_LFC3D_pos_psig'], 
        names=['negative', 'positive'], 
        pthr_cutoff='p<0.05', score_type='LFC3D', 
        clustering_kwargs = {"n_clusters": None, "metric": "euclidean", "linkage": "single", }, 
        subplots_kwargs={'figsize':(20,15)}, 
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
    df_hits_clust[["x_coord", "y_coord", "z_coord"]] = df_struc

    colnames = dict(zip(names, columns))
    # OPEN CLUSTERING FILE #
    df_pvals_clust = pd.read_csv(edits_filedir / hits_clust_filename, sep = '\t')

    for name, colname in colnames.items(): 
        df_pvals_temp = df_hits_clust.loc[(df_pvals[colname] == pthr_cutoff), ].reset_index(drop=True)
        
        np_META_hits_coord = np.array(df_pvals_temp[["x_coord", "y_coord", "z_coord"]]).copy()
        if np_META_hits_coord.shape[0] < 2: 
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            raise None

        func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=thr_distance)

        clustering = func_clustering.fit(np_META_hits_coord)
        clus_lbl = clustering.labels_
        n_c_output = int(max(clus_lbl)+1)
        print(f'Number of clusters of {name} hits: {n_c_output}')

        dendogram_filename = edits_filedir / f"plots/{input_gene}_{input_uniprot}_{name}_Dendogram_{str(int(thr_distance))}A.png"
        fig = plot_dendrogram(clustering, df_pvals_temp, dendogram_filename, name, input_gene, subplots_kwargs)

        # CLUSTER INDEX AND LENGTH
        df_pvals_clust = df_pvals_clust.loc[(df_pvals_clust[colname] == pthr_cutoff), ].reset_index(drop=True)
        clust_indices = df_pvals_clust[f'{colname}_hit_clust_{str(int(thr_distance))}'].unique()

        for c in clust_indices: 
            this_c_data = df_pvals_clust.loc[df_pvals_clust[f'{colname}_hit_clust_{str(int(thr_distance))}'] == c, ].reset_index(drop=True)
            print(f'{c} : {len(this_c_data)} : ', this_c_data.at[0, 'unipos'], '-', this_c_data.at[len(this_c_data)-1, 'unipos'])

    return clust_indices

def plot_dendrogram(
        clustering, df, dendogram_filename, name, input_gene, 
        subplots_kwargs={'figsize':(20,15)}, 
):
    
    fig, ax = plt.subplots(**subplots_kwargs)
    # create the counts of samples under each node
    counts = np.zeros(clustering.children_.shape[0])
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
    xlbl = list(df['unipos'])

    # PLOT CORRESPONDING DENDROGRAM #
    dendrogram(linkage_matrix, color_threshold=6.0, labels=xlbl, leaf_rotation=90.)
    plt.title(f'{input_gene} {name} Clusters')
    plt.savefig(dendogram_filename, dpi=300)
    plt.show()
