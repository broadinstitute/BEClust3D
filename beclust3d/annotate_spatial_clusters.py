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
        df_struc_consvr, df_META, 
        workdir, 
        input_gene, structureid, 
        i_affn='euclidean', i_link='single', n_clusters=None, 
        max_distances=20, pthr=0.05, 
):

    edits_filedir = Path(workdir + '/' + input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'cluster_LFC3D'):
        os.mkdir(edits_filedir / 'cluster_LFC3D')

    df_hits_clust = df_META.copy()
    df_hits_clust["x_coord"] = df_struc_consvr["x_coord"]
    df_hits_clust["y_coord"] = df_struc_consvr["y_coord"]
    df_hits_clust["z_coord"] = df_struc_consvr["z_coord"]

    # CLUSTERING #
    arr_d_thr = [float(i+1) for i in range(max_distances)]
    hits = {'sensitizing': {'name':'SUM_LFC3D_neg_psig', 'arr':[]},
            'resistant':   {'name':'SUM_LFC3D_pos_psig', 'arr':[]}, }
    
    for key, val in hits.items(): 
        df_META_temp = df_hits_clust.loc[(df_hits_clust[val['name']] == 'p<'+str(pthr)), ]
        df_META_temp = df_META_temp.reset_index(drop=True)
        dict_hits = {}
        dict_hits['unipos'] = list(df_META_temp['unipos'])

        for thr_distance in arr_d_thr: # for every radius value 1.0 to 20.0
            colname_hits = f"{key}_hit_clust_{str(round(thr_distance))}"
            
            df_META_hits_coord = pd.DataFrame()
            df_META_hits_coord['x_coord'] = df_META_temp['x_coord']
            df_META_hits_coord['y_coord'] = df_META_temp['y_coord']
            df_META_hits_coord['z_coord'] = df_META_temp['z_coord']
            np_META_hits_coord = np.array(df_META_hits_coord)
            
            func_clustering = AgglomerativeClustering(n_clusters=n_clusters, metric=i_affn, 
                                                      linkage=i_link, distance_threshold=thr_distance)
            
            if np_META_hits_coord.shape[0] < 2: 
                warnings.warn(f"Not enough data to perform agglomerative clustering")
                return None, None
            clustering = func_clustering.fit(np_META_hits_coord)
            clus_lbl = clustering.labels_
            n_c_output = int(max(clus_lbl)+1)
            print(f'Number of clusters of {key} hits:', '@ d = ', thr_distance, n_c_output)
            val['arr'].append(n_c_output)
            df_META_temp[colname_hits] = clus_lbl
            dict_hits[colname_hits] = clus_lbl

        df_hits_clust = df_hits_clust.merge(pd.DataFrame(dict_hits), how='left', on=['unipos'])

    df_hits_clust.fillna('-')
    hits_clust_filename = edits_filedir / f"cluster_LFC3D/{structureid}_MetaAggr_Hits_Clust_p_l001.tsv"
    df_hits_clust.to_csv(hits_clust_filename, sep='\t', index=False)

    # PLOT #
    cluster_distance_filename = edits_filedir / f"cluster_LFC3D/{structureid}_MetaAggr_Hits_Clust_Dist_Stat_pl001.tsv"
    plot_cluster_distance(arr_d_thr, [hits['sensitizing']['arr'], hits['resistant']['arr']], 
                          cluster_distance_filename, edits_filedir)

    return arr_d_thr, (hits['sensitizing']['arr'], hits['resistant']['arr'])


def plot_cluster_distance(
        x, ys, out_filename, 
        edits_filedir, 
): 
    dist_stat = pd.DataFrame()
    dist_stat['cluster_distance'] = x
    dist_stat['n_sens_clusters'] = ys[0]
    dist_stat['n_resi_clusters'] = ys[1]

    fig, ax = plt.subplots()

    sns.lineplot(data=dist_stat, x="cluster_distance", y="n_sens_clusters")
    sns.lineplot(data=dist_stat, x="cluster_distance", y="n_resi_clusters")

    plt.xlabel('cluster radius')
    plt.ylabel('number of clusters')
    plt.title('sensitizing vs resistant clusters')
    plt.savefig(edits_filedir / f"cluster_LFC3D/cluster_distance.png") 
    dist_stat.to_csv(out_filename, sep = '\t', index=False)

def clustering_distance(
        df_struc_consvr, df_META, 
        thr_distance, 
        workdir, input_gene, input_uniprot, structureid, 
        i_affn='euclidean', i_link='single', n_clusters=None, pthr=0.05, 
):
    
    edits_filedir = Path(workdir + '/' + input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    df_hits_clust = df_META.copy()
    df_hits_clust["x_coord"] = df_struc_consvr["x_coord"]
    df_hits_clust["y_coord"] = df_struc_consvr["y_coord"]
    df_hits_clust["z_coord"] = df_struc_consvr["z_coord"]

    # SENSITIZING
    names = {'sensitizing':'SUM_LFC3D_neg_psig', 'resistant':'SUM_LFC3D_pos_psig'}
    for name, col in names.items(): # for sensitizing and resistant

        df_META_hits_coord = pd.DataFrame()
        df_META_temp = df_hits_clust.loc[(df_META[col] == 'p<'+str(pthr)), ]
        df_META_temp = df_META_temp.reset_index(drop=True)
        df_META_hits_coord['x_coord'] = df_META_temp['x_coord']
        df_META_hits_coord['y_coord'] = df_META_temp['y_coord']
        df_META_hits_coord['z_coord'] = df_META_temp['z_coord']
        np_META_hits_coord = np.array(df_META_hits_coord)

        func_clustering = AgglomerativeClustering(n_clusters=n_clusters, metric=i_affn, 
                                                  linkage=i_link, distance_threshold=thr_distance)

        if np_META_hits_coord.shape[0] < 2: 
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            raise None
        clustering = func_clustering.fit(np_META_hits_coord)
        clus_lbl = clustering.labels_
        n_c_output = int(max(clus_lbl) + 1)
        print(f'Number of clusters of {name} hits:', n_c_output)

        # fig = figure(figsize=(12, 8), dpi=300)
        dendogram_filename = edits_filedir / f"plots/{input_gene}_{input_uniprot}_{name}_hits_Dendogram_p_l001_{str(int(thr_distance))}A.png"
        fig = plot_dendrogram(clustering, df_META_temp, dendogram_filename, name)
        plt.savefig(dendogram_filename, dpi = 300)
        plt.show()

        # CLUSTER INDEX AND LENGTH
        hits_clust_filename = edits_filedir / f"cluster_LFC3D/{structureid}_MetaAggr_Hits_Clust_p_l001.tsv"
        df_META_hits_clust = pd.read_csv(hits_clust_filename, sep = '\t')
        df_META_clust = df_META_hits_clust.loc[(df_META_hits_clust[col] == 'p<'+str(pthr)), ]
        df_META_clust = df_META_clust.reset_index(drop=True)
        clust_indices = df_META_clust[f'{name}_hit_clust_{str(int(thr_distance))}'].unique()

        for c in clust_indices: 
            this_c_data = df_META_clust.loc[df_META_clust[f'{name}_hit_clust_{str(int(thr_distance))}'] == c, ]
            this_c_data = this_c_data.reset_index(drop=True)
            this_c_len = len(this_c_data)
            print(c, ':', this_c_len, ':', this_c_data.at[0, 'unipos'], '-', this_c_data.at[len(this_c_data)-1, 'unipos'])

    return clust_indices

def plot_dendrogram(
        clustering, df, 
        dendogram_filename, 
        name
):
    
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
        [clustering.children_, clustering.distances_, counts]
    ).astype(float)
    xlbl = list(df['unipos'])

    # plot the corresponding dendrogram
    dendrogram(linkage_matrix, color_threshold=6.0, 
               labels=xlbl, leaf_rotation=90.)
    plt.title(f'{name} clusters')
    plt.savefig(dendogram_filename, dpi=300) 
