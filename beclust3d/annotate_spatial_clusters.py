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
from matplotlib.pyplot import fig
from pathlib import Path

def clustering(
        df_struc_consvr, df_META, 
        workdir, input_gene, input_uniprot, structureid, 
        i_affn='euclidean', i_link='single', n_clusters=None,
):

    filedir = Path(workdir + input_gene)

    df_META["x_coord"] = df_struc_consvr["x_coord"]
    df_META["y_coord"] = df_struc_consvr["y_coord"]
    df_META["z_coord"] = df_struc_consvr["z_coord"]

    # CLUSTERING #
    arr_d_thr = [float(i+1) for i in range(20)]

    # SENSITIZING
    df_sens_hits_clust = pd.DataFrame()
    df_META_sens_hits = df_META.loc[(df_META['SUM_LFC3D_neg_psig'] == 'p<0.001'), ]
    df_META_sens_hits = df_META_sens_hits.reset_index(drop=True)
    arr_d_nc_sens = []
    df_sens_hits_clust['unipos'] = df_META_sens_hits['unipos']

    for thr_distance in arr_d_thr: # for every radius value 1.0 to 20.0
        colname_sens = "Sensitizing_hit_clust" + "_" + str(round(thr_distance))
        
        df_META_sens_hits_coord = pd.DataFrame()
        df_META_sens_hits_coord['x_coord'] = df_META_sens_hits['x_coord']
        df_META_sens_hits_coord['y_coord'] = df_META_sens_hits['y_coord']
        df_META_sens_hits_coord['z_coord'] = df_META_sens_hits['z_coord']
        np_META_sens_hits_coord = np.array(df_META_sens_hits_coord)
        
        func_clustering = AgglomerativeClustering(n_clusters=n_clusters, metric=i_affn, 
                                                  linkage=i_link, distance_threshold=thr_distance)
        clustering = func_clustering.fit(np_META_sens_hits_coord)
        clus_lbl = clustering.labels_
        n_c_output = int(max(clus_lbl)+1)
        print('Number of clusters of sensitizing hits:', '@ d = ', thr_distance, n_c_output)
        arr_d_nc_sens.append(n_c_output)
        df_META_sens_hits[colname_sens] = clus_lbl
        df_sens_hits_clust[colname_sens] = df_META_sens_hits[colname_sens]

    # RESISTANT
    df_resi_hits_clust = pd.DataFrame()
    df_META_resi_hits = df_META.loc[(df_META['SUM_LFC3D_pos_psig'] == 'p<0.001'), ]
    df_META_resi_hits = df_META_resi_hits.reset_index(drop=True)
    arr_d_nc_resi = [] 
    df_resi_hits_clust['unipos'] = df_META_resi_hits['unipos']

    for thr_distance in arr_d_thr:
        colname_resi = "Resistant_hit_clust" + "_" + str(round(thr_distance))
        
        df_META_resi_hits_coord = pd.DataFrame()
        df_META_resi_hits_coord['x_coord'] = df_META_resi_hits['x_coord']
        df_META_resi_hits_coord['y_coord'] = df_META_resi_hits['y_coord']
        df_META_resi_hits_coord['z_coord'] = df_META_resi_hits['z_coord']
        np_META_resi_hits_coord = np.array(df_META_resi_hits_coord)

        func_clustering = AgglomerativeClustering(n_clusters=n_clusters, metric=i_affn, 
                                                  linkage=i_link, distance_threshold=thr_distance)
        clustering = func_clustering.fit(np_META_resi_hits_coord)
        clus_lbl = clustering.labels_
        n_c_output = int(max(clus_lbl)+1)
        print('Number of clusters of resistant hits:', '@ d = ', thr_distance, n_c_output)
        arr_d_nc_resi.append(n_c_output)
        df_META_resi_hits[colname_resi] = clus_lbl
        df_resi_hits_clust[colname_resi] = df_META_resi_hits[colname_resi]

    # COMBINE
    df_META_sens_hits_clust = df_META.merge(df_sens_hits_clust, how='left', on=['unipos'])
    df_META_sens_resi_hits_clust = df_META_sens_hits_clust.merge(df_resi_hits_clust, how='left', on=['unipos'])
    df_META_sens_resi_hits_clust = df_META_sens_resi_hits_clust.fillna('-')

    meta_clust_filename = filedir / f"cluster_LFC3D/{input_gene}_{input_uniprot}_{structureid}_MetaAggr_Hits_Clust_p_l001.tsv"
    df_META_sens_resi_hits_clust.to_csv(meta_clust_filename, sep='\t', index=False)

    # PLOT
    cluster_distance_filename = filedir / f"cluster_LFC3D/{input_gene}_{input_uniprot}_{structureid}_MetaAggr_Hits_Clust_Dist_Stat_pl001.tsv"
    plot_cluster_distance(arr_d_thr, [arr_d_nc_sens, arr_d_nc_resi], 
                          cluster_distance_filename)

    return arr_d_thr, arr_d_nc_sens, arr_d_nc_resi

def plot_cluster_distance(
        x, ys, out_filename,
): 
    dist_stat = pd.DataFrame()
    dist_stat['cluster_distance'] = x
    dist_stat['n_sens_clusters'] = ys[0]
    dist_stat['n_resi_clusters'] = ys[1]

    sns.lineplot(data=dist_stat, x="cluster_distance", y="n_sens_clusters")
    sns.lineplot(data=dist_stat, x="cluster_distance", y="n_resi_clusters")
    ### rename x and y axis as cluster radius and number of clusters
    fig.savefig("cluster_distance.png") 
    dist_stat.to_csv(out_filename, sep = '\t', index=False)


def clustering_distance(
        df_META, 
        arr_d_nc_sens, arr_d_nc_resi, 
        thr_distance, 
        workdir, input_gene, input_uniprot, structureid, 
        i_affn='euclidean', i_link='single', n_clusters=None,
):
    
    filedir = Path(workdir + input_gene)

    # SENSITIZING

    df_META_sens_hits_coord = pd.DataFrame()
    df_META_sens_hits = df_META.loc[(df_META['SUM_LFC3D_neg_psig'] == 'p<0.001'), ]
    df_META_sens_hits = df_META_sens_hits.reset_index(drop=True)
    df_META_sens_hits_coord['x_coord'] = df_META_sens_hits['x_coord']
    df_META_sens_hits_coord['y_coord'] = df_META_sens_hits['y_coord']
    df_META_sens_hits_coord['z_coord'] = df_META_sens_hits['z_coord']
    np_META_sens_hits_coord = np.array(df_META_sens_hits_coord)

    func_clustering = AgglomerativeClustering(n_clusters=n_clusters, metric=i_affn, 
                                              linkage=i_link, distance_threshold=thr_distance)
    clustering = func_clustering.fit(np_META_sens_hits_coord)
    clus_lbl = clustering.labels_
    n_c_output = int(max(clus_lbl) + 1)
    print('Number of clusters of sensitizing hits:', n_c_output)
    arr_d_nc_sens.append(n_c_output)

    figure = fig(figsize=(12, 8), dpi=300)
    figure = plot_dendrogram(clustering, df_META_sens_hits)

    sens_clust_dendogram_filename = filedir / f"plots/{input_gene}_{input_uniprot}_SensHits_Dendogram_p_l001_6A.png"
    plt.savefig(sens_clust_dendogram_filename, dpi = 300) 

    # CLUSTER INDEX AND LENGTH

    meta_clust_filename = filedir / f"cluster_LFC3D/{input_gene}_{input_uniprot}_{structureid}_MetaAggr_Hits_Clust_p_l001.tsv"
    df_META_sens_resi_hits_clust = pd.read_csv(meta_clust_filename, sep = '\t')
    df_META_sens_clust = df_META_sens_resi_hits_clust.loc[(df_META_sens_resi_hits_clust['SUM_LFC3D_neg_psig'] == 'p<0.001'), ]
    df_META_sens_clust = df_META_sens_clust.reset_index(drop=True)
    sens_clust_indices = df_META_sens_clust['Sensitizing_hit_clust'].unique()

    print('Sensitizing cluster index', ':', 'length')
    for c in sens_clust_indices:
        this_c_data = df_META_sens_clust.loc[df_META_sens_clust['Sensitizing_hit_clust'] == c, ]
        this_c_data = this_c_data.reset_index(drop=True)
        this_c_len = len(this_c_data)
        print(c, ':', this_c_len, ':', this_c_data.at[0, 'unipos'], '-', this_c_data.at[len(this_c_data)-1, 'unipos'])

    # RESISTANT

    df_META_resi_hits_coord = pd.DataFrame()
    df_META_resi_hits = df_META.loc[(df_META['SUM_LFC3D_pos_psig'] == 'p<0.001'), ]
    df_META_resi_hits = df_META_resi_hits.reset_index(drop=True)
    df_META_resi_hits_coord['x_coord'] = df_META_resi_hits['x_coord']
    df_META_resi_hits_coord['y_coord'] = df_META_resi_hits['y_coord']
    df_META_resi_hits_coord['z_coord'] = df_META_resi_hits['z_coord']
    np_META_resi_hits_coord = np.array(df_META_resi_hits_coord)

    func_clustering = AgglomerativeClustering(n_clusters=n_clusters, metric=i_affn, 
                                              linkage=i_link, distance_threshold=thr_distance)
    clustering = func_clustering.fit(np_META_resi_hits_coord)
    clus_lbl = clustering.labels_
    n_c_output = int(max(clus_lbl) + 1)
    print('Number of clusters of resistant hits:', n_c_output)
    arr_d_nc_resi.append(n_c_output)

    figure = fig(figsize=(12, 8), dpi=300)
    figure = plot_dendrogram(clustering, df_META_resi_hits)

    resi_clust_dendogram_filename = filedir / f"plots/{input_gene}_{input_uniprot}_ResiHits_Dendogram_p_l001_6A.png"
    plt.savefig(resi_clust_dendogram_filename, dpi = 300) 

    # CLUSTER INDEX AND START/END

    meta_clust_filename = filedir / f"cluster_LFC3D/{input_gene}_{input_uniprot}_{structureid}_MetaAggr_Hits_Clust_p_l001.tsv"
    df_META_sens_resi_hits_clust = pd.read_csv(meta_clust_filename, sep = '\t')    
    df_META_resi_clust = df_META_sens_resi_hits_clust.loc[(df_META_sens_resi_hits_clust['SUM_LFC3D_pos_psig'] == 'p<0.001'), ]
    df_META_resi_clust = df_META_resi_clust.reset_index(drop=True)
    resi_clust_indices = df_META_resi_clust['Resistant_hit_clust'].unique()

    print('Resistant cluster index', ':', 'length')
    for c in resi_clust_indices:
        this_c_data = df_META_resi_clust.loc[df_META_resi_clust['Resistant_hit_clust'] == c, ]
        this_c_data = this_c_data.reset_index(drop=True)
        this_c_len = len(this_c_data)
        print(c, ':', this_c_len, ':', this_c_data.at[0, 'unipos'], '-', this_c_data.at[len(this_c_data)-1, 'unipos'])


def plot_dendrogram(clustering, df, **kwargs):
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
               labels=xlbl, leaf_rotation=90., **kwargs)
