import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings

from variables import *
from beclust3d.annotate_spatial_clusters import clustering, clustering_distance


@pytest.mark.parametrize(("gene", "uniprot", "structid"), zip(all_genes, all_uniprots, all_structureids))
def test_annotatespatialclusters_integration_human(gene, uniprot, structid): 

    filename_meta = f'{workdir}/{gene}/metaaggregation/{structid}_MetaAggr_LFC3D_and_randomized_background.tsv'
    if not os.path.exists(filename_meta): 
        warnings.warn(f"{filename_meta} does not exist")
        return True
    filename_coordstruc = f'{workdir}/{gene}/{structid}_coord_struc_features.tsv'
    if not os.path.exists(filename_coordstruc): 
        warnings.warn(f"{filename_coordstruc} does not exist")
        return True
    
    df_str_cons = pd.read_csv(filename_coordstruc, sep = "\t")
    df_META = pd.read_csv(filename_meta, sep = "\t")
    dists, arrs = clustering(
        df_struc_consvr=df_str_cons, 
        df_META=df_META, 
        workdir=workdir, 
        input_gene=gene, 
        structureid=structid,
        )
    
    if dists is None and arrs is None: 
        return True

    assert f'{structid}_MetaAggr_Hits_Clust_Dist_Stat_pl001.tsv' in os.listdir(f'{workdir}/{gene}/cluster_LFC3D')
    
    clust_dist = 6.0
    res = clustering_distance(
        df_struc_consvr=df_str_cons, 
        df_META=df_META, 
        thr_distance=clust_dist, 
        workdir=workdir, 
        input_gene=gene, 
        input_uniprot=uniprot, 
        structureid=structid, 
        )
    
    if res is None: 
        return True

    assert f"{gene}_{uniprot}_sensitizing_hits_Dendogram_p_l001_6A.png" in os.listdir(f'{workdir}/{gene}/plots')
    assert f"{gene}_{uniprot}_resistant_hits_Dendogram_p_l001_6A.png" in os.listdir(f'{workdir}/{gene}/plots')
