import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings

from variables import *
from beclust3d.metaaggregation import metaaggregation


@pytest.mark.parametrize(("gene", "uniprot", "structid", "mouse_gene"), zip(all_genes, all_uniprots, all_structureids, all_mouse_genes))
@pytest.mark.parametrize("screen", all_human_screens)
def test_afstructuralfeatures_human(gene, uniprot, structid, mouse_gene, screen): 

    screen_name = screen.split('.')[0]
    filename = f'{workdir}/{gene}/LFC3D/{gene}_{screen_name}_LFC_LFC3D_per_Random_LFC3Dr.tsv'
    if not os.path.exists(filename): 
        warnings.warn(f"{filename} does not exist")
        return True
    
    df_LFC_LFC3D = pd.read_csv(filename, sep = "\t")
    metaaggregation(
            df_LFC_LFC3D=df_LFC_LFC3D, 
            workdir=workdir, 
            input_gene=gene, 
            input_uniprot=uniprot, 
            structureid=structid, 
            input_screen=screen,
            )

    assert f'{structid}_{screen_name}_MetaAggr_LFC3D_and_randomized_background.tsv' in os.listdir(f'{workdir}/{gene}/metaaggregation')
    plot_dir = os.listdir(f'{workdir}/{gene}/plots')
    assert f'{gene}_{screen_name}_signal_vs_background_sensitizing.png' in plot_dir
    assert f'{gene}_{screen_name}_signal_vs_background_resistant.png' in plot_dir
    assert f'{gene}_{screen_name}_Aggr_LFC3D_neg_pvalue_histogram.png' in plot_dir
    assert f'{gene}_{screen_name}_Aggr_LFC3D_pos_pvalue_histogram.png' in plot_dir
    assert f'{gene}_{screen_name}_Aggr_LFC3D_neg_dis_dot_per_residue.png' in plot_dir
    assert f'{gene}_{screen_name}_Aggr_LFC3D_pos_dis_dot_per_residue.png' in plot_dir
