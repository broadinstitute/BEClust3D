import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings

from variables import *
from beclust3d.prioritize_by_conservation import prioritize_by_conservation

@pytest.mark.parametrize(("gene", "structid", "mouse_gene"), zip(all_genes, all_structureids, all_mouse_genes))
@pytest.mark.parametrize("screen", all_human_screens)
def test_afstructuralfeatures_human(gene, structid, mouse_gene, screen): 

    df_struc = pd.read_csv(f"{workdir}/{gene}/{structid}_coord_struc_features.tsv", sep = "\t")
    df_consrv = pd.read_csv(f"{workdir}/{gene}/Human{gene}_Mouse{mouse_gene}_residuemap_conservation.tsv", sep = '\t')

    prioritize_by_conservation(
        df_struc     =df_struc, 
        df_consrv    =df_consrv, 
        workdir      =workdir, 
        input_gene   =gene, 
        input_screen =screen, 
        structureid  =structid, 
    )

    screen_name = screen.split('.')[0]
    for f in os.listdir(f'{workdir}/{gene}/screendata'): 
        new_file = f'{workdir}/{gene}/screendata/{f}'
        ref_file = f'{refdir}/{gene}/screendata/{f}'
        if os.path.exists(new_file) and os.path.exists(ref_file): 
            if 'struc_consrv' in new_file and 'struc_consrv' in ref_file: 
                if screen_name in new_file and screen_name in ref_file: 
                    similarity_ratio = fuzzy_compare(new_file, ref_file)
                    if 0.9 > similarity_ratio > 0.8: 
                        warnings.warn(f"Similarity Ratio between {new_file} and {ref_file} are moderate")
                    assert similarity_ratio >= 0.8
                    # 0.9 is the arbitrary threshold of similarity

                    ### validate if this is a good threshold
