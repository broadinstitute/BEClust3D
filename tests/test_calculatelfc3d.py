import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings

from variables import *
from beclust3d.calculate_lfc3d import calculate_lfc3d


@pytest.mark.parametrize(("gene", "uniprot", "structid"), zip(all_genes, all_uniprots, all_structureids))
@pytest.mark.parametrize("screen", all_human_screens)
def test_afstructuralfeatures_human(gene, uniprot, structid, screen): 

    df_str_cons = pd.read_csv(f"{workdir}/{gene}/{structid}_coord_struc_features.tsv", sep = "\t")

    res = calculate_lfc3d(
        df_str_cons  =df_str_cons, 
        workdir      =workdir, 
        input_gene   =gene, 
        input_screen =screen, 
    )

    if res is not None: 
        screen_name = screen.split('.')[0]
        assert f'{gene}_{screen_name}_LFC_LFC3D_per_Random_LFC3Dr.tsv' in os.listdir(f'{workdir}/{gene}/LFC3D')

### MOUSE screens
