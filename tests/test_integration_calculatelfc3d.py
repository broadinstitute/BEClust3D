import pandas as pd
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from variables import *
from beclust3d.calculate_lfc3d import calculate_lfc3d


@pytest.mark.parametrize(("gene", "structid"), zip(all_genes, all_structureids))
def test_afstructuralfeatures_human(gene, structid): 

    df_str_cons = pd.read_csv(f"{workdir}/{gene}/{structid}_coord_struc_features.tsv", sep = "\t")

    res = calculate_lfc3d(
        df_str_cons  =df_str_cons, 
        workdir      =workdir, 
        input_gene   =gene, 
        input_screens=all_human_screens, 
    )

    if res is not None: 
        assert f'{gene}_LFC_LFC3D_per_Random_LFC3Dr.tsv' in os.listdir(f'{workdir}/{gene}/LFC3D')

### MOUSE screens