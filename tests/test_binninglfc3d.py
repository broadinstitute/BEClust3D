import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings

from variables import *
from beclust3d.binning_lfc3d import binning_lfc3d


@pytest.mark.parametrize("gene", all_genes)
@pytest.mark.parametrize("screen", all_human_screens)
def test_binninglfc3d_human(gene, screen): 

    filename = f"{workdir}/{gene}/LFC3D/{gene}_LFC_LFC3D_LFC3Dr_bidirectional.tsv"
    if not os.path.exists(filename): 
        warnings.warn(f"{filename} does not exist")
        return True
    
    df_bidir = pd.read_csv(filename, sep = "\t")
    res = binning_lfc3d(
                 df_LFC_LFC3D=df_bidir, 
                 workdir=workdir, 
                 input_gene=gene,
                 input_screens=[screen], 
                 )

    if res is not None: 
        assert f'{gene}_LFC_LFC3D_dis_wght_Signal_Only_per_Screen.tsv' in os.listdir(f'{workdir}/{gene}/LFC3D')
