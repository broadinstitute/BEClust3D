import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings
from variables import *
from beclust3d.average_split_lfc3d import average_and_split




@pytest.mark.parametrize(("gene", "uniprot"), zip(all_genes, all_uniprots))
@pytest.mark.parametrize("screen", all_human_screens)
def test_afstructuralfeatures_human(gene, uniprot, screen): 

    screen_name = screen.split('.')[0]
    filename = f"{workdir}/{gene}/LFC3D/{gene}_{screen_name}_LFC_LFC3D_per_Random_LFC3Dr.tsv"
    if not os.path.exists(filename): 
        warnings.warn(f"{filename} does not exist")
        return True
    
    df_str_cons_3daggr = pd.read_csv(filename, 
                                    sep = "\t")
                                    
    res = average_and_split(
                df_LFC_LFCrN_LFC3D_LFC3DrN=df_str_cons_3daggr, 
                workdir=workdir, 
                input_gene=gene,
                input_screen=screen, 
                )

    if res is not None: 
        assert f'{gene}_{screen_name}_LFC_LFC3D_LFC3Dr_bidirectional.tsv' in os.listdir(f'{workdir}/{gene}/LFC3D')
