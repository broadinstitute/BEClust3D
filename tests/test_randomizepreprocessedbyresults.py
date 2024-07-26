import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from variables import *
from beclust3d.randomize_preprocessed_be_results import randomize_be_results

@pytest.mark.parametrize("gene", all_genes)
@pytest.mark.parametrize("screen", all_human_screens)
def test_afstructuralfeatures_human(gene, screen): 
    screen_name = screen.split('.')[0]
    df_filename = f'{workdir}/{gene}/screendata/{gene}_{screen_name}_Missense_edits_list.tsv'
    df_missense = pd.read_csv(df_filename, sep='\t')

    randomize_be_results(
        df_missense  = df_missense, 
        workdir      = workdir, 
        input_gene   = gene, 
        input_screen = screen
        )

    assert f'{gene}_{screen_name}_missense_edits_randomized.tsv' in os.listdir(f'{workdir}/{gene}/randomized_screendata')
