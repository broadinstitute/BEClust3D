import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from variables import *
from beclust3d.randomize_by_conservation import randomize_by_conservation

@pytest.mark.parametrize(("gene", "structid"), zip(all_genes, all_structureids))
@pytest.mark.parametrize("screen", all_human_screens)
def test_afstructuralfeatures_human(gene, structid, screen): 
    res = randomize_by_conservation(
        workdir      =workdir, 
        input_gene   =gene, 
        input_screen =screen, 
        structureid  =structid, 
    )

    if res is not None: 
        screen_name = screen.split('.')[0]
        assert f'{gene}_{screen_name}_struc_consrv_missenseedits_randomized.tsv' in os.listdir(f'{workdir}/{gene}/randomized_screendata')
