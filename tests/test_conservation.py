import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from variables import *
from beclust3d.conservation import conservation

email='xiaohu@g.harvard.edu'
title='samplejob'

@pytest.mark.parametrize(("gene", "uniprot", "mouse_gene", "mouse_uniprot"), 
                         zip(all_genes, all_uniprots, all_mouse_genes, all_mouse_uniprots))
def test_afstructuralfeatures_human(gene, uniprot, mouse_gene, mouse_uniprot): 
    conservation(
        workdir = workdir, 
        input_human_gene = gene, input_mouse_gene = mouse_gene, 
        input_human_uniid = uniprot, input_mouse_uniid = mouse_uniprot, 
        email = email, title = title, 
        )

    check_fileformats = ['tsv', 'fasta']
    for f in os.listdir(f'{workdir}/{gene}'): 
        new_file = f'{workdir}/{gene}/{f}'
        ref_file = f'{refdir}/{gene}/{f}'
        if os.path.exists(new_file) and os.path.exists(ref_file): 
            new_file_suff = new_file.split('.')[-1]
            ref_file_suff = ref_file.split('.')[-1]
            if new_file_suff in check_fileformats and ref_file_suff in check_fileformats: 
                if filecmp.cmp(new_file, ref_file): 
                    assert True
                else: 
                    assert False
