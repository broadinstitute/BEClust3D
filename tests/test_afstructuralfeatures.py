import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from variables import *
from beclust3d.af_structural_features import af_structural_features

@pytest.mark.parametrize(("gene", "uniprot", "structid"), zip(all_genes, all_uniprots, all_structureids))
def test_afstructuralfeatures(gene, uniprot, structid): 
    af_structural_features(
        workdir       = workdir, 
        input_gene    = gene, 
        input_uniprot = uniprot, 
        structureid   = structid, 
    )

    check_fileformats = ['pdb', 'tsv', 'dssp', 'fasta']
    for f in os.listdir(f'{workdir}/{gene}'): 
        new_file = f'{workdir}/{gene}/{f}'
        ref_file = f'{refdir}/{gene}/{f}'
        if os.path.exists(new_file) and os.path.exists(ref_file): 
            new_file_suff = new_file.split('.')[-1]
            ref_file_suff = ref_file.split('.')[-1]
            if new_file_suff in check_fileformats and ref_file_suff in check_fileformats: 
                assert filecmp.cmp(new_file, ref_file)
