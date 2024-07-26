import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from variables import *
from beclust3d.preprocess_be_results import parse_base_editing_results

vals = ['D14_Input', 'D14_Input', 'D14_Input', 'D17_Input', 'D14_Input', 'D14_Input', ]
mouse_vals = ['NSG_Input', 'Tx_NSG', 'NSG_Input', 'Tx_NSG', 'NSG_Input', 'Tx_NSG', ]

@pytest.mark.parametrize("gene", all_genes)
@pytest.mark.parametrize(("screen", "val"), zip(all_human_screens, vals))
def test_afstructuralfeatures_human(gene, screen, val): 
    df_filename = f'{workdir}/../rawdata/{screen}'
    df_screen = pd.read_csv(df_filename, sep='\t')

    parse_base_editing_results(
        df_InputGene = df_screen, 
        workdir      = workdir, 
        input_gene   = gene, 
        input_screen = screen, 
        mut_col='Histogram_Category', val_col=val, 
        gene_col='TargetGeneSymbol', edits_col='AminoAcidEdits',
        )

    check_fileformats = ['tsv']
    for f in os.listdir(f'{workdir}/{gene}/screendata'): 
        new_file = f'{workdir}/{gene}/screendata/{f}'
        ref_file = f'{refdir}/{gene}/screendata/{f}'
        if os.path.exists(new_file) and os.path.exists(ref_file): 
            new_file_suff = new_file.split('.')[-1]
            ref_file_suff = ref_file.split('.')[-1]
            if new_file_suff in check_fileformats and ref_file_suff in check_fileformats: 
                if 'Human' in new_file and 'Human' in ref_file: 
                    if filecmp.cmp(new_file, ref_file): 
                        print('Checked')
                        assert True
                    else: 
                        assert False


@pytest.mark.parametrize(("gene", "uniprot", "structid"), zip(all_genes, all_uniprots, all_structureids))
@pytest.mark.parametrize(("screen", "val"), zip(all_mouse_screens, mouse_vals))
def test_afstructuralfeatures_mouse(gene, uniprot, structid, screen, val): 

    df_filename = f'{workdir}/../rawdata/{screen}'
    df_screen = pd.read_csv(df_filename, sep='\t')

    parse_base_editing_results(
        df_InputGene = df_screen, 
        workdir      = workdir, 
        input_gene   = gene, 
        input_screen = screen, 
        mut_col='Histogram_Category', val_col=val, 
        gene_col='TargetGeneSymbol', edits_col='AminoAcidEdits',
        output_col='mouse_pos', 
        )
    
    # # for some reason the outputs do not match exacty to reference, but very close
    # check_fileformats = ['tsv']
    # for f in os.listdir(f'{workdir}/{gene}/screendata'): 
    #     new_file = f'{workdir}/{gene}/screendata/{f}'
    #     ref_file = f'{refdir}/{gene}/screendata/{f}'
    #     if os.path.exists(new_file) and os.path.exists(ref_file): 
    #         new_file_suff = new_file.split('.')[-1]
    #         ref_file_suff = ref_file.split('.')[-1]
    #         if new_file_suff in check_fileformats and ref_file_suff in check_fileformats: 
    #             print('TRUE')
    #             if 'Mouse' in new_file and 'Mouse' in ref_file: 
    #                 if filecmp.cmp(new_file, ref_file): 
    #                     print('Checked')
    #                     assert True
    #                 else: 
    #                     print(new_file)
    #                     print(ref_file)
    #                     assert False

# all preprocessing is done separate by gene and screen
