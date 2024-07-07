"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2

"""

import pandas as pd
from pathlib import Path
import math
import os

def prioritize_by_conservation(
        df_struc, df_consrv, 
        workdir, 
        input_gene, input_screen, structureid, 
): 
    """
    Description
        Takes in results across multiple categories for a screen, and
        prioritizes residue scores by conservation. 

    Params
        df_struc: pandas dataframe, required
            DataFrame output from conservation()
        df_consrv: pandas dataframe, required
            DataFrame output from conservation()
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        structureid: str, required
            the name of the AF and uniprot input

    Returns
        df_struc_consvr: pandas dataframe
            DataFrame
    """

    df_struc_consvr = df_struc
    df_struc_consvr['human_res_pos'] = df_consrv['human_res_pos']
    df_struc_consvr['mouse_res_pos'] = df_consrv['mouse_res_pos']
    df_struc_consvr['mouse_res']     = df_consrv['mouse_res']
    df_struc_consvr['conservation']  = df_consrv['conservation']

    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')

    edits_filename = edits_filedir / f"screendata/{input_gene}_{structureid}_struc_consrv.tsv"
    df_struc_consvr.to_csv(edits_filename, 
                           sep = "\t", index=False)

    # for Missense, Silent, Nonsense
    for edit_type in ['Missense', 'Silent', 'Nonsense']: 

        edits_filename = edits_filedir / f"screendata/{input_gene}_{screen_name}_{edit_type}_edits_list.tsv"
        df_edit = pd.read_csv(edits_filename, sep='\t')
        
        arr_unique_LFC = []
        arr_all_edits = []

        for i in range(0, len(df_struc_consvr)):
            human_res_pos = df_struc_consvr.at[i, 'human_res_pos']
            unique_LFC_per_residue_human = '-'
            all_edits_per_residue_human = '-'
            
            df_pos_edits = df_edit.loc[df_edit['human_pos'] == int(human_res_pos), ]
            df_pos_edits = df_pos_edits.reset_index()

            if len(df_pos_edits) > 1:
                # edits
                pos_edits_list = df_pos_edits['edit'].tolist()
                all_edits_per_residue_human = ';'.join(list(set(pos_edits_list)))
                # scores
                pos_LFCscore_list = df_pos_edits['LFC'].tolist()
                pos_LFCscore = round(sum(pos_LFCscore_list) / len(pos_LFCscore_list), 3)
                unique_LFC_per_residue_human = str(pos_LFCscore)
            elif len(df_pos_edits) == 1:   
                all_edits_per_residue_human = df_pos_edits.at[0,'edit']
                unique_LFC_per_residue_human = str(round(df_pos_edits.at[0,'LFC'], 3))
            else:
                all_edits_per_residue_human = '-'
                unique_LFC_per_residue_human = '-'

            arr_unique_LFC.append(unique_LFC_per_residue_human)
            arr_all_edits.append(all_edits_per_residue_human)

        df_struc_consvr[f'all_{edit_type}_edits'] = arr_all_edits
        df_struc_consvr[f'mean_{edit_type}_LFC'] = arr_unique_LFC

    out_filename = edits_filedir / f"screendata/{input_gene}_{screen_name}_struc_consrv_proteinedits.tsv"
    df_struc_consvr.to_csv(out_filename, sep = '\t', index=False)

    return df_struc_consvr
