"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2

"""

import pandas as pd
from pathlib import Path
import os

def prioritize_by_conservation(
        df_struc, df_consrv, 
        workdir, 
        input_gene, input_screen, structureid, 
): 
    """
    Description
        Takes in results across multiple edit types for a screen, and
        aggregates the edits for each residue with conservation information. 

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
    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')

    df_struc_consvr = df_struc
    df_struc_consvr['human_res_pos'] = df_consrv['human_res_pos']
    df_struc_consvr['mouse_res_pos'] = df_consrv['mouse_res_pos']
    df_struc_consvr['mouse_res']     = df_consrv['mouse_res']
    df_struc_consvr['conservation']  = df_consrv['conservation']
    del df_struc, df_consrv

    struc_consrv_filename =  f"screendata/{input_gene}_{structureid}_struc_consrv.tsv"
    df_struc_consvr.to_csv(edits_filedir / struc_consrv_filename, sep = "\t", index=False)

    # FOR EACH EDIT TYPE, AGGREGATE LFC AND EDITS WITH CONSERVATION #
    ### this can be done without conservation, just with an input dataframe that only has the original sequence
    for edit_type in ['Missense', 'Silent', 'Nonsense']: 

        in_filename = f"screendata/{input_gene}_{screen_name}_{edit_type}_edits_list.tsv"
        df_edit = pd.read_csv(edits_filedir / in_filename, sep='\t')
        
        arr_unique_LFC = []
        arr_all_edits = []

        # FOR EACH RESIDUE #
        for i in range(len(df_struc_consvr)): 
            human_res_pos = df_struc_consvr.at[i, 'human_res_pos'] ###
            df_pos_edits = df_edit.loc[df_edit['edit_pos'] == int(human_res_pos), ].reset_index() ###

            if len(df_pos_edits) > 1: 
                pos_LFCscore_list = df_pos_edits['LFC'].tolist()
                pos_LFCscore = round(sum(pos_LFCscore_list) / len(pos_LFCscore_list), 3)
                unique_LFC_res = str(pos_LFCscore)

                pos_edits_list = df_pos_edits['this_edit'].tolist()
                all_edits_res = ';'.join(list(set(pos_edits_list)))
            elif len(df_pos_edits) == 1:   
                unique_LFC_res = str(round(df_pos_edits.at[0, 'LFC'], 3))
                all_edits_res = df_pos_edits.at[0, 'this_edit']
            else:
                unique_LFC_res, all_edits_res = '-', '-'

            arr_unique_LFC.append(unique_LFC_res)
            arr_all_edits.append(all_edits_res)

        df_struc_consvr[f'mean_{edit_type}_LFC'] = arr_unique_LFC
        df_struc_consvr[f'all_{edit_type}_edits'] = arr_all_edits

    strcons_edits_filename = f"screendata/{input_gene}_{screen_name}_struc_consrv_proteinedits.tsv"
    df_struc_consvr.to_csv(edits_filedir / strcons_edits_filename, sep = '\t', index=False)

    return df_struc_consvr
