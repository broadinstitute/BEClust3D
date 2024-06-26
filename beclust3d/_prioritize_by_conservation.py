"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2

"""

import pandas as pd
from pathlib import Path
import math


def prioritize_by_conservation(df_struc, df_consrv, 
                               workdir, input_gene, structureid, 
                               screenid, 
                               ): 

    df_struc_consvr = df_struc
    df_struc_consvr['human_res_pos'] = df_consrv['human_res_pos']
    df_struc_consvr['mouse_res_pos'] = df_consrv['mouse_res_pos']
    df_struc_consvr['mouse_res']     = df_consrv['mouse_res']
    df_struc_consvr['conservation']  = df_consrv['conservation']

    edits_filedir = Path(workdir + input_gene)
    edits_filename = edits_filedir / f"screendata/{input_gene}_{structureid}_struc_consrv.tsv"
    df_struc_consvr.to_csv(edits_filename, sep = "\t", index=False)

    # for missense, silent, nonsense
    for edit_type in ['missense', 'silent', 'nonsense']: 

        edits_filename = edits_filedir / f"screendata/{input_gene}_{screenid}_{edit_type}_edits_list.tsv"
        df_edit = pd.read_csv(edits_filename, sep = '\t')
        
        arr_unique_LFC = []
        arr_all_edits = []

        for i in range(0, len(df_struc_consvr)):
            human_res_pos = df_struc_consvr.at[i, 'human_res_pos']
            
            unique_LFC_per_residue_human = '-'
            all_edits_per_residue_human = '-'
            
            if not math.isnan(human_res_pos): ### 

                df_this_humanpos_edits = df_edit.loc[df_edit['human_pos'] == int(human_res_pos), ]
                df_this_humanpos_edits = df_this_humanpos_edits.reset_index()

                if len(df_this_humanpos_edits) > 1:
                    #edit
                    this_humanpos_edits_list = df_this_humanpos_edits['edit'].tolist()
                    this_humanpos_edits_str = ';'.join(list(set(this_humanpos_edits_list)))
                    #score
                    this_humanpos_LFCscore_list = df_this_humanpos_edits['LFC'].tolist()
                    this_humanpos_LFCscore = round(sum(this_humanpos_LFCscore_list) / len(this_humanpos_LFCscore_list), 3)
                    this_humanpos_LFCscore_str = str(this_humanpos_LFCscore)

                elif len(df_this_humanpos_edits) == 1:   
                    this_humanpos_edits_str = df_this_humanpos_edits.at[0,'edit']
                    this_humanpos_LFCscore_str = str(round(df_this_humanpos_edits.at[0,'LFC'], 3))
                else:
                    this_humanpos_edits_str = '-'
                    this_humanpos_LFCscore_str = '-'

                unique_LFC_per_residue_human = this_humanpos_LFCscore_str
                all_edits_per_residue_human = this_humanpos_edits_str

            arr_unique_LFC.append(unique_LFC_per_residue_human)
            arr_all_edits.append(all_edits_per_residue_human)

        df_struc_consvr[f'all_{edit_type}_edits'] = arr_all_edits
        df_struc_consvr[f'mean_{edit_type}_LFC'] = arr_unique_LFC

    return df_struc_consvr
