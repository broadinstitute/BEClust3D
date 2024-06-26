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

                df_pos_edits = df_edit.loc[df_edit['human_pos'] == int(human_res_pos), ]
                df_pos_edits = df_pos_edits.reset_index()

                if len(df_pos_edits) > 1:
                    # edit
                    pos_edits_list = df_pos_edits['edit'].tolist()
                    pos_edits_str = ';'.join(list(set(pos_edits_list)))
                    # score
                    pos_LFCscore_list = df_pos_edits['LFC'].tolist()
                    pos_LFCscore = round(sum(pos_LFCscore_list) / len(pos_LFCscore_list), 3)
                    pos_LFCscore_str = str(pos_LFCscore)

                elif len(df_pos_edits) == 1:   
                    pos_edits_str = df_pos_edits.at[0,'edit']
                    pos_LFCscore_str = str(round(df_pos_edits.at[0,'LFC'], 3))
                else:
                    pos_edits_str = '-'
                    pos_LFCscore_str = '-'

                unique_LFC_per_residue_human = pos_LFCscore_str
                all_edits_per_residue_human = pos_edits_str

            arr_unique_LFC.append(unique_LFC_per_residue_human)
            arr_all_edits.append(all_edits_per_residue_human)

        df_struc_consvr[f'all_{edit_type}_edits'] = arr_all_edits
        df_struc_consvr[f'mean_{edit_type}_LFC'] = arr_unique_LFC

    edits_filedir = Path(workdir + input_gene)
    out_filename = edits_filedir / f"screendata/{input_gene}_{structureid}_struc_consrv_proteinedits.tsv"
    df_struc_consvr.to_csv(out_filename, sep = '\t', index=False)

    return df_struc_consvr
