"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2.5

"""

import pandas as pd
from pathlib import Path
import math

# instead of parse loop and randomizing loop, can probably be one loop

def randomize_by_conservation(df_struc, df_consrv, 
        workdir, input_human_gene, input_human_uniprot, 
        screenid, 
        nRandom=200, 
): 
    
    structureid = f"AF-{input_human_uniprot}-F1-model_v4"

    df_struc_consvr = df_struc
    df_struc_consvr['human_res_pos'] = df_consrv['human_res_pos']
    df_struc_consvr['mouse_res_pos'] = df_consrv['mouse_res_pos']
    df_struc_consvr['mouse_res']     = df_consrv['mouse_res']
    df_struc_consvr['conservation']  = df_consrv['conservation']

    edits_filedir = Path(workdir + input_human_gene)
    edits_filename = edits_filedir / f"screendata/{input_human_gene}_{structureid}_struc_consrv.tsv"
    df_struc_consvr.to_csv(edits_filename, sep = "\t", index=False)
                
    # missense
    missense_filename = edits_filedir / f"randomized_screendata/{input_human_gene}_{screenid}_missense_edits_randomized.tsv"
    df_missense = pd.read_csv(missense_filename, sep = '\t');

    arr_unique_LFC = []
    arr_all_edits = []

    for i in range(0, len(df_struc_consvr)):
        human_res_pos = df_struc_consvr.at[i, 'human_res_pos']

        unique_LFC_per_residue_human = '-'
        all_edits_per_residue_human = '-'

        if not math.isnan(human_res_pos): ### 

            df_pos_edits = df_missense.loc[df_missense['human_pos'] == int(human_res_pos), ]
            df_pos_edits = df_pos_edits.reset_index()

            if len(df_pos_edits) > 1:
                #edit
                pos_edits_list = df_pos_edits['edit'].tolist()
                pos_edits_str = ';'.join(list(set(pos_edits_list)))
                #score
                pos_LFCscore_list = df_pos_edits['LFC'].tolist()
                pos_LFCscore = round(sum(pos_LFCscore_list) / len(pos_LFCscore_list),3)
                pos_LFCscore_str = str(pos_LFCscore)

            elif len(df_pos_edits) == 1:   
                pos_edits_str = df_pos_edits.at[0,'edit']
                pos_LFCscore_str = str(round(df_pos_edits.at[0,'LFC'],3))
            else:
                pos_edits_str = '-'
                pos_LFCscore_str = '-'

            unique_LFC_per_residue_human = pos_LFCscore_str
            all_edits_per_residue_human = pos_edits_str

        arr_unique_LFC.append(unique_LFC_per_residue_human)
        arr_all_edits.append(all_edits_per_residue_human)

    df_struc_consvr['all_missense_edits'] = arr_all_edits
    df_struc_consvr['mean_missense_LFC'] = arr_unique_LFC
    
    for r in range(0, nRandom):
        arr_unique_LFC = []
        arr_all_edits = []

        for i in range(0, len(df_struc_consvr)):
            human_res_pos = df_struc_consvr.at[i, 'human_res_pos']

            unique_LFC_per_residue_human = '-'
            all_edits_per_residue_human = '-'

            if not math.isnan(human_res_pos): ### 
                
                df_pos_edits = df_missense.loc[df_missense['human_pos'] == int(human_res_pos), ]
                df_pos_edits = df_pos_edits.reset_index()

                header_LFC_rN = "LFCr" + str(r+1);
                if len(df_pos_edits) > 1:
                    # edit
                    pos_edits_list = df_pos_edits['edit'].tolist()
                    pos_edits_str = ';'.join(list(set(pos_edits_list)))
                    # score
                    pos_LFCscore_list = df_pos_edits[header_LFC_rN].tolist()
                    pos_LFCscore = round(sum(pos_LFCscore_list) / len(pos_LFCscore_list), 3)
                    pos_LFCscore_str = str(pos_LFCscore)

                elif len(df_pos_edits) == 1:   
                    pos_edits_str = df_pos_edits.at[0,'edit']
                    pos_LFCscore_str = str(round(df_pos_edits.at[0,header_LFC_rN], 3))
                else:
                    pos_edits_str = '-'
                    pos_LFCscore_str = '-'

                unique_LFC_per_residue_human = pos_LFCscore_str
                all_edits_per_residue_human = pos_edits_str

            arr_unique_LFC.append(unique_LFC_per_residue_human)
            arr_all_edits.append(all_edits_per_residue_human)

        col_head = "mean_missense_LFCr" + str(r+1)
        df_struc_consvr[col_head] = arr_unique_LFC

    edits_filedir = Path(workdir + input_human_gene)
    out_filename = edits_filedir / f"screendata/{input_human_gene}_{structureid}_struc_consrv_missenseedits_randomized.tsv"
    df_struc_consvr.to_csv(out_filename, sep = '\t', index=False)

    return df_struc_consvr
