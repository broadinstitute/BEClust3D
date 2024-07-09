"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2.5

"""

import pandas as pd
from pathlib import Path
import os
from functools import lru_cache

def randomize_by_conservation(
        workdir, 
        input_gene, input_screen, structureid, 
        nRandom=1000, 
): 

    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')

    edits_filename = edits_filedir / f"screendata/{input_gene}_{structureid}_struc_consrv.tsv"
    df_struc_consvr = pd.read_csv(edits_filename, sep = "\t")
    
    # for Missense load output from prioritize_by_conservation
    missense_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screen_name}_missense_edits_randomized.tsv"
    df_missense = pd.read_csv(missense_filename, sep = '\t')

    arr_unique_LFC = [0.0] * len(df_struc_consvr)
    arr_all_edits = [0.0] * len(df_struc_consvr)
    df_pos_edits_list = []

    for i in range(0, len(df_struc_consvr)):
        human_res_pos = df_struc_consvr.at[i, 'human_res_pos']
        header_LFC_rN = "LFC"
        unique_LFC, all_edits, df_pos_edits = calc_lfc_edits_helper(df_missense, human_res_pos, header_LFC_rN)
        arr_unique_LFC[i], arr_all_edits[i] = unique_LFC, all_edits
        df_pos_edits_list.append(df_pos_edits)

    df_struc_consvr['all_missense_edits'] = arr_all_edits
    df_struc_consvr['mean_missense_LFC'] = arr_unique_LFC
    
    dict_temp = {}
    for r in range(0, nRandom):
        arr_unique_LFC = [0.0] * len(df_struc_consvr)
        # arr_all_edits = [0.0] * len(df_struc_consvr)
        header_LFC_rN = f"LFCr{str(r+1)}"
        header = f"mean_missense_LFCr{str(r+1)}"

        for i in range(0, len(df_struc_consvr)):
            human_res_pos = df_struc_consvr.at[i, 'human_res_pos']
            unique_LFC = calc_lfc_helper(header_LFC_rN, df_pos_edits_list[i])
            arr_unique_LFC[i] = unique_LFC
        dict_temp[header] = arr_unique_LFC

    df_struc_consvr = pd.concat((df_struc_consvr, pd.DataFrame(dict_temp)), axis=1)
    out_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screen_name}_struc_consrv_missenseedits_randomized.tsv"
    df_struc_consvr.to_csv(out_filename, sep = '\t', index=False)

    return df_struc_consvr

def calc_lfc_edits_helper(df_missense, human_res_pos, header_LFC_rN): 
    unique_LFC_per_residue_human, all_edits_per_residue_human = '-', '-'
    df_pos_edits = df_missense.loc[df_missense['human_pos'] == int(human_res_pos), ].reset_index()

    if len(df_pos_edits) > 1:
        # edit
        pos_edits_list = df_pos_edits['edit'].tolist()
        all_edits_per_residue_human = ';'.join(list(set(pos_edits_list)))
        # score
        pos_LFCscore_list = df_pos_edits[header_LFC_rN].tolist()
        unique_LFC_per_residue_human = str_round_div(sum(pos_LFCscore_list), len(pos_LFCscore_list))
    elif len(df_pos_edits) == 1:
        all_edits_per_residue_human = df_pos_edits.at[0,'edit']
        unique_LFC_per_residue_human = str_round(df_pos_edits.at[0, header_LFC_rN])
        
    return unique_LFC_per_residue_human, all_edits_per_residue_human, df_pos_edits

def calc_lfc_helper(header_LFC_rN, df_pos_edits): 
    unique_LFC_per_residue_human = '-'

    if len(df_pos_edits) > 1:
        pos_LFCscore_list = df_pos_edits[header_LFC_rN].tolist()
        unique_LFC_per_residue_human = str_round_div(sum(pos_LFCscore_list), len(pos_LFCscore_list))
    elif len(df_pos_edits) == 1:
        unique_LFC_per_residue_human = str_round(df_pos_edits.at[0, header_LFC_rN])
        
    return unique_LFC_per_residue_human

@lru_cache
def str_round_div(x, y): 
    return str(round(x/y, 3))

@lru_cache
def str_round(x): 
    return str(round(x, 3))
