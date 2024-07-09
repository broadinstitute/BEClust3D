"""
File: average_split_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.3.5

"""

import pandas as pd
from pathlib import Path

def average_and_split(
        df_LFC_LFCrN_LFC3D_LFC3DrN, 
        workdir, input_gene,
        input_screen, 
        nRandom=1000, 
): 
    
    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir + '/' + input_gene)
    df_bidir = pd.DataFrame()
    df_bidir['unipos'] = df_LFC_LFCrN_LFC3D_LFC3DrN['unipos']
    df_bidir['unires'] = df_LFC_LFCrN_LFC3D_LFC3DrN['unires']
    
    header_LFC = f"{screen_name}_LFC"
    header_LFC3D = f"{screen_name}_LFC3D"
    df_bidir[header_LFC] = df_LFC_LFCrN_LFC3D_LFC3DrN[header_LFC] # LFC/screen
    df_bidir[header_LFC3D] = df_LFC_LFCrN_LFC3D_LFC3DrN[header_LFC3D] # LFC3D/screen
    
    taa_wise_LFC3D_pos, taa_wise_LFC3D_neg = [], []
    taa_wise_AVG_LFC3Dr_pos, taa_wise_AVG_LFC3Dr_neg, taa_wise_AVG_LFC3Dr = [], [], []
    
    for aa in range(0, len(df_LFC_LFCrN_LFC3D_LFC3DrN)):
        taa_LFC3D_raw = df_LFC_LFCrN_LFC3D_LFC3DrN.at[aa, header_LFC3D] # target LFC3D/aa/screen 
        taa_LFC3D_pos, taa_LFC3D_neg = 0.0, 0.0
        
        if taa_LFC3D_raw != '-':
            taa_LFC3D = float(taa_LFC3D_raw)
            if taa_LFC3D < 0.0: taa_LFC3D_neg = taa_LFC3D
            if taa_LFC3D > 0.0: taa_LFC3D_pos = taa_LFC3D
        taa_wise_LFC3D_neg.append(taa_LFC3D_neg) # either the value or 0.0
        taa_wise_LFC3D_pos.append(taa_LFC3D_pos) # either the value or 0.0
        
        # For randomized LFC3Dr - find AVG(neg_LFC3Dr) and AVG(pos_LFC3Dr) 
        taa_SUM_LFC3Dr, taa_SUM_LFC3Dr_neg, taa_SUM_LFC3Dr_pos = 0.0, 0.0, 0.0
        taa_AVG_LFC3Dr, taa_AVG_LFC3Dr_neg, taa_AVG_LFC3Dr_pos = 0.0, 0.0, 0.0
        
        for r in range(0, nRandom):
            col_head =  f'{screen_name}_LFC3D{str(r+1)}'
            taa_LFC3Dr_raw = df_LFC_LFCrN_LFC3D_LFC3DrN.at[aa, col_head] # target LFC3D (randomized)
            
            taa_LFC3Dr_pos, taa_LFC3Dr_neg = 0.0, 0.0
            if taa_LFC3Dr_raw == '-':
                taa_LFC3Dr = 0.0
            else:
                taa_LFC3Dr = float(taa_LFC3Dr_raw)
                if taa_LFC3Dr < 0.0: taa_LFC3Dr_neg = taa_LFC3Dr
                if taa_LFC3Dr > 0.0: taa_LFC3Dr_pos = taa_LFC3Dr
                    
            taa_SUM_LFC3Dr     = taa_SUM_LFC3Dr + taa_LFC3Dr
            taa_SUM_LFC3Dr_neg = taa_SUM_LFC3Dr_neg + taa_LFC3Dr_neg
            taa_SUM_LFC3Dr_pos = taa_SUM_LFC3Dr_pos + taa_LFC3Dr_pos
        
        taa_AVG_LFC3Dr_neg = round(taa_SUM_LFC3Dr_neg / nRandom, 3)
        taa_wise_AVG_LFC3Dr_neg.append(taa_AVG_LFC3Dr_neg)
            
        taa_AVG_LFC3Dr_pos = round(taa_SUM_LFC3Dr_pos / nRandom, 3)        
        taa_wise_AVG_LFC3Dr_pos.append(taa_AVG_LFC3Dr_pos)

        taa_AVG_LFC3Dr = round(taa_SUM_LFC3Dr / nRandom, 3)
        taa_wise_AVG_LFC3Dr.append(taa_AVG_LFC3Dr)

    df_bidir[f"{screen_name}_LFC3D_neg"]      = taa_wise_LFC3D_neg ## LFC3D_neg/s
    df_bidir[f"{screen_name}_LFC3D_pos"]      = taa_wise_LFC3D_pos ## LFC3D_pos/s
    df_bidir[f"{screen_name}_AVG_LFC3Dr"]     = taa_wise_AVG_LFC3Dr ## AVG_LFC3Dr/s
    df_bidir[f"{screen_name}_AVG_LFC3Dr_neg"] = taa_wise_AVG_LFC3Dr_neg ## AVG_LFC3Dr_neg/s
    df_bidir[f"{screen_name}_AVG_LFC3Dr_pos"] = taa_wise_AVG_LFC3Dr_pos ## AVG_LFC3Dr_pos/s   
    
    out_filename = edits_filedir / f"LFC3D/{input_gene}_{screen_name}_LFC_LFC3D_LFC3Dr_bidirectional.tsv"
    df_bidir.to_csv(out_filename, 
                    sep='\t', index=False)

    return df_bidir
