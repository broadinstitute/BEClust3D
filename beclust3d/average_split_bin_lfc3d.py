"""
File: average_split_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.3.5 and 3.4

"""

import pandas as pd
from pathlib import Path
import warnings
import os

def average_and_split(
        df_LFC_LFCrN_LFC3D_LFC3DrN, 
        workdir, 
        input_gene, input_screens, 
        nRandom=1000, 
): 
    """
    Description
        Averages the LFC3D scores, splits into positive and negative, 
        retrieves randomized LFC3D scores, bins LFC 3D scores into percentiles

    Params
        df_LFC_LFCrN_LFC3D_LFC3DrN: pandas dataframe
            from previous step calculate_lfc3d()
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        nRandom: int, optional
            the number of randomize iterations

    Returns
        df_bidir: pandas dataframe
            a dataframe listing the positive and negative components of df_LFC_LFCrN_LFC3D_LFC3DrN
        df_LFC_LFC3D_dis: pandas dataframe
            a dataframe listing how scores portion into quantiles
    """
    
    edits_filedir = Path(workdir + '/' + input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'LFC3D'):
        os.mkdir(edits_filedir / 'LFC3D')

    df_bidir = pd.DataFrame()
    df_bidir['unipos'] = df_LFC_LFCrN_LFC3D_LFC3DrN['unipos']
    df_bidir['unires'] = df_LFC_LFCrN_LFC3D_LFC3DrN['unires']
    
    for input_screen in input_screens: # for every screen individually

        screen_name = input_screen.split('.')[0]
        header_LFC = f"{screen_name}_LFC"
        header_LFC3D = f"{screen_name}_LFC3D"
        if header_LFC not in df_LFC_LFCrN_LFC3D_LFC3DrN.columns: # make sure screen gene combo exists
            warnings.warn(f'{screen_name} screen not found for {input_gene}')
            continue

        df_bidir[header_LFC] = df_LFC_LFCrN_LFC3D_LFC3DrN[header_LFC] # LFC/screen
        df_bidir[header_LFC3D] = df_LFC_LFCrN_LFC3D_LFC3DrN[header_LFC3D] # LFC3D/screen
        
        taa_wise_LFC3D_pos, taa_wise_LFC3D_neg = [], []
        taa_wise_AVG_LFC3Dr_pos, taa_wise_AVG_LFC3Dr_neg, taa_wise_AVG_LFC3Dr = [], [], []
        
        for aa in range(0, len(df_LFC_LFCrN_LFC3D_LFC3DrN)): # for each residue calculate positive and negative signals
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

            for r in range(0, nRandom):
                col_head =  f'{screen_name}_LFC3Dr{str(r+1)}'
                taa_LFC3Dr_raw = df_LFC_LFCrN_LFC3D_LFC3DrN.at[aa, col_head] # target LFC3D (randomized)
                
                taa_LFC3Dr_pos, taa_LFC3Dr_neg = 0.0, 0.0
                if taa_LFC3Dr_raw == '-': 
                    taa_LFC3Dr = 0.0
                else: 
                    taa_LFC3Dr = float(taa_LFC3Dr_raw)
                    if taa_LFC3Dr < 0.0: taa_LFC3Dr_neg = taa_LFC3Dr
                    if taa_LFC3Dr > 0.0: taa_LFC3Dr_pos = taa_LFC3Dr
                        
                taa_SUM_LFC3Dr     += taa_LFC3Dr
                taa_SUM_LFC3Dr_neg += taa_LFC3Dr_neg
                taa_SUM_LFC3Dr_pos += taa_LFC3Dr_pos
            
            taa_wise_AVG_LFC3Dr_neg.append(round(taa_SUM_LFC3Dr_neg/nRandom, 3))
            taa_wise_AVG_LFC3Dr_pos.append(round(taa_SUM_LFC3Dr_pos/nRandom, 3))
            taa_wise_AVG_LFC3Dr.append(round(taa_SUM_LFC3Dr/nRandom, 3))

        df_bidir[f"{screen_name}_LFC3D_neg"]      = taa_wise_LFC3D_neg ## LFC3D_neg/s
        df_bidir[f"{screen_name}_LFC3D_pos"]      = taa_wise_LFC3D_pos ## LFC3D_pos/s
        df_bidir[f"{screen_name}_AVG_LFC3Dr"]     = taa_wise_AVG_LFC3Dr ## AVG_LFC3Dr/s
        df_bidir[f"{screen_name}_AVG_LFC3Dr_neg"] = taa_wise_AVG_LFC3Dr_neg ## AVG_LFC3Dr_neg/s
        df_bidir[f"{screen_name}_AVG_LFC3Dr_pos"] = taa_wise_AVG_LFC3Dr_pos ## AVG_LFC3Dr_pos/s   
    
        # BINNING #
        df_LFC_LFC3D_dis = df_bidir[['unipos', 'unires', header_LFC, header_LFC3D]].copy()

        # generate bin thresholds
        df_LFC_LFC3D_nodash = df_bidir.loc[df_bidir[header_LFC3D] != '-', ]
        df_LFC_LFC3D_nodash = df_LFC_LFC3D_nodash.reset_index(drop=True)
        df_LFC_LFC3D_nodash[header_LFC3D] = df_LFC_LFC3D_nodash[header_LFC3D].astype(float)

        df_LFC3D_neg = df_LFC_LFC3D_nodash.loc[df_LFC_LFC3D_nodash[header_LFC3D] < 0, ]
        df_LFC3D_neg = df_LFC3D_neg.reset_index(drop=True)
        print("length of neg: " + str(len(df_LFC3D_neg)))
        df_LFC3D_pos = df_LFC_LFC3D_nodash.loc[df_LFC_LFC3D_nodash[header_LFC3D] > 0, ]
        df_LFC3D_pos = df_LFC3D_pos.reset_index(drop=True)
        print("length of pos: " + str(len(df_LFC3D_pos)))

        df_LFC3D_neg_stats = df_LFC3D_neg[header_LFC3D].describe()
        df_LFC3D_pos_stats = df_LFC3D_pos[header_LFC3D].describe()

        quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}
        quantile_values = {}
        for name, q in quantiles.items(): 
            quantile_values[name] = round(df_LFC3D_neg[header_LFC3D].quantile(q), 4)
            print(name, q)

        arr_LFC3D_discrete, arr_LFC3D_weight = binning(df_bidir, df_LFC3D_neg_stats, df_LFC3D_pos_stats, 
                                                       quantile_values.values(), header_LFC3D )
        df_LFC_LFC3D_dis[f"{screen_name}_LFC3D_dis"]  = arr_LFC3D_discrete
        df_LFC_LFC3D_dis[f"{screen_name}_LFC3D_wght"] = arr_LFC3D_weight

    out_filename = edits_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_LFC3Dr_bidirectional.tsv"
    df_bidir.to_csv(out_filename, sep='\t', index=False)
    out_filename = edits_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_dis_wght_Signal_Only_per_Screen.tsv"
    df_LFC_LFC3D_dis.to_csv(out_filename, sep = '\t', index=False)

    return df_bidir, df_LFC_LFC3D_dis


def binning(
        df_LFC_LFC3D, df_LFC3D_neg_stats, df_LFC3D_pos_stats, 
        quantile_vals, LFC3D_header
):
    NEG_10p_v, POS_90p_v, NEG_05p_v, POS_95p_v = quantile_vals
    # binning and weighting 
    arr_LFC3D_discrete = []
    arr_LFC3D_weight = []

    for i in range(0, len(df_LFC_LFC3D)): 
        LFC3D = df_LFC_LFC3D.at[i, LFC3D_header]
        if LFC3D == '-' or LFC3D == 0.0:
            LFC3D_discrete = '-'
            LFC3D_weight = 0.0
        else: 
            LFC3Df = round(float(LFC3D), 3)
            # aligned for better readability
            if                               LFC3Df <= NEG_05p_v:                 LFC3D_discrete, LFC3D_weight = 'NEG_05p', -0.95
            elif                 NEG_05p_v < LFC3Df <= NEG_10p_v:                 LFC3D_discrete, LFC3D_weight = 'NEG_10p', -0.9
            elif                 NEG_10p_v < LFC3Df <= df_LFC3D_neg_stats['25%']: LFC3D_discrete, LFC3D_weight = 'NEG_25p', -0.75
            elif df_LFC3D_neg_stats['25%'] < LFC3Df <= df_LFC3D_neg_stats['50%']: LFC3D_discrete, LFC3D_weight = 'NEG_50p', -0.5
            elif df_LFC3D_neg_stats['50%'] < LFC3Df <= df_LFC3D_neg_stats['75%']: LFC3D_discrete, LFC3D_weight = 'NEG_75p', -0.25
            elif df_LFC3D_neg_stats['75%'] < LFC3Df <= df_LFC3D_neg_stats['max']: LFC3D_discrete, LFC3D_weight = 'NEG_100p', -0.05
            elif df_LFC3D_pos_stats['25%'] > LFC3Df >= df_LFC3D_pos_stats['min']: LFC3D_discrete, LFC3D_weight = 'POS_0p', 0.05
            elif df_LFC3D_pos_stats['50%'] > LFC3Df >= df_LFC3D_pos_stats['25%']: LFC3D_discrete, LFC3D_weight = 'POS_25p', 0.25
            elif df_LFC3D_pos_stats['75%'] > LFC3Df >= df_LFC3D_pos_stats['50%']: LFC3D_discrete, LFC3D_weight = 'POS_50p', 0.50
            elif                 POS_90p_v > LFC3Df >= df_LFC3D_pos_stats['75%']: LFC3D_discrete, LFC3D_weight = 'POS_75p', 0.75
            elif                 POS_95p_v > LFC3Df >= POS_90p_v:                 LFC3D_discrete, LFC3D_weight = 'POS_90p', 0.90
            elif                             LFC3Df >= POS_95p_v:                 LFC3D_discrete, LFC3D_weight = 'POS_95p', 0.95

        arr_LFC3D_discrete.append(LFC3D_discrete)
        arr_LFC3D_weight.append(LFC3D_weight)

    return arr_LFC3D_discrete, arr_LFC3D_weight
