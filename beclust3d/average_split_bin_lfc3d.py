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
from _average_split_bin_plots_ import *

def average_split_bin(
        df_LFC_LFC3D_rand, 
        workdir, input_gene, structureid, input_screens, screen_names=[], 
        nRandom=1000, pthr=0.05, score_type='LFC3D', 
): 
    """
    Description
        Averages the LFC3D scores, splits into positive and negative, 
        retrieves randomized LFC3D scores, bins LFC 3D scores into percentiles

    Params
        df_LFC_LFC3D_rand: pandas dataframe
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
            a dataframe listing the positive and negative components of df_LFC_LFC3D_rand
        df_LFC_LFC3D_dis: pandas dataframe
            a dataframe listing how scores portion into quantiles
    """
    
    edits_filedir = Path(workdir + '/' + input_gene)
    if not screen_names: 
        screen_names = [input_screen.split('.')[0] for input_screen in input_screens]
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / score_type):
        os.mkdir(edits_filedir / score_type)

    df_bidir = pd.DataFrame()
    df_bidir['unipos'] = df_LFC_LFC3D_rand['unipos']
    df_bidir['unires'] = df_LFC_LFC3D_rand['unires']
    
    for screen_name in screen_names: # FOR EVERY SCREEN INDIVIDUALLY #

        header_LFC = f"{screen_name}_LFC"
        header_LFC3D = f"{screen_name}_{score_type}"
        if header_LFC not in df_LFC_LFC3D_rand.columns: # make sure screen gene combo exists
            warnings.warn(f'{screen_name} screen not found for {input_gene}')
            continue

        df_bidir[header_LFC] = df_LFC_LFC3D_rand[header_LFC] # LFC/screen
        df_bidir[header_LFC3D] = df_LFC_LFC3D_rand[header_LFC3D] # LFC3D/screen
        
        taa_wise_LFC3D_pos, taa_wise_LFC3D_neg = [], []
        taa_wise_AVG_LFC3Dr_pos, taa_wise_AVG_LFC3Dr_neg, taa_wise_AVG_LFC3Dr = [], [], []
        
        # CALCULATE THE AVG OF ALL THE RANDOMIZATIONS #
        for aa in range(len(df_LFC_LFC3D_rand)): # for each residue calculate positive and negative signals
            taa_LFC3D_raw = df_LFC_LFC3D_rand.at[aa, header_LFC3D] # target LFC3D/aa/screen 
            taa_LFC3D_pos, taa_LFC3D_neg = 0.0, 0.0
            
            if taa_LFC3D_raw != '-':
                taa_LFC3D = float(taa_LFC3D_raw)
                if taa_LFC3D < 0.0: taa_LFC3D_neg = taa_LFC3D
                if taa_LFC3D > 0.0: taa_LFC3D_pos = taa_LFC3D
            taa_wise_LFC3D_neg.append(taa_LFC3D_neg) # either the value or 0.0
            taa_wise_LFC3D_pos.append(taa_LFC3D_pos) # either the value or 0.0
            
            # For randomized LFC3Dr - find AVG(neg_LFC3Dr) and AVG(pos_LFC3Dr)
            taa_SUM_LFC3Dr, taa_SUM_LFC3Dr_neg, taa_SUM_LFC3Dr_pos = 0.0, 0.0, 0.0
            taas_LFC3Dr_raw = df_LFC_LFC3D_rand.loc[aa, [f'{screen_name}_{score_type}r{str(r+1)}' for r in range(nRandom)]]

            for taa_LFC3Dr_raw in taas_LFC3Dr_raw: 
                taa_LFC3Dr_pos, taa_LFC3Dr_neg = 0.0, 0.0
                if taa_LFC3Dr_raw == '-': 
                    break
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

        df_bidir[f"{screen_name}_{score_type}_neg"]      = taa_wise_LFC3D_neg ## LFC3D_neg/s
        df_bidir[f"{screen_name}_{score_type}_pos"]      = taa_wise_LFC3D_pos ## LFC3D_pos/s
        df_bidir[f"{screen_name}_AVG_{score_type}r"]     = taa_wise_AVG_LFC3Dr ## AVG_LFC3Dr/s
        df_bidir[f"{screen_name}_AVG_{score_type}r_neg"] = taa_wise_AVG_LFC3Dr_neg ## AVG_LFC3Dr_neg/s
        df_bidir[f"{screen_name}_AVG_{score_type}r_pos"] = taa_wise_AVG_LFC3Dr_pos ## AVG_LFC3Dr_pos/s   
    
        # BINNING #
        df_LFC_LFC3D_dis = df_bidir[['unipos', 'unires', header_LFC, header_LFC3D]].copy()

        # generate bin thresholds
        df_nodash = df_bidir.loc[df_bidir[header_LFC3D] != '-', ]
        df_nodash = df_nodash.reset_index(drop=True)
        df_nodash[header_LFC3D] = df_nodash[header_LFC3D].astype(float)

        df_LFC3D_neg = df_nodash.loc[df_nodash[header_LFC3D] < 0, ].reset_index(drop=True)
        df_LFC3D_pos = df_nodash.loc[df_nodash[header_LFC3D] > 0, ].reset_index(drop=True)
        df_neg_stats = df_LFC3D_neg[header_LFC3D].describe()
        df_pos_stats = df_LFC3D_pos[header_LFC3D].describe()

        quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}
        quantile_values = {}
        for name, q in quantiles.items(): 
            quantile_values[name] = round(df_LFC3D_neg[header_LFC3D].quantile(q), 4)

        arr_disc, arr_weight = binning(df_bidir, df_neg_stats, df_pos_stats, quantile_values.values(), header_LFC3D)
        df_LFC_LFC3D_dis[f"{screen_name}_{score_type}_dis"]  = arr_disc
        df_LFC_LFC3D_dis[f"{screen_name}_{score_type}_wght"] = arr_weight

    out_filename = edits_filedir / f"{score_type}/{input_gene}_{score_type}_bidirectional.tsv"
    df_bidir.to_csv(out_filename, sep='\t', index=False)
    out_filename = edits_filedir / f"{score_type}/{input_gene}_{score_type}_dis_wght.tsv"
    df_LFC_LFC3D_dis.to_csv(out_filename, sep = '\t', index=False)

    df_z = pd.DataFrame()
    for screen_name in screen_names: # FOR EVERY SCREEN INDIVIDUALLY #

        df_z['unipos'] = df_bidir['unipos']
        df_z['unires'] = df_bidir['unires']
        df_z[f'{screen_name}_SUM_{score_type}_neg'] = df_bidir[f'{screen_name}_{score_type}_neg']
        df_z[f'{screen_name}_SUM_{score_type}_pos'] = df_bidir[f'{screen_name}_{score_type}_pos']
        df_z[f'{screen_name}_AVG_{score_type}r_neg'] = df_bidir[f'{screen_name}_AVG_{score_type}r_neg']
        df_z[f'{screen_name}_AVG_{score_type}r_pos'] = df_bidir[f'{screen_name}_AVG_{score_type}r_pos']

        # CALCULATE Z SCORE #
        colnames = [f'{screen_name}_SUM_{score_type}_{sign}' for sign in ['neg', 'pos']]
        params = [{'mu': df_bidir[f'{screen_name}_AVG_{score_type}r_{sign}'].mean(),
                   's': df_bidir[f'{screen_name}_AVG_{score_type}r_{sign}'].std()} 
                   for sign in ['neg', 'pos']]
        result_data = {f'{screen_name}_SUM_{score_type}_{sign}_{suffix}': [] 
                       for sign in ['neg', 'pos'] for suffix in ['z', 'p', 'psig']}

        for i in range(len(df_z)):
            for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
                signal = float(df_z.at[i, colname])
                signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
                
                # Append results to the dictionary
                result_data[f'{screen_name}_SUM_{score_type}_{sign}_z'].append(signal_z)
                result_data[f'{screen_name}_SUM_{score_type}_{sign}_p'].append(signal_p)
                result_data[f'{screen_name}_SUM_{score_type}_{sign}_psig'].append(signal_plabel)

        df_z = pd.concat([df_z, pd.DataFrame(result_data)], axis=1).round(4)
        # PLOTS #
        LFC3D_plots(df_z, edits_filedir, input_gene, pthr, screen_name+'_', score=score_type, )

    filename = edits_filedir / f"{score_type}/{structureid}_NonAggr_{score_type}.tsv"
    df_z.to_csv(filename, "\t", index=False)

    return df_bidir, df_LFC_LFC3D_dis, df_z
