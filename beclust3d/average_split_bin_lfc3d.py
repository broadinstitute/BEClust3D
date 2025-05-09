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
from _average_split_bin_helpers_ import *
import warnings
warnings.filterwarnings('ignore')

def average_split_bin(
        df_LFC_LFC3D_rand, 
        workdir, input_gene, screen_names, 
        pthr=0.05, score_type='LFC3D', print_statements=True,
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
    
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / score_type):
        os.mkdir(edits_filedir / score_type)

    quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}

    df_bidir = pd.DataFrame()
    df_bidir['unipos'] = df_LFC_LFC3D_rand['unipos']
    df_bidir['unires'] = df_LFC_LFC3D_rand['unires']
    df_LFC_LFC3D_dis = df_bidir[['unipos', 'unires']].copy()
    
    for screen_name in screen_names: # FOR EVERY SCREEN INDIVIDUALLY #
        header_LFC = f"{screen_name}_LFC"
        df_bidir[f"{screen_name}_LFC"] = df_LFC_LFC3D_rand[header_LFC] # LFC per screen

        header_LFC3D = f"{screen_name}_{score_type}"
        if header_LFC3D not in df_LFC_LFC3D_rand.columns: # make sure screen gene combo exists
            warnings.warn(f'{screen_name} screen not found for {input_gene}')
            continue
        df_bidir[header_LFC3D] = df_LFC_LFC3D_rand[header_LFC3D] # LFC or LFC3D per screen
        
        taa_wise_LFC3D_pos, taa_wise_LFC3D_neg = [], []
        
        # CALCULATE THE AVG OF ALL THE RANDOMIZATIONS AND SPLIT POS NEG #
        for aa in range(len(df_LFC_LFC3D_rand)): 
            taa_LFC3D_raw = df_LFC_LFC3D_rand.at[aa, header_LFC3D] # TARGET SCORE per AA per SCREEN 
            
            # SEPARATE INTO POS AND NEG #
            taa_LFC3D = float(taa_LFC3D_raw) if taa_LFC3D_raw != '-' else None
            taa_wise_LFC3D_neg.append(taa_LFC3D if taa_LFC3D is not None and taa_LFC3D < 0 else '-')
            taa_wise_LFC3D_pos.append(taa_LFC3D if taa_LFC3D is not None and taa_LFC3D > 0 else '-')

        df_bidir[f"{screen_name}_{score_type}_neg"]      = taa_wise_LFC3D_neg # LFC3D_neg per SCREEN
        df_bidir[f"{screen_name}_{score_type}_pos"]      = taa_wise_LFC3D_pos # LFC3D_pos per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r"]     = df_LFC_LFC3D_rand[f"{screen_name}_AVG_{score_type}r"] # AVG_LFC3Dr per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r_neg"] = df_LFC_LFC3D_rand[f"{screen_name}_AVG_{score_type}r_neg"] # AVG_LFC3Dr_neg per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r_pos"] = df_LFC_LFC3D_rand[f"{screen_name}_AVG_{score_type}r_pos"] # AVG_LFC3Dr_pos per SCREEN

        # BINNING #
        df_LFC_LFC3D_dis[header_LFC] = df_bidir[header_LFC]
        df_LFC_LFC3D_dis[header_LFC3D] = df_bidir[header_LFC3D]

        # GENERATE THRESHOLDS FOR BINNING #
        df_nodash = df_bidir.loc[df_bidir[header_LFC3D] != '-', ].reset_index(drop=True)
        df_nodash[header_LFC3D] = df_nodash[header_LFC3D].astype(float)
        df_LFC3D_neg = df_nodash.loc[df_nodash[header_LFC3D] < 0, ].reset_index(drop=True)
        df_LFC3D_pos = df_nodash.loc[df_nodash[header_LFC3D] > 0, ].reset_index(drop=True)
        df_neg_stats = df_LFC3D_neg[header_LFC3D].describe()
        df_pos_stats = df_LFC3D_pos[header_LFC3D].describe()

        # CALCULATE BINS #
        quantile_values = {}
        for name, q in quantiles.items():
            quantile_values[name] = df_LFC3D_neg[header_LFC3D].replace('_', np.nan).astype(float).quantile(q)

        arr_disc, arr_weight = binning_neg_pos(df_bidir, df_neg_stats, df_pos_stats, 
                                               quantile_values.values(), header_LFC3D)
        df_LFC_LFC3D_dis[f"{screen_name}_{score_type}_dis"]  = arr_disc
        df_LFC_LFC3D_dis[f"{screen_name}_{score_type}_wght"] = arr_weight

    out_filename_bidir = edits_filedir / f"{score_type}/{input_gene}_{score_type}_bidirectional.tsv"
    df_bidir.to_csv(out_filename_bidir, sep='\t', index=False)
    out_filename_dis = edits_filedir / f"{score_type}/{input_gene}_{score_type}_dis_wght.tsv"
    df_LFC_LFC3D_dis.to_csv(out_filename_dis, sep = '\t', index=False)

    df_z = pd.DataFrame()
    df_z['unipos'] = df_bidir['unipos']
    df_z['unires'] = df_bidir['unires']

    for screen_name in screen_names: # FOR EVERY SCREEN INDIVIDUALLY #

        df_z[f'{screen_name}_{score_type}_neg'] = df_bidir[f'{screen_name}_{score_type}_neg']
        df_z[f'{screen_name}_{score_type}_pos'] = df_bidir[f'{screen_name}_{score_type}_pos']
        df_z[f'{screen_name}_AVG_{score_type}r_neg'] = df_bidir[f'{screen_name}_AVG_{score_type}r_neg']
        df_z[f'{screen_name}_AVG_{score_type}r_pos'] = df_bidir[f'{screen_name}_AVG_{score_type}r_pos']

        # CALCULATE Z SCORE #
        colnames = [f'{screen_name}_{score_type}_{sign}' for sign in ['neg', 'pos']]
        params = [{'mu': df_bidir[f'{screen_name}_AVG_{score_type}r_{sign}'].replace('-', np.nan).astype(float).mean(), 
                    's': df_bidir[f'{screen_name}_AVG_{score_type}r_{sign}'].replace('-', np.nan).astype(float).std()}  
                    ### can we fix the default format 250120
                   for sign in ['neg', 'pos']]
        result_data = {f'{screen_name}_{score_type}_{sign}_{suffix}': [] 
                       for sign in ['neg', 'pos'] for suffix in ['z', 'p', 'psig', '005_psig', '001_psig','0001_psig']}

        for i in range(len(df_z)):
            for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
                if df_z.at[i, colname] != '-' and df_z.at[i, colname] != None:
                    signal = float(df_z.at[i, colname])
                    signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
                    _, _, signal_plabel_005 = calculate_stats(signal, param, 0.05)
                    _, _, signal_plabel_001 = calculate_stats(signal, param, 0.01)
                    _, _, signal_plabel_0001 = calculate_stats(signal, param, 0.001)
                else:
                    signal,signal_z,signal_p, signal_plabel,signal_plabel_005, signal_plabel_001, signal_plabel_0001 = '-','-','-','-','-','-','-'
                # Append results to the dictionary
                result_data[f'{screen_name}_{score_type}_{sign}_z'].append(signal_z)
                result_data[f'{screen_name}_{score_type}_{sign}_p'].append(signal_p)
                result_data[f'{screen_name}_{score_type}_{sign}_psig'].append(signal_plabel)
                result_data[f'{screen_name}_{score_type}_{sign}_005_psig'].append(signal_plabel_005)
                result_data[f'{screen_name}_{score_type}_{sign}_001_psig'].append(signal_plabel_001)
                result_data[f'{screen_name}_{score_type}_{sign}_0001_psig'].append(signal_plabel_0001)

        df_z = pd.concat([df_z, pd.DataFrame(result_data)], axis=1)

    filename = edits_filedir / f"{score_type}/{input_gene}_NonAggr_{score_type}.tsv"
    df_z.to_csv(filename, "\t", index=False)

    return df_bidir, df_LFC_LFC3D_dis, df_z


def average_split_score(
        df_LFC_LFC3D_rand, 
        workdir, input_gene, screen_names, 
        score_type='LFC3D', 
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
    """
    
    edits_filedir = Path(workdir)
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
        df_bidir[header_LFC] = df_LFC_LFC3D_rand[header_LFC] # LFC per screen
        df_bidir[header_LFC3D] = df_LFC_LFC3D_rand[header_LFC3D] # LFC or LFC3D per screen
        
        taa_wise_LFC3D_pos, taa_wise_LFC3D_neg = [], []
        
        # CALCULATE THE AVG OF ALL THE RANDOMIZATIONS AND SPLIT POS NEG #
        taa_LFC3D_raws_dict = df_LFC_LFC3D_rand[header_LFC3D].to_dict()
        for aa in range(len(df_LFC_LFC3D_rand)): 
            taa_LFC3D_raw = taa_LFC3D_raws_dict[aa] # TARGET SCORE per AA per SCREEN 
            
            # SEPARATE INTO POS AND NEG #
            taa_LFC3D = float(taa_LFC3D_raw) if taa_LFC3D_raw != '-' else 0.0
            taa_wise_LFC3D_neg.append(taa_LFC3D if taa_LFC3D < 0 else '-') # either the value or 0.0
            taa_wise_LFC3D_pos.append(taa_LFC3D if taa_LFC3D > 0 else '-') # either the value or 0.0

        df_bidir[f"{screen_name}_{score_type}_neg"]      = taa_wise_LFC3D_neg # LFC3D_neg per SCREEN
        df_bidir[f"{screen_name}_{score_type}_pos"]      = taa_wise_LFC3D_pos # LFC3D_pos per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r"]     = df_LFC_LFC3D_rand[f"{screen_name}_AVG_{score_type}r"] # AVG_LFC3Dr per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r_neg"] = df_LFC_LFC3D_rand[f"{screen_name}_AVG_{score_type}r_neg"] # AVG_LFC3Dr_neg per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r_pos"] = df_LFC_LFC3D_rand[f"{screen_name}_AVG_{score_type}r_pos"] # AVG_LFC3Dr_pos per SCREEN

    out_filename_bidir = edits_filedir / f"{score_type}/{input_gene}_{score_type}_bidirectional.tsv"
    df_bidir.to_csv(out_filename_bidir, sep='\t', index=False)

    return df_bidir

def bin_score(
        df_bidir, 
        workdir, input_gene, screen_names, 
        score_type='LFC3D', 
): 
    """
    Description
        Averages the LFC3D scores, splits into positive and negative, 
        retrieves randomized LFC3D scores, bins LFC 3D scores into percentiles

    Params
        df_bidir: pandas dataframe
            from previous step average_split_lfc3d()
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
    
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / score_type):
        os.mkdir(edits_filedir / score_type)

    quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}

    headers_LFC3D = [f"{sn}_{score_type}" for sn in screen_names]
    df_LFC_LFC3D_dis = df_bidir[['unipos', 'unires'] + headers_LFC3D].copy()
    
    df_neg_stats_list, df_pos_stats_list = [], []
    for screen_name in screen_names: # FOR EVERY SCREEN INDIVIDUALLY #
        header_LFC3D = f"{screen_name}_{score_type}"

        # BINNING #
        df_LFC_LFC3D_dis[header_LFC3D] = df_bidir[header_LFC3D]

        # GENERATE THRESHOLDS FOR BINNING #
        # df_nodash = df_bidir.loc[df_bidir[header_LFC3D] != '-', ].reset_index(drop=True)
        # df_nodash[header_LFC3D] = df_nodash[header_LFC3D].astype(float)
        # df_LFC3D_neg = df_nodash.loc[df_nodash[header_LFC3D] < 0, ].reset_index(drop=True)
        # df_LFC3D_pos = df_nodash.loc[df_nodash[header_LFC3D] > 0, ].reset_index(drop=True)
        # df_neg_stats = df_LFC3D_neg[header_LFC3D].describe()
        # df_pos_stats = df_LFC3D_pos[header_LFC3D].describe()
        df_temp = df_LFC_LFC3D_dis[header_LFC3D].replace('-', np.nan).astype(float)
        mask_neg = df_temp < 0.0
        mask_pos = df_temp > 0.0
        df_neg_stats = df_temp[mask_neg].describe()
        df_pos_stats = df_temp[mask_pos].describe()
        df_neg_stats_list.append(df_neg_stats)
        df_pos_stats_list.append(df_pos_stats)

        # CALCULATE BINS #
        quantile_values = {}
        for name, q in quantiles.items(): 
            quantile_values[name] = df_LFC_LFC3D_dis[header_LFC3D].replace('-', np.nan).astype(float).quantile(q)

        arr_disc, arr_weight = binning_neg_pos(df_bidir, df_neg_stats, df_pos_stats, 
                                               quantile_values.values(), header_LFC3D)
        df_LFC_LFC3D_dis[f"{screen_name}_{score_type}_dis"]  = arr_disc
        df_LFC_LFC3D_dis[f"{screen_name}_{score_type}_wght"] = arr_weight

    out_filename_dis = edits_filedir / f"{score_type}/{input_gene}_{score_type}_dis_wght.tsv"
    df_LFC_LFC3D_dis.to_csv(out_filename_dis, sep = '\t', index=False)

    return df_LFC_LFC3D_dis, df_neg_stats_list, df_pos_stats_list


def znorm_score(
        df_bidir, df_neg_stats_list, df_pos_stats_list, 
        workdir, input_gene, screen_names, 
        pthrs=[0.05, 0.01, 0.001], score_type='LFC3D', 
): 
    """
    Description
        Averages the LFC3D scores, splits into positive and negative, 
        retrieves randomized LFC3D scores, bins LFC 3D scores into percentiles

    Params
        df_bidir: pandas dataframe
            from previous step average_split_lfc3d()
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
    
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / score_type):
        os.mkdir(edits_filedir / score_type)
    
    df_z = pd.DataFrame()
    df_z['unipos'] = df_bidir['unipos']
    df_z['unires'] = df_bidir['unires']

    for sn, df_neg, df_pos in zip(screen_names, df_neg_stats_list, df_pos_stats_list): # FOR EVERY SCREEN INDIVIDUALLY #

        df_z[f'{sn}_{score_type}_neg'] = df_bidir[f'{sn}_{score_type}_neg']
        df_z[f'{sn}_{score_type}_pos'] = df_bidir[f'{sn}_{score_type}_pos']
        df_z[f'{sn}_AVG_{score_type}r_neg'] = df_bidir[f'{sn}_AVG_{score_type}r_neg']
        df_z[f'{sn}_AVG_{score_type}r_pos'] = df_bidir[f'{sn}_AVG_{score_type}r_pos']

        header_main = f'{sn}_{score_type}'
        pthrs_str = [str(pthr).split('.')[1] for pthr in pthrs]
        # CALCULATE Z SCORE #
        colnames = [f'{sn}_{score_type}_{sign}' for sign in ['neg', 'pos']]
        params = [{'mu': df_neg['mean'], 's': df_neg['std']}, 
                  {'mu': df_pos['mean'], 's': df_pos['std']}, ]
        # params = [{'mu': df_bidir[f'{sn}_AVG_{score_type}r_{sign}'].replace('-', np.nan).astype(float).mean(), 
        #             's': df_bidir[f'{sn}_AVG_{score_type}r_{sign}'].replace('-', np.nan).astype(float).std()}  
        #            for sign in ['neg', 'pos']]

        result_data = {f'{header_main}_{sign}_{pthr_str}_{suffix}': [] 
                       for sign in ['neg', 'pos'] for suffix in ['z', 'p', 'psig'] for pthr_str in pthrs_str}

        for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
            signals_dict = df_z[colname].replace('-', np.nan).to_dict()

            for pthr, pthr_str in zip(pthrs, pthrs_str): 
                for i in range(len(df_z)):
                    signal = float(signals_dict[i])
                    signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
                    
                    # Append results to the dictionary
                    result_data[f'{sn}_{score_type}_{sign}_{pthr_str}_z'].append(signal_z)
                    result_data[f'{sn}_{score_type}_{sign}_{pthr_str}_p'].append(signal_p)
                    result_data[f'{sn}_{score_type}_{sign}_{pthr_str}_psig'].append(signal_plabel)

        df_z = pd.concat([df_z, pd.DataFrame(result_data)], axis=1)

    filename = edits_filedir / f"{score_type}/{input_gene}_NonAggr_{score_type}.tsv"
    df_z.to_csv(filename, "\t", index=False)

    return df_z
