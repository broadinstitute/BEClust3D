"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 4

"""

import pandas as pd
import statistics
import scipy.stats as stats
from pathlib import Path
import os
from _average_split_bin_plots_ import *
import warnings

def metaaggregation(
    df_LFC_LFC3D, workdir, 
    input_gene, structureid, input_screens, screen_names=[], 
    nRandom=1000, pthr=0.05, score_type='LFC3D', aggr_func=np.sum, 
): 
    """
    Description
        A point to meta aggregate across multiple screens or just one screen, 
        calculate signal vs background, bin this new metaaggregated signal 
        into top and bottom 10 %, and plot QC graphics

    Params
        df_LFC_LFC3D: pandas dataframe
            from previous step average_split_lfc3d()
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        structureid: str, required
            the name of the AF and uniprot input
        nRandom: int, optional
            the number of randomize iterations
        score_type: str, optional
            'LFC' or 'LFC3D'
        pthr: float, optional
            the p-value threshold
        aggr_func: function, optional
            the function to apply ie np.sum np.min, np.max, np.median, np.mean

    Returns
        df_meta: pandas DataFrame
            results of meta aggregated screen data
    """
    
    edits_filedir = Path(workdir + '/' + input_gene)
    if not screen_names: 
        screen_names = [input_screen.split('.')[0] for input_screen in input_screens]
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'metaaggregation'):
        os.mkdir(edits_filedir / 'metaaggregation')
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    # AGGREGATE LFC3D #
    df_meta = pd.DataFrame()
    df_meta['unipos'] = df_LFC_LFC3D['unipos']
    df_meta['unires'] = df_LFC_LFC3D['unires']

    # SUM LFC3D VALUES ACROSS SCREENS #
    list_LFC3D_neg, list_LFC3D_pos = [], []
    for i in range(len(df_LFC_LFC3D)): 
        values_LFC3D_neg, values_LFC3D_pos = [], []

        for screen_name in screen_names: 
            header_LFC3D = f"{screen_name}_{score_type}"
            if header_LFC3D not in df_LFC_LFC3D.columns: 
                warnings.warn(f'{header_LFC3D} not in input df_LFC_LFC3D')
                continue

            LFC3D = df_LFC_LFC3D.at[i, header_LFC3D]
            if LFC3D != '-':
                LFC3D_value = float(LFC3D)
                if LFC3D_value < 0.0:   values_LFC3D_neg.append(LFC3D_value)
                elif LFC3D_value > 0.0: values_LFC3D_pos.append(LFC3D_value)

        # APPLY AGGR FUNCTION #
        list_LFC3D_neg.append(aggr_func(values_LFC3D_neg) if values_LFC3D_neg else 0.0)
        list_LFC3D_pos.append(aggr_func(values_LFC3D_pos) if values_LFC3D_pos else 0.0)

    df_meta[f'{aggr_func.__name__.upper()}_{score_type}_neg'] = list_LFC3D_neg
    df_meta[f'{aggr_func.__name__.upper()}_{score_type}_pos'] = list_LFC3D_pos
    del list_LFC3D_neg, list_LFC3D_pos

    # RANDOMIZE #
    dict_META = {}
    for r in range(nRandom):
        list_sum_LFC3D_neg, list_sum_LFC3D_pos = [], []
        for i in range(len(df_LFC_LFC3D)): 
            res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0

            for screen_name in screen_names: # sum across screens
                header_LFC3D = f"{screen_name}_{score_type}r{str(r+1)}"
                if header_LFC3D not in df_LFC_LFC3D.columns: 
                    continue

                LFC3D = df_LFC_LFC3D.at[i, header_LFC3D]
                if (LFC3D != '-') and (float(LFC3D) < 0.0):
                    res_sum_LFC3D_neg += float(LFC3D)
                if (LFC3D != '-') and (float(LFC3D) > 0.0):
                    res_sum_LFC3D_pos += float(LFC3D)

            list_sum_LFC3D_neg.append(res_sum_LFC3D_neg)
            list_sum_LFC3D_pos.append(res_sum_LFC3D_pos)
            
        dict_META[f'SUM_{score_type}r{str(r+1)}_neg'] = list_sum_LFC3D_neg
        dict_META[f'SUM_{score_type}r{str(r+1)}_pos'] = list_sum_LFC3D_pos
        del list_sum_LFC3D_neg, list_sum_LFC3D_pos

    df_meta = pd.concat([df_meta, pd.DataFrame(dict_META)], axis=1)
    del dict_META

    # AVG SUM OF RANDOMIZED SIGNAL FOR BACKGROUND #
    list_avg_LFC3D_neg, list_avg_LFC3D_pos = [], []
    for i in range(len(df_meta)): 
        res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0
        for r in range(nRandom): 
            LFC3D_neg = df_meta.at[i, f'SUM_{score_type}r{str(r+1)}_neg']
            res_sum_LFC3D_neg += float(LFC3D_neg)
            LFC3D_pos = df_meta.at[i, f'SUM_{score_type}r{str(r+1)}_pos']
            res_sum_LFC3D_pos += float(LFC3D_pos)
        
        list_avg_LFC3D_neg.append(res_sum_LFC3D_neg / nRandom)
        list_avg_LFC3D_pos.append(res_sum_LFC3D_pos / nRandom)

    df_meta[f'AVG_{score_type}r_neg'] = list_avg_LFC3D_neg 
    df_meta[f'AVG_{score_type}r_pos'] = list_avg_LFC3D_pos
    del list_avg_LFC3D_neg, list_avg_LFC3D_pos

    # DELETE ALL RANDOM COLUMNS #
    df_meta = df_meta.drop(columns=[f'SUM_{score_type}r{str(r+1)}_neg' for r in range(0, nRandom)])
    df_meta = df_meta.drop(columns=[f'SUM_{score_type}r{str(r+1)}_pos' for r in range(0, nRandom)])

    # CONVERT SIGNAL TO Z SCORE #
    colnames = [f'SUM_{score_type}_{sign}' for sign in ['neg', 'pos']]
    params = [{'mu': df_meta[f'AVG_{score_type}r_{sign}'].mean(),
                's': df_meta[f'AVG_{score_type}r_{sign}'].std()} 
                for sign in ['neg', 'pos']]
    result_data = {f'SUM_{score_type}_{sign}_{suffix}': [] 
                    for sign in ['neg', 'pos'] for suffix in ['z', 'p', 'psig']}

    for i in range(len(df_meta)):
        for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
            signal = float(df_meta.at[i, colname])
            signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
            
            # Append results to the dictionary
            result_data[f'SUM_{score_type}_{sign}_z'].append(signal_z)
            result_data[f'SUM_{score_type}_{sign}_p'].append(signal_p)
            result_data[f'SUM_{score_type}_{sign}_psig'].append(signal_plabel)

    df_meta_Z = pd.concat([df_meta, pd.DataFrame(result_data)], axis=1).round(4)

    filename = edits_filedir / f"metaaggregation/{structureid}_MetaAggr_{score_type}.tsv"
    df_meta_Z.to_csv(filename, "\t", index=False)
    
    # PLOTS #
    LFC3D_plots(df_meta_Z, edits_filedir, input_gene, pthr, func=aggr_func.__name__.upper(), score=score_type)

    return df_meta_Z
