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
    input_gene, structureid, input_screens, 
    nRandom=1000, pthr=0.05, score_type='LFC3D', 
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

    Returns
        df_meta: pandas DataFrame
            results of meta aggregated screen data
    """
    
    edits_filedir = Path(workdir + '/' + input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'metaaggregation'):
        os.mkdir(edits_filedir / 'metaaggregation')
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    # Sum LFC3D #
    df_meta = pd.DataFrame()
    df_meta['unipos'] = df_LFC_LFC3D['unipos']
    df_meta['unires'] = df_LFC_LFC3D['unires']
    
    # SUM LFC3D VALUES ACROSS SCREENS #
    list_sum_LFC3D_neg, list_sum_LFC3D_pos = [], []
    for i in range(0, len(df_LFC_LFC3D)):
        res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0

        for input_screen in input_screens: 
            screen_name = input_screen.split('.')[0]
            header_LFC3D = f"{screen_name}_{score_type}"
            if header_LFC3D not in df_LFC_LFC3D.columns: 
                warnings.warn(f'{header_LFC3D} not in input df_LFC_LFC3D')
                continue
            
            LFC3D = df_LFC_LFC3D.at[i, header_LFC3D]
            if (LFC3D != '-') and (float(LFC3D) < 0.0):
                res_sum_LFC3D_neg += float(LFC3D)
            if (LFC3D != '-') and (float(LFC3D) > 0.0):
                res_sum_LFC3D_pos += float(LFC3D)
        
        list_sum_LFC3D_neg.append(res_sum_LFC3D_neg)
        list_sum_LFC3D_pos.append(res_sum_LFC3D_pos)

    df_meta[f'SUM_{score_type}_neg'] = list_sum_LFC3D_neg 
    df_meta[f'SUM_{score_type}_pos'] = list_sum_LFC3D_pos
    del list_sum_LFC3D_neg, list_sum_LFC3D_pos

    # RANDOMIZE #
    dict_META = {}
    for r in range(0, nRandom):
        list_sum_LFC3D_neg, list_sum_LFC3D_pos = [], []

        for i in range(0, len(df_LFC_LFC3D)): 
            res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0

            for input_screen in input_screens: # sum across screens
                screen_name = input_screen.split('.')[0]
                header_LFC3D = f"{screen_name}_{score_type}r{str(r+1)}"
                if header_LFC3D not in df_LFC_LFC3D.columns: 
                    # warnings.warn(f'{header_LFC3D} not in input df_LFC_LFC3D')
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
    for i in range(0, len(df_meta)): 
        res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0

        for r in range(0, nRandom): 
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
    colnames = [f'SUM_{score_type}_neg', f'SUM_{score_type}_pos']
    params = [{'mu':df_meta[f'AVG_{score_type}r_neg'].mean(), 's':df_meta[f'AVG_{score_type}r_neg'].std()}, 
              {'mu':df_meta[f'AVG_{score_type}r_pos'].mean(), 's':df_meta[f'AVG_{score_type}r_pos'].std()} ]
    lists_LFC3D = [{'z':[], 'p':[], 'lab':[]}, {'z':[], 'p':[], 'lab':[]}, ]

    for i in range(0, len(df_meta)):
        for colname, param, lists in zip(colnames, params, lists_LFC3D): 
            signal = float(df_meta.at[i, colname])
            if param['s'] == 0: 
                lists['z'].append(0)
                lists['p'].append(0)
                lists['lab'].append(0)
            else: 
                signal_z = statistics.NormalDist(mu=param['mu'], sigma=param['s']).zscore(signal)
                signal_p = stats.norm.sf(abs(signal_z))
                signal_plabel = f'p<{str(pthr)}' if signal_p < pthr else f'p>={str(pthr)}'                
                lists['z'].append(signal_z)
                lists['p'].append(signal_p)
                lists['lab'].append(signal_plabel)
    
    dict_META = {f'SUM_{score_type}_neg_z':lists_LFC3D[0]['z'], f'SUM_{score_type}_neg_p':lists_LFC3D[0]['p'], 
                 f'SUM_{score_type}_neg_psig':lists_LFC3D[0]['lab'], f'SUM_{score_type}_pos_z':lists_LFC3D[1]['z'], 
                 f'SUM_{score_type}_pos_p':lists_LFC3D[1]['p'], f'SUM_{score_type}_pos_psig':lists_LFC3D[1]['lab'], }
    df_meta_Z = pd.concat([df_meta, pd.DataFrame(dict_META)], axis=1)
    df_meta_Z.round(4)

    filename = edits_filedir / f"metaaggregation/{structureid}_MetaAggr_{score_type}_and_randomized_background.tsv"
    df_meta_Z.to_csv(filename, "\t", index=False)
    
    # PLOTS #
    LFC3D_plots(df_meta_Z, edits_filedir, input_gene, pthr, score_type=score_type, )

    return df_meta_Z
    