"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 4

"""

import pandas as pd
from pathlib import Path
import os
from _average_split_bin_helpers_ import *
import warnings
# warnings.filterwarnings('ignore')

def metaaggregation(
    df_LFC_LFC3D, workdir, 
    input_gene, structureid, screen_names, 
    pthr=0.05, score_type='LFC3D', aggr_func=np.sum, aggr_func_name='SUM', nRandom=1000, 
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
        df_bidir_meta: pandas DataFrame
            results of meta aggregated screen data
    """
    
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'metaaggregation'):
        os.mkdir(edits_filedir / 'metaaggregation')

    # AGGREGATE LFC3D #
    df_bidir_meta = pd.DataFrame()
    df_bidir_meta['unipos'] = df_LFC_LFC3D['unipos']
    df_bidir_meta['unires'] = df_LFC_LFC3D['unires']
    header_main = f'{aggr_func_name}_{score_type}'

    # AGGR LFC3D VALUES ACROSS SCREENS FOR EACH RESIDUE #
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
                if LFC3D_value < 0.0: 
                    values_LFC3D_neg.append(LFC3D_value)
                elif LFC3D_value > 0.0: 
                    values_LFC3D_pos.append(LFC3D_value)

        # APPLY AGGR FUNCTION #
        list_LFC3D_neg.append(aggr_func(values_LFC3D_neg) if values_LFC3D_neg else 0.0)
        list_LFC3D_pos.append(aggr_func(values_LFC3D_pos) if values_LFC3D_pos else 0.0)

    df_bidir_meta[f'{aggr_func_name}_{score_type}_neg'] = list_LFC3D_neg
    df_bidir_meta[f'{aggr_func_name}_{score_type}_pos'] = list_LFC3D_pos
    df_bidir_meta[header_main] = [sum(x) for x in zip(list_LFC3D_neg, list_LFC3D_pos)]
    del list_LFC3D_neg, list_LFC3D_pos

    # PULL RANDOMIZED DATA #
    for n in range(1, nRandom+1): 
        for sn in screen_names: # SPLIT INTO POS NEG #
            pref = f"{sn}_{score_type}r{str(n)}"
            df_LFC_LFC3D[f"{pref}_neg"] = df_LFC_LFC3D[f"{pref}"].apply(lambda x: float(x) if x != '-' and float(x) < 0 else np.nan)
            df_LFC_LFC3D[f"{pref}_pos"] = df_LFC_LFC3D[f"{pref}"].apply(lambda x: float(x) if x != '-' and float(x) > 0 else np.nan)
        
        # SUM ACROSS ALL SCREENS #
        headers_neg = [f"{sn}_{score_type}r{str(n)}_neg" for sn in screen_names]
        headers_pos = [f"{sn}_{score_type}r{str(n)}_pos" for sn in screen_names]
        df_bidir_meta[f"SUM_{score_type}r{str(n)}_neg"] = df_LFC_LFC3D[headers_neg].sum(axis=1)
        df_bidir_meta[f"SUM_{score_type}r{str(n)}_pos"] = df_LFC_LFC3D[headers_pos].sum(axis=1)

    # AVG ACROSS ALL RANDOMIZATIONS #
    headers_neg = [f"SUM_{score_type}r{str(n)}_neg" for n in range(1, nRandom+1)]
    headers_pos = [f"SUM_{score_type}r{str(n)}_pos" for n in range(1, nRandom+1)]
    df_bidir_meta[f"SUM_{score_type}r_neg"] = df_bidir_meta[headers_neg].mean(axis=1)
    df_bidir_meta[f"SUM_{score_type}r_pos"] = df_bidir_meta[headers_pos].mean(axis=1)

    # lfc3dr_colnames_neg = [f"{screen_name}_AVG_{score_type}r_neg" for screen_name in screen_names]
    # df_bidir_meta[f'SUM_{score_type}r_neg'] = df_LFC_LFC3D[lfc3dr_colnames_neg].sum(axis=1)
    # lfc3dr_colnames_pos = [f"{screen_name}_AVG_{score_type}r_pos" for screen_name in screen_names]
    # df_bidir_meta[f'SUM_{score_type}r_pos'] = df_LFC_LFC3D[lfc3dr_colnames_pos].sum(axis=1)

    # BINNING #
    df_LFC_LFC3D_dis = df_bidir_meta[['unipos', 'unires', header_main, 
                                      f'SUM_{score_type}r_neg', f'SUM_{score_type}r_pos']].copy()
    quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}

    # GENERATE THRESHOLDS FOR BINNING #
    # df_nodash = df_bidir_meta.loc[df_bidir_meta[header_main] != 0.0, ].reset_index(drop=True)
    # df_nodash[header_main] = df_nodash[header_main].astype(float)
    # df_LFC3D_neg = df_nodash.loc[df_nodash[header_main] < 0, ].reset_index(drop=True)
    # df_LFC3D_pos = df_nodash.loc[df_nodash[header_main] > 0, ].reset_index(drop=True)
    df_neg_stats = df_bidir_meta[f'SUM_{score_type}r_neg'].describe()
    df_pos_stats = df_bidir_meta[f'SUM_{score_type}r_pos'].describe()
    print(df_neg_stats)
    print(df_pos_stats)

    # CALCULATE BINS #
    quantile_values = {}
    for name, q in quantiles.items(): 
        quantile_values[name] = round(df_bidir_meta[header_main].quantile(q), 4)

    arr_disc, arr_weight = binning_neg_pos(df_bidir_meta, df_neg_stats, df_pos_stats, 
                                           quantile_values.values(), header_main)
    df_LFC_LFC3D_dis[f"{aggr_func_name}_{score_type}_dis"]  = arr_disc
    df_LFC_LFC3D_dis[f"{aggr_func_name}_{score_type}_wght"] = arr_weight

    out_filename_bidir = edits_filedir / f"metaaggregation/{input_gene}_{score_type}_bidirectional.tsv"
    df_bidir_meta.to_csv(out_filename_bidir, sep='\t', index=False)
    out_filename_dis = edits_filedir / f"metaaggregation/{input_gene}_{score_type}_dis_wght.tsv"
    df_LFC_LFC3D_dis.to_csv(out_filename_dis, sep = '\t', index=False)

    # CONVERT SIGNAL TO Z SCORE #
    colnames = [f'{aggr_func_name}_{score_type}_{sign}' for sign in ['neg', 'pos']]
    params = [{'mu': df_neg_stats['mean'], 's': df_neg_stats['std']}, 
              {'mu': df_pos_stats['mean'], 's': df_pos_stats['std']}, ]
    result_data = {f'{aggr_func_name}_{score_type}_{sign}_{suffix}': [] 
                   for sign in ['neg', 'pos'] for suffix in ['z', 'p', 'psig']}

    for i in range(len(df_bidir_meta)):
        for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
            signal = float(df_bidir_meta.at[i, colname])
            signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
            
            # Append results to the dictionary
            result_data[f'{aggr_func_name}_{score_type}_{sign}_z'].append(signal_z)
            result_data[f'{aggr_func_name}_{score_type}_{sign}_p'].append(signal_p)
            result_data[f'{aggr_func_name}_{score_type}_{sign}_psig'].append(signal_plabel)

    df_meta_Z = pd.concat([df_bidir_meta, pd.DataFrame(result_data)], axis=1).round(4)

    filename = edits_filedir / f"metaaggregation/{structureid}_MetaAggr_{score_type}.tsv"
    df_meta_Z.to_csv(filename, "\t", index=False)

    return df_bidir_meta, df_LFC_LFC3D_dis, df_meta_Z
