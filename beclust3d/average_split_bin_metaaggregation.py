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
warnings.filterwarnings('ignore')

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
        list_LFC3D_neg.append(aggr_func(values_LFC3D_neg) if values_LFC3D_neg else None)
        list_LFC3D_pos.append(aggr_func(values_LFC3D_pos) if values_LFC3D_pos else None)

    df_bidir_meta[f'{aggr_func_name}_{score_type}_neg'] = list_LFC3D_neg
    df_bidir_meta[f'{aggr_func_name}_{score_type}_pos'] = list_LFC3D_pos
    list_LFC3D_neg = [np.nan if x == '-' else x for x in list_LFC3D_neg]
    list_LFC3D_pos = [np.nan if x == '-' else x for x in list_LFC3D_pos]
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
        new_col_neg = df_LFC_LFC3D[headers_neg].sum(axis=1).rename(f"SUM_{score_type}r{str(n)}_neg")
        new_col_pos = df_LFC_LFC3D[headers_pos].sum(axis=1).rename(f"SUM_{score_type}r{str(n)}_pos")
        df_bidir_meta = pd.concat([df_bidir_meta, new_col_neg, new_col_pos], axis=1)


    # AVG ACROSS ALL RANDOMIZATIONS #
    headers_neg = [f"SUM_{score_type}r{str(n)}_neg" for n in range(1, nRandom+1)]
    headers_pos = [f"SUM_{score_type}r{str(n)}_pos" for n in range(1, nRandom+1)]
    df_bidir_meta[f"SUM_{score_type}r_neg"] = df_bidir_meta[headers_neg].mean(axis=1)
    df_bidir_meta[f"SUM_{score_type}r_pos"] = df_bidir_meta[headers_pos].mean(axis=1)
    df_bidir_meta = df_bidir_meta

    # BINNING #
    df_LFC_LFC3D_dis = df_bidir_meta[['unipos', 'unires', header_main, 
                                      f'SUM_{score_type}r_neg', f'SUM_{score_type}r_pos']].copy()
    quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}

    # GENERATE THRESHOLDS FOR BINNING #
    df_LFC3D_neg = df_bidir_meta.loc[df_bidir_meta[header_main] < 0.0, ].reset_index(drop=True)
    df_LFC3D_pos = df_bidir_meta.loc[df_bidir_meta[header_main] > 0.0, ].reset_index(drop=True)
    df_neg_stats = df_LFC3D_neg[header_main].describe()
    df_pos_stats = df_LFC3D_pos[header_main].describe()
    # print(df_neg_stats)
    # print(df_pos_stats)

    # CALCULATE BINS #
    quantile_values = {}
    for name, q in quantiles.items(): 
        quantile_values[name] = df_bidir_meta[header_main].quantile(q)

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

    df_meta_Z = pd.concat([df_bidir_meta, pd.DataFrame(result_data)], axis=1)

    filename = edits_filedir / f"metaaggregation/{structureid}_MetaAggr_{score_type}.tsv"
    df_meta_Z.to_csv(filename, "\t", index=False)

    return df_bidir_meta, df_LFC_LFC3D_dis, df_meta_Z


def average_split_meta(
    df_LFC_LFC3D, 
    workdir, input_gene, screen_names, 
    score_type='LFC3D', aggr_func=np.sum, aggr_func_name='SUM', nRandom=1000, 
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

    # AGGR LFC3D VALUES ACROSS SCREENS FOR EACH RESIDUE #
    list_LFC3D_neg, list_LFC3D_pos = [], []
    header_scores = [f"{sn}_{score_type}" for sn in screen_names]
    screen_name_dicts = [df_LFC_LFC3D[header].to_dict() for header in header_scores]

    # AGGR A VALUE OR '_' FOR EACH RESIDUE #
    for i in range(len(df_LFC_LFC3D)): 
        values_LFC3D_neg, values_LFC3D_pos = [], []

        # ADD POS AND NEG VALS SEPARATELY FOR EACH RESIDUE #
        for screen_dict in screen_name_dicts: 
            LFC3D = screen_dict[i]
            if LFC3D != '-': 
                LFC3D_value = float(LFC3D)
                if LFC3D_value < 0.0: 
                    values_LFC3D_neg.append(LFC3D_value)
                elif LFC3D_value > 0.0: 
                    values_LFC3D_pos.append(LFC3D_value)

        # APPLY AGGR FUNCTION #
        list_LFC3D_neg.append(aggr_func(values_LFC3D_neg) if values_LFC3D_neg else '-')
        list_LFC3D_pos.append(aggr_func(values_LFC3D_pos) if values_LFC3D_pos else '-')

    df_bidir_meta[f'{aggr_func_name}_{score_type}_neg'] = list_LFC3D_neg
    df_bidir_meta[f'{aggr_func_name}_{score_type}_pos'] = list_LFC3D_pos
    df_bidir_meta[f'{aggr_func_name}_{score_type}'] = [sum_dash(x) for x in zip(list_LFC3D_neg, list_LFC3D_pos)]
    del list_LFC3D_neg, list_LFC3D_pos

    # PULL RANDOMIZED DATA #
    # because we randomize the data then split into pos and neg, 
    # there may be neg randomized data for rows where there is no neg value, 
    # and same for positive 
    for n in range(1, nRandom+1): 
        new_col_neg_list, new_col_pos_list = [], []
        for sn in screen_names: # SPLIT INTO POS NEG #
            colnam = f"{sn}_{score_type}r{str(n)}"
            new_col_neg_list.append(df_LFC_LFC3D[f"{colnam}"].apply(lambda x: process_col(x, 'neg')).rename(f"{colnam}_neg"))
            new_col_pos_list.append(df_LFC_LFC3D[f"{colnam}"].apply(lambda x: process_col(x, 'pos')).rename(f"{colnam}_pos"))
        df_temp = pd.concat(new_col_neg_list + new_col_pos_list, axis=1)
        
        # SUM ACROSS ALL SCREENS #
        headers_neg = [f"{sn}_{score_type}r{str(n)}_neg" for sn in screen_names]
        headers_pos = [f"{sn}_{score_type}r{str(n)}_pos" for sn in screen_names]
        aggr_col_neg = df_temp[headers_neg].replace('-', np.nan).sum(axis=1)
        aggr_col_pos = df_temp[headers_pos].replace('-', np.nan).sum(axis=1)
        aggr_col_neg = aggr_col_neg.rename(f"SUM_{score_type}r{str(n)}_neg").replace(0.0, '-')
        aggr_col_pos = aggr_col_pos.rename(f"SUM_{score_type}r{str(n)}_pos").replace(0.0, '-')
        df_bidir_meta = pd.concat([df_bidir_meta, aggr_col_neg, aggr_col_pos], axis=1)

    # AVG ACROSS ALL RANDOMIZATIONS #
    headers_neg = [f"SUM_{score_type}r{str(n)}_neg" for n in range(1, nRandom+1)]
    headers_pos = [f"SUM_{score_type}r{str(n)}_pos" for n in range(1, nRandom+1)]
    new_col_neg = df_bidir_meta[headers_neg].replace({'-':np.nan, 0.0:np.nan}).mean(axis=1)
    new_col_pos = df_bidir_meta[headers_pos].replace({'-':np.nan, 0.0:np.nan}).mean(axis=1)
    new_col_neg = new_col_neg.rename(f"SUM_{score_type}r_neg").replace(0.0, '-')
    new_col_pos = new_col_pos.rename(f"SUM_{score_type}r_pos").replace(0.0, '-')
    df_bidir_meta = pd.concat([df_bidir_meta, new_col_neg, new_col_pos], axis=1)

    # headers_neg = [f"SUM_{score_type}r{str(n)}_neg" for n in range(1, nRandom+1)]
    # headers_pos = [f"SUM_{score_type}r{str(n)}_pos" for n in range(1, nRandom+1)]
    # df_bidir_meta[headers_neg] = df_bidir_meta[headers_neg].replace({'-':np.nan, 0.0:np.nan}).apply(pd.to_numeric, errors='coerce')
    # df_bidir_meta[headers_pos] = df_bidir_meta[headers_pos].replace({'-':np.nan, 0.0:np.nan}).apply(pd.to_numeric, errors='coerce')
    # sum_pos = df_bidir_meta[headers_pos].where(df_bidir_meta[headers_pos] > 0).sum(axis=1, skipna=True)
    # count_pos = df_bidir_meta[headers_pos].where(df_bidir_meta[headers_pos] > 0).count(axis=1)
    # sum_neg = df_bidir_meta[headers_neg].where(df_bidir_meta[headers_neg] < 0).sum(axis=1, skipna=True)
    # count_neg = df_bidir_meta[headers_neg].where(df_bidir_meta[headers_neg] < 0).count(axis=1)
    # new_col_pos = (sum_pos / count_pos).fillna(0).rename(f"SUM_{score_type}r_pos").replace({0.0: '-'})
    # new_col_neg = (sum_neg / count_neg).fillna(0).rename(f"SUM_{score_type}r_neg").replace({0.0: '-'})
    # df_bidir_meta = pd.concat([df_bidir_meta, new_col_pos, new_col_neg], axis=1)

    # new_col_neg = df_bidir_meta[headers_neg].replace('-', np.nan).std(axis=1) ###
    # new_col_pos = df_bidir_meta[headers_pos].replace('-', np.nan).std(axis=1) ###
    # new_col_neg = new_col_neg.rename(f"SUM_{score_type}r_neg_stdev").replace(0.0, '-') ###
    # new_col_pos = new_col_pos.rename(f"SUM_{score_type}r_pos_stdev").replace(0.0, '-') ###
    # df_bidir_meta = pd.concat([df_bidir_meta, new_col_neg, new_col_pos], axis=1)
    
    out_filename_bidir = edits_filedir / f"metaaggregation/{input_gene}_{score_type}_bidirectional.tsv"
    df_bidir_meta.to_csv(out_filename_bidir, sep='\t', index=False)

    return df_bidir_meta

def bin_meta(
    df_bidir_meta, 
    workdir, input_gene, 
    score_type='LFC3D', aggr_func_name='SUM', 
): 
    """
    Description
        A point to meta aggregate across multiple screens or just one screen, 
        calculate signal vs background, bin this new metaaggregated signal 
        into top and bottom 10 %, and plot QC graphics

    Params
        df_bidir_meta: pandas dataframe
            from previous step average_split_meta()
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

    header_main = f'{aggr_func_name}_{score_type}'
    random_neg, random_pos = f'SUM_{score_type}r_neg', f'SUM_{score_type}r_pos'
    headers = [header_main, f'{header_main}_neg', f'{header_main}_pos', random_neg, random_pos]
    # BINNING #
    df_LFC_LFC3D_dis = df_bidir_meta[['unipos', 'unires'] + headers].copy()
    quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}

    # GENERATE THRESHOLDS FOR BINNING #
    mask_neg = df_LFC_LFC3D_dis[random_neg] != 0.0
    mask_pos = df_LFC_LFC3D_dis[random_pos] != 0.0
    df_neg_stats = df_LFC_LFC3D_dis[random_neg][mask_neg].describe()
    df_pos_stats = df_LFC_LFC3D_dis[random_pos][mask_pos].describe()

    # CALCULATE BINS #
    quantile_values = {}
    for name, q in quantiles.items(): 
        quantile_values[name] = df_LFC_LFC3D_dis[header_main].replace('-', np.nan).astype(float).quantile(q)

    arr_disc, arr_weight = binning_neg_pos(df_bidir_meta, df_neg_stats, df_pos_stats, 
                                           quantile_values.values(), header_main)
    df_LFC_LFC3D_dis[f"{header_main}_dis"]  = arr_disc
    df_LFC_LFC3D_dis[f"{header_main}_wght"] = arr_weight

    out_filename_dis = edits_filedir / f"metaaggregation/{input_gene}_{score_type}_dis_wght.tsv"
    df_LFC_LFC3D_dis.to_csv(out_filename_dis, sep = '\t', index=False)

    return df_LFC_LFC3D_dis, df_neg_stats, df_pos_stats

def znorm_meta(
    df_bidir_meta, df_neg_stats, df_pos_stats, 
    workdir, input_gene, 
    pthrs=[0.05, 0.01, 0.001], score_type='LFC3D', aggr_func_name='SUM', 
): 
    """
    Description
        A point to meta aggregate across multiple screens or just one screen, 
        calculate signal vs background, bin this new metaaggregated signal 
        into top and bottom 10 %, and plot QC graphics

    Params
        df_bidir_meta: pandas dataframe
            from previous step average_split_meta()
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

    header_main = f'{aggr_func_name}_{score_type}'
    pthrs_str = [str(pthr).replace('.','') for pthr in pthrs]
    # CONVERT SIGNAL TO Z SCORE #
    colnames = [f'{header_main}_{sign}' for sign in ['neg', 'pos']]
    params = [{'mu': df_neg_stats['mean'], 's': df_neg_stats['std']}, 
              {'mu': df_pos_stats['mean'], 's': df_pos_stats['std']}, ]
    # params = [[{'mu': df_bidir_meta[f"SUM_{score_type}r_neg"].iloc[i], 's': df_bidir_meta[f"SUM_{score_type}r_neg_stdev"].iloc[i]} for i in range(len(df_bidir_meta))], 
    #           [{'mu': df_bidir_meta[f"SUM_{score_type}r_pos"].iloc[i], 's': df_bidir_meta[f"SUM_{score_type}r_pos_stdev"].iloc[i]} for i in range(len(df_bidir_meta))], ]
    result_data = {f'{header_main}_{sign}_{pthr_str}_{suffix}': [] 
                   for sign in ['neg', 'pos'] for suffix in ['z', 'p', 'psig'] for pthr_str in pthrs_str}

    for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
        signals_dict = df_bidir_meta[colname].replace('-', np.nan).to_dict()

        for pthr, pthr_str in zip(pthrs, pthrs_str): 
            for i in range(len(df_bidir_meta)):
                signal = float(signals_dict[i])
                signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
                # signal_z, signal_p, signal_plabel = calculate_stats(signal, param[i], pthr)
                
                # Append results to the dictionary
                result_data[f'{header_main}_{sign}_{pthr_str}_z'].append(signal_z)
                result_data[f'{header_main}_{sign}_{pthr_str}_p'].append(signal_p)
                result_data[f'{header_main}_{sign}_{pthr_str}_psig'].append(signal_plabel)

    df_meta_Z = pd.concat([df_bidir_meta, pd.DataFrame(result_data).replace(0,'-')], axis=1)
    df_meta_Z = df_meta_Z

    filename = edits_filedir / f"metaaggregation/{input_gene}_MetaAggr_{score_type}.tsv"
    df_meta_Z.to_csv(filename, "\t", index=False)

    return df_meta_Z

def sum_dash(values): 
    new_values = [x for x in values if x != '-']
    if len(new_values) == 0: return '-'
    else: return sum(new_values)

def process_col(x, mode): 
    if mode == 'neg': 
        return float(x) if x != '-' and float(x) < 0 else np.nan
    if mode == 'pos': 
        return float(x) if x != '-' and float(x) > 0 else np.nan
