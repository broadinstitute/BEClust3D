"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 4

"""

import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
import statistics
import scipy.stats as stats
from pathlib import Path
from scipy.stats import mannwhitneyu
import os
from binning_lfc3d import binning
import warnings

def metaaggregation(
    df_LFC_LFC3D,
    workdir, 
    input_gene, structureid, 
    input_screens, 
    nRandom=1000, 
    pthr=0.05, 
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
        df_METAggregation: pandas DataFrame
            results of meta aggregated screen data
    """
    
    # Sum LFC3D #

    df_METAggregation = pd.DataFrame()
    df_METAggregation['unipos'] = df_LFC_LFC3D['unipos']
    df_METAggregation['unires'] = df_LFC_LFC3D['unires']
    
    list_sum_LFC3D_neg, list_sum_LFC3D_pos = [], []
    for i in range(0, len(df_LFC_LFC3D)):
        res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0

        for input_screen in input_screens: 
            screen_name = input_screen.split('.')[0]
            header_LFC3D = f"{screen_name}_LFC3D"
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

    df_METAggregation['SUM_LFC3D_neg'] = list_sum_LFC3D_neg 
    df_METAggregation['SUM_LFC3D_pos'] = list_sum_LFC3D_pos
    del list_sum_LFC3D_neg, list_sum_LFC3D_pos

    dict_META = {}
    for r in range(0, nRandom):
        list_sum_LFC3D_neg, list_sum_LFC3D_pos = [], []

        for i in range(0, len(df_LFC_LFC3D)): 
            res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0

            for input_screen in input_screens: # sum across screens
                screen_name = input_screen.split('.')[0]
                header_LFC3D = f"{screen_name}_LFC3Dr{str(r+1)}"
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

        dict_META[f'SUM_LFC3Dr{str(r+1)}_neg'] = list_sum_LFC3D_neg 
        dict_META[f'SUM_LFC3Dr{str(r+1)}_pos'] = list_sum_LFC3D_pos
        del list_sum_LFC3D_neg, list_sum_LFC3D_pos

    df_METAggregation = pd.concat([df_METAggregation, pd.DataFrame(dict_META)], axis=1)
    del dict_META

    # Average of SUM of randomized signal, as the background noise #

    list_avg_LFC3D_neg, list_avg_LFC3D_pos = [], []
    for i in range(0, len(df_METAggregation)): 
        res_sum_LFC3D_neg, res_sum_LFC3D_pos = 0.0, 0.0

        for r in range(0, nRandom): 
            LFC3D_neg = df_METAggregation.at[i, f'SUM_LFC3Dr{str(r+1)}_neg']
            res_sum_LFC3D_neg += float(LFC3D_neg)
            LFC3D_pos = df_METAggregation.at[i, f'SUM_LFC3Dr{str(r+1)}_pos']
            res_sum_LFC3D_pos += float(LFC3D_pos)
        
        list_avg_LFC3D_neg.append(res_sum_LFC3D_neg / nRandom)
        list_avg_LFC3D_pos.append(res_sum_LFC3D_pos / nRandom)

    df_METAggregation['AVG_LFC3Dr_neg'] = list_avg_LFC3D_neg 
    df_METAggregation['AVG_LFC3Dr_pos'] = list_avg_LFC3D_pos
    del list_avg_LFC3D_neg, list_avg_LFC3D_pos

    # convert signal to Z score #

    lists_LFC3D = [{'z':[], 'p':[], 'lab':[]}, 
                    {'z':[], 'p':[], 'lab':[]}, 
                    ]
    params = [{'mu':df_METAggregation['AVG_LFC3Dr_neg'].mean(), 
                's':df_METAggregation['AVG_LFC3Dr_neg'].std()}, 
                {'mu':df_METAggregation['AVG_LFC3Dr_pos'].mean(), 
                's':df_METAggregation['AVG_LFC3Dr_pos'].std()}]
    colnames = ['SUM_LFC3D_neg', 'SUM_LFC3D_pos']

    for i in range(0, len(df_METAggregation)):
        # sensitizing and resistant
        for colname, param, lists in zip(colnames, params, lists_LFC3D): 
            signal = float(df_METAggregation.at[i, colname])
            if param['s'] == 0: 
                lists['z'].append(0)
                lists['p'].append(0)
                lists['lab'].append(0)
                continue
            
            signal_z = statistics.NormalDist(mu=param['mu'], sigma=param['s']).zscore(signal)
            signal_p = stats.norm.sf(abs(signal_z))
            signal_plabel = f'p<{str(pthr)}' if signal_p < pthr else f'p>={str(pthr)}'                
            lists['z'].append(signal_z)
            lists['p'].append(signal_p)
            lists['lab'].append(signal_plabel)
    
    dict_META = {'SUM_LFC3D_neg_z':lists_LFC3D[0]['z'], 'SUM_LFC3D_neg_p':lists_LFC3D[0]['p'], 'SUM_LFC3D_neg_psig':lists_LFC3D[0]['lab'], 
                 'SUM_LFC3D_pos_z':lists_LFC3D[1]['z'], 'SUM_LFC3D_pos_p':lists_LFC3D[1]['p'], 'SUM_LFC3D_pos_psig':lists_LFC3D[1]['lab'], }
    df_METAggregation = pd.concat([df_METAggregation, pd.DataFrame(dict_META)], axis=1)

    edits_filedir = Path(workdir + '/' + input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'metaaggregation'):
        os.mkdir(edits_filedir / 'metaaggregation')
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    filename = edits_filedir / f"metaaggregation/{structureid}_MetaAggr_LFC3D_and_randomized_background.tsv"
    df_METAggregation.to_csv(filename, "\t", index=False)

    # HISTOGRAMS #
    res_neg = metaaggregation_histogram(
        df_METAggregation, edits_filedir, input_gene, 
        avg_colname='AVG_LFC3Dr_neg', sum_colname='SUM_LFC3D_neg', out_keyword='sensitizing', 
    )
    res_pos = metaaggregation_histogram(
        df_METAggregation, edits_filedir, input_gene, 
        avg_colname='AVG_LFC3Dr_pos', sum_colname='SUM_LFC3D_pos', out_keyword='resistant', 
        )
    if res_neg is None or res_pos is None: 
        return None
    print(res_neg)
    print(res_pos)
        
    df_METAggregation = binning_lfc3d(df_METAggregation)
    # DISPLOTS #
    displots = [('SUM_LFC3D_neg_dis', 'SUM_LFC3D_neg', 'SUM_LFC3D_neg_psig', 'neg_pvalue'), 
                ('SUM_LFC3D_neg_dis', 'SUM_LFC3D_neg', 'SUM_LFC3D_neg_dis', 'neg_pvalue'), 
                ('SUM_LFC3D_pos_dis', 'SUM_LFC3D_pos', 'SUM_LFC3D_pos_psig', 'pos_pvalue'), 
                ('SUM_LFC3D_pos_dis', 'SUM_LFC3D_pos', 'SUM_LFC3D_pos_dis', 'pos_pvalue'), 
                ]
    for filter, x, hue, name in displots: 
        metaaggregation_displot(
            df_METAggregation, 
            filter_col=filter, x_col=x, hue_col=hue, 
            filedir=edits_filedir, input_gene=input_gene, out_keyword=name
        )
    # SCATTERPLOT #
    scatterplots = [('SUM_LFC3D_neg_dis', 'SUM_LFC3D_neg_psig', 'SUM_LFC3D_neg', 'neg_dis'), 
                    ('SUM_LFC3D_pos_dis', 'SUM_LFC3D_pos_psig', 'SUM_LFC3D_pos', 'pos_dis')]
    for dis, pval, y, out in scatterplots: 
        metaaggregation_scatterplot(
                df_METAggregation, 
                dis_col=dis, pval_col=pval, y_col=y,
                filedir=edits_filedir, input_gene=input_gene, out_keyword=out, pthr=pthr, 
        )
    
    return df_METAggregation


def metaaggregation_histogram(
        df_METAggregation, filedir, input_gene, 
        avg_colname, sum_colname, out_keyword, 
): 
    """
    Description
        Helper function to plot histograms of the values along the length of the gene
    """
    
    plt.figure(figsize=(16, 8), dpi=300)

    results = {}
    df_METAggregation_plot = pd.DataFrame()
    df_METAggregation_plot['unipos'] = df_METAggregation['unipos']
    df_METAggregation_plot[sum_colname] = df_METAggregation[sum_colname]
    df_METAggregation_plot[avg_colname] = df_METAggregation[avg_colname]

    U1, p = mannwhitneyu(df_METAggregation_plot[sum_colname], 
                         df_METAggregation_plot[avg_colname], method="asymptotic")
    results['mannwhitneyu U1'], results['mannwhitneyu p'] = U1, p
    r, p = stats.pearsonr(df_METAggregation_plot[sum_colname], 
                          df_METAggregation_plot[avg_colname])
    results['pearsonr r'], results['pearsonr p'] = r, p

    # SUM #
    temp = df_METAggregation_plot[sum_colname]
    results['sum min'], results['sum mean'] = temp.min(), temp.mean()
    results['sum med'], results['sum std'] = temp.median(), temp.std()

    if results['sum std'] == 0: 
        return None
    z = statistics.NormalDist(mu=results['sum mean'], sigma=results['sum std']).zscore(-4.6)
    results['z'], results['p cdf'], results['p sf'] = z, stats.norm.cdf(z), stats.norm.sf(abs(z))

    # AVG #
    temp = df_METAggregation_plot[avg_colname]
    results['avg min'], results['avg mean'] = temp.min(), temp.mean()
    results['avg med'], results['avg std'] = temp.median(), temp.std()

    # PLOT #
    df_METAggregation_plot.plot.area(x='unipos', alpha=0.55, stacked = False)
    plt.axhline(y = results['sum mean'], color = 'r', linestyle = '-')
    plt.axhline(y = results['sum mean']-results['sum std'], color = 'r', linestyle = '--')

    plt.legend(loc='lower left', borderaxespad=0)
    plt.xticks(np.arange(0,len(df_METAggregation), 100))
    plt.title(out_keyword)

    out_filename = filedir / f"plots/{input_gene}_signal_vs_background_{out_keyword}.png"
    plt.savefig(out_filename, dpi = 500)
    del df_METAggregation_plot

    return results

def binning_lfc3d(
        df_METAggregation, 
        colnames = ['SUM_LFC3D_neg', 'SUM_LFC3D_pos'], 
): 
    """
    Description
        Helper function to bin the top 10 and bottom 10 % of points
    """
    
    df_3daggr_list = [pd.DataFrame() for _ in colnames]
    quantile_numbers = {colnames[0]: (0.1, 0.05), 
                        colnames[1]: (0.9, 0.95),
                        }
    results = {}

    for colname, df in zip(colnames, df_3daggr_list): 
        res = {}
        df_3daggr_clean = df_METAggregation.loc[df_METAggregation[colname] != 0.0, ]
        df_3daggr_clean = df_3daggr_clean.reset_index(drop=True)
        df_3daggr_clean[colname] = df_3daggr_clean[colname].astype(float)
        df = df_3daggr_clean.loc[df_3daggr_clean[colname] != 0.0, ]
        df = df.reset_index(drop=True)
        print(f"length of {colname}: " + str(len(df)))
        res['dfstats'] = df[colname].describe()
        res['p1'] = round(df.SUM_LFC3D_neg.quantile(quantile_numbers[colname][0]), 4) # (bottom 10th percentile)
        res['p2'] = round(df.SUM_LFC3D_neg.quantile(quantile_numbers[colname][1]), 4) # (bottom 5th percentile)
        results[colname] = res

    df_3daggr_neg_stats, df_3daggr_pos_stats = results[colnames[0]]['dfstats'], results[colnames[1]]['dfstats']
    NEG_10p_v, NEG_05p_v = results[colnames[0]]['p1'], results[colnames[0]]['p2']
    POS_90p_v, POS_95p_v = results[colnames[1]]['p1'], results[colnames[1]]['p2']

    for colname, df in zip(colnames, df_3daggr_list): 
        arr_LFC3D_discrete, _ = binning(df_METAggregation, 
                                        df_3daggr_neg_stats, df_3daggr_pos_stats, 
                                        NEG_10p_v, POS_90p_v, NEG_05p_v, POS_95p_v, 
                                        colname)
        df_METAggregation[colname+'_dis'] = arr_LFC3D_discrete

    return df_METAggregation

def metaaggregation_displot(
        df_METAggregation, 
        filter_col, x_col, hue_col, 
        filedir, input_gene, out_keyword, 
): 
    """
    Description
        Helper function to plot the distributions for the top 10 and bottom 10 % of points
    """
    
    df_combined_clean = df_METAggregation.loc[df_METAggregation[filter_col] != '-', ]
    df_combined_clean = df_combined_clean.reset_index(drop=True)
    df_combined_clean[x_col] = df_combined_clean[x_col].astype(float)
    plt.figure(figsize=(20, 20), dpi=300)
    sns.displot(df_combined_clean, x=x_col, hue=hue_col, bins=50, palette='tab10')
    plt.title(out_keyword)

    out_name = filedir / f"plots/{input_gene}_Aggr_LFC3D_{out_keyword}_histogram.png"
    plt.savefig(out_name, dpi=300) 

def metaaggregation_scatterplot(
        df_METAggregation, 
        dis_col, pval_col, y_col, 
        filedir, input_gene, out_keyword, pthr, 
): 
    """
    Description
        Helper function to plot scatterplots for the top 10 and bottom 10 % of points
    """

    df_combined_clean = df_METAggregation.loc[df_METAggregation[dis_col] != '-', ]
    df_combined_psig = df_METAggregation.loc[df_METAggregation[pval_col] == 'p<'+str(pthr), ]
    v_combined_psig_SUM_LFC3D_neg_max = max(df_combined_psig[y_col])

    plt.figure(figsize=(12, 8), dpi=300)
    sns.scatterplot(data=df_combined_clean, x="unipos", y=y_col, hue=pval_col, palette='tab10')
    plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    plt.axhline(y = v_combined_psig_SUM_LFC3D_neg_max, color = 'r', linestyle = '--')
    plt.xticks(np.arange(0,len(df_METAggregation), 100))
    plt.title(out_keyword)

    outname = filedir / f"plots/{input_gene}_Aggr_LFC3D_{out_keyword}_dot_per_residue.png"
    plt.savefig(outname, dpi=300)
