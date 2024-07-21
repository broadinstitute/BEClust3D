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
# from warnings import simplefilter
# simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

def metaaggregation(
    df_LFC_LFC3D,
    workdir, 
    input_gene, input_uniprot, structureid, 
    input_screen, 
    nRandom=1000, 
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
        input_uniprot: str, required
            the Uniprot ID for a particular protein
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

    screen_name = input_screen.split('.')[0]
    df_METAggregation = pd.DataFrame()

    df_METAggregation['unipos'] = df_LFC_LFC3D['unipos']
    df_METAggregation['unires'] = df_LFC_LFC3D['unires']
    
    LFC3Ds = []
    list_sum_LFC3D_neg, list_sum_LFC3D_pos = [], []
    for i in range(0, len(df_LFC_LFC3D)):
        this_res_sum_LFC3D_neg, this_res_sum_LFC3D_pos = 0.0, 0.0
        header_LFC3D = f"{screen_name}_LFC3D"
        
        LFC3D = df_LFC_LFC3D.at[i, header_LFC3D]
        LFC3Ds.append(LFC3D)
        if (LFC3D != '-') and (float(LFC3D) < 0.0):
            this_res_sum_LFC3D_neg += float(LFC3D)
        if (LFC3D != '-') and (float(LFC3D) > 0.0):
            this_res_sum_LFC3D_pos += float(LFC3D)
        
        list_sum_LFC3D_neg.append(this_res_sum_LFC3D_neg)
        list_sum_LFC3D_pos.append(this_res_sum_LFC3D_pos)

    df_METAggregation['SUM_LFC3D_neg'] = list_sum_LFC3D_neg 
    df_METAggregation['SUM_LFC3D_pos'] = list_sum_LFC3D_pos

    dict_META = {}
    for r in range(0, nRandom):
        list_sum_LFC3D_neg, list_sum_LFC3D_pos = [], []
        for i in range(0, len(df_LFC_LFC3D)):
            this_res_sum_LFC3D_neg, this_res_sum_LFC3D_pos = 0.0, 0.0
            header_LFC3D = f"{screen_name}_LFC3Dr{str(r+1)}"

            LFC3D = LFC3Ds[i]
            if (LFC3D != '-') and (float(LFC3D) < 0.0):
                this_res_sum_LFC3D_neg = this_res_sum_LFC3D_neg + float(LFC3D)
            if (LFC3D != '-') and (float(LFC3D) > 0.0):
                this_res_sum_LFC3D_pos = this_res_sum_LFC3D_pos + float(LFC3D)

            list_sum_LFC3D_neg.append(this_res_sum_LFC3D_neg)
            list_sum_LFC3D_pos.append(this_res_sum_LFC3D_pos)

        dict_META[f'SUM_LFC3Dr{str(r+1)}_neg'] = list_sum_LFC3D_neg 
        dict_META[f'SUM_LFC3Dr{str(r+1)}_pos'] = list_sum_LFC3D_pos
    df_METAggregation = pd.concat([df_METAggregation, pd.DataFrame(dict_META)], axis=1)

    # Average of SUM of randomized signal, as the background noise #

    list_avg_LFC3D_neg, list_avg_LFC3D_pos = [], []
    for i in range(0, len(df_METAggregation)):
        this_res_sum_LFC3D_neg, this_res_sum_LFC3D_pos = 0.0, 0.0
        this_res_avg_LFC3D_neg, this_res_avg_LFC3D_pos = 0.0, 0.0

        for r in range(0, nRandom):
            LFC3D = df_METAggregation.at[i, f'SUM_LFC3Dr{str(r+1)}_neg']
            this_res_sum_LFC3D_neg += float(LFC3D)
            this_res_avg_LFC3D_neg = this_res_sum_LFC3D_neg / nRandom
                
            LFC3D = df_METAggregation.at[i, f'SUM_LFC3Dr{str(r+1)}_pos']
            this_res_sum_LFC3D_pos += float(LFC3D)
            this_res_avg_LFC3D_pos = this_res_sum_LFC3D_pos / nRandom
        
        list_avg_LFC3D_neg.append(this_res_avg_LFC3D_neg)
        list_avg_LFC3D_pos.append(this_res_avg_LFC3D_pos)

    df_METAggregation['AVG_LFC3Dr_neg'] = list_avg_LFC3D_neg 
    df_METAggregation['AVG_LFC3Dr_pos'] = list_avg_LFC3D_pos

    # convert signal to Z score #

    list_LFC3D_neg_zscore, list_LFC3D_neg_pval, list_LFC3D_neg_pval_label = [], [], []
    list_LFC3D_pos_zscore, list_LFC3D_pos_pval, list_LFC3D_pos_pval_label = [], [], []
    lists_LFC3D = [{'z':list_LFC3D_neg_zscore, 'p':list_LFC3D_neg_pval, 'lab':list_LFC3D_neg_pval_label}, 
                   {'z':list_LFC3D_pos_zscore, 'p':list_LFC3D_pos_pval, 'lab':list_LFC3D_pos_pval_label}
                   ]

    mu_neg = df_METAggregation['AVG_LFC3Dr_neg'].mean()
    sigma_neg = df_METAggregation['AVG_LFC3Dr_neg'].std()
    mu_pos = df_METAggregation['AVG_LFC3Dr_pos'].mean()
    sigma_pos = df_METAggregation['AVG_LFC3Dr_pos'].std()
    params = [{'mu':mu_neg, 's':sigma_neg}, {'mu':mu_pos, 's':sigma_pos}]

    colnames = ['SUM_LFC3D_neg', 'SUM_LFC3D_pos']
    pthr = 0.001

    for i in range(0, len(df_METAggregation)):
        # sensitizing and resistant
        for colname, param, lists in zip(colnames, params, lists_LFC3D): 
            signal = float(df_METAggregation.at[i, colname])
            signal_z = statistics.NormalDist(mu=param['mu'], 
                                             sigma=param['s']).zscore(signal)
            signal_p = stats.norm.sf(abs(signal_z))
            
            if signal_p < pthr: signal_plabel = f'p<{str(pthr)}'
            else:               signal_plabel = f'p>={str(pthr)}'

            lists['z'].append(signal_z)
            lists['p'].append(signal_p)
            lists['lab'].append(signal_plabel)
        
    dict_META = {'SUM_LFC3D_neg_z':list_LFC3D_neg_zscore, 'SUM_LFC3D_neg_p':list_LFC3D_neg_pval, 'SUM_LFC3D_neg_psig':list_LFC3D_neg_pval_label, 
                 'SUM_LFC3D_pos_z':list_LFC3D_pos_zscore, 'SUM_LFC3D_pos_p':list_LFC3D_pos_pval, 'SUM_LFC3D_pos_psig':list_LFC3D_pos_pval_label, }
    df_METAggregation = pd.concat([df_METAggregation, pd.DataFrame(dict_META)], axis=1)

    edits_filedir = Path(workdir + '/' + input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'metaaggregation'):
        os.mkdir(edits_filedir / 'metaaggregation')

    filename = edits_filedir / f"metaaggregation/{structureid}_{screen_name}_MetaAggr_LFC3D_and_randomized_background.tsv"
    df_METAggregation.to_csv(filename, "\t", index=False)

    # CALL PLOTS #
    end_plots = [('SUM_LFC3D_neg', 'AVG_LFC3Dr_neg', 'sensitizing'), 
                 ('SUM_LFC3D_pos', 'AVG_LFC3Dr_pos', 'resistant')]
    for avg, sum, out in end_plots: 
        res = metaaggregation_histogram(
            df_METAggregation, 
            edits_filedir, input_gene, screen_name, 
            avg_colname=avg, sum_colname=sum, out_keyword=out, 
        )
        print(res)
    df_METAggregation = binning_lfc3d(df_METAggregation)
    # DISPLOTS #
    displots = [('SUM_LFC3D_neg_dis', 'SUM_LFC3D_neg', 'SUM_LFC3D_neg_psig', 'neg_pvalue'), 
                ('SUM_LFC3D_neg_dis', 'SUM_LFC3D_neg', 'SUM_LFC3D_neg_dis', 'neg_pvalue'), 
                ('SUM_LFC3D_pos_dis', 'SUM_LFC3D_pos', 'SUM_LFC3D_pos_psig', 'pos_pvalue'), 
                ('SUM_LFC3D_pos_dis', 'SUM_LFC3D_pos', 'SUM_LFC3D_pos_dis', 'pos_pvalue'), 
                ]
    for filter, x, hue, name in displots: 
        metaaggregation_displot(
            df_METAggregation, screen_name, 
            filter_col=filter, x_col=x, hue_col=hue, 
            filedir=edits_filedir, input_gene=input_gene, out_keyword=name
        )
    # SCATTERPLOT #
    scatterplots = [('SUM_LFC3D_neg_dis', 'SUM_LFC3D_neg_psig', 'SUM_LFC3D_neg', 'neg_dis'), 
                    ('SUM_LFC3D_pos_dis', 'SUM_LFC3D_pos_psig', 'SUM_LFC3D_pos', 'pos_dis')]
    for dis, pval, y, out in scatterplots: 
        metaaggregation_scatterplot(
                df_METAggregation, screen_name, 
                dis_col=dis, pval_col=pval, y_col=y,
                filedir=edits_filedir, input_gene=input_gene, out_keyword=out,
        )
    
    return df_METAggregation


def metaaggregation_histogram(
        df_METAggregation, 
        filedir, input_gene, screen_name, 
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

    temp = df_METAggregation_plot[sum_colname]
    mu, sigma = temp.mean(), temp.std()
    results['sum min'], results['sum mean'] = temp.min(), mu
    results['sum med'], results['sum std'] = temp.median(), sigma

    z = statistics.NormalDist(mu=mu, sigma=sigma).zscore(-4.6)
    p_cdf = stats.norm.cdf(z)
    p_sf = stats.norm.sf(abs(z))
    results['z'], results['p cdf'], results['p sf'] = z, p_cdf, p_sf

    temp = df_METAggregation_plot[avg_colname]
    results['sum min'], results['sum mean'] = temp.min(), temp.mean()
    results['sum med'], results['sum std'] = temp.median(), temp.std()

    df_METAggregation_plot.plot.area(x='unipos', alpha=0.55, stacked = False)
    plt.legend(loc='lower left', borderaxespad=0)
    plt.xticks(np.arange(0,len(df_METAggregation), 100))
    plt.axhline(y = mu, color = 'r', linestyle = '-')
    plt.axhline(y = mu-sigma, color = 'r', linestyle = '--')
    plt.title(out_keyword)

    out_filename = filedir / f"plots/{input_gene}_{screen_name}_signal_vs_background_{out_keyword}.png"
    plt.savefig(out_filename, dpi = 300)

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
        screen_name, 
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

    out_name = filedir / f"plots/{input_gene}_{screen_name}_Aggr_LFC3D_{out_keyword}_histogram.png"
    plt.savefig(out_name, dpi=300) 

def metaaggregation_scatterplot(
        df_METAggregation, 
        screen_name, 
        dis_col, pval_col, y_col, 
        filedir, input_gene, out_keyword,
): 
    """
    Description
        Helper function to plot scatterplots for the top 10 and bottom 10 % of points
    """

    df_combined_clean = df_METAggregation.loc[df_METAggregation[dis_col] != '-', ]
    df_combined_psig = df_METAggregation.loc[df_METAggregation[pval_col] == 'p<0.001', ]
    v_combined_psig_SUM_LFC3D_neg_max = max(df_combined_psig[y_col])

    plt.figure(figsize=(12, 8), dpi=300)
    sns.scatterplot(data=df_combined_clean, x="unipos", y=y_col, hue=pval_col, palette='tab10')
    plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    plt.axhline(y = v_combined_psig_SUM_LFC3D_neg_max, color = 'r', linestyle = '--')
    plt.xticks(np.arange(0,len(df_METAggregation), 100))
    plt.title(out_keyword)

    outname = filedir / f"plots/{input_gene}_{screen_name}_Aggr_LFC3D_{out_keyword}_dot_per_residue.png"
    plt.savefig(outname, dpi=300)
