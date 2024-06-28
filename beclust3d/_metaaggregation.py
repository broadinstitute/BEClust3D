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
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

def metaaggregation(df_LFC_LFC3D,
    workdir, input_gene, input_uniid, structureid, 
    screenid, nRandom, 
): 
    
    # SUM LFC3D #
    
    df_METAggregation = pd.DataFrame()
    df_METAggregation['unipos'] = df_LFC_LFC3D['unipos']
    df_METAggregation['unires'] = df_LFC_LFC3D['unires']
    
    list_sum_LFC3D_neg = []
    list_sum_LFC3D_pos = []
    for i in range(0, len(df_LFC_LFC3D)):
        
        this_res_sum_LFC3D_neg = 0.0
        this_res_sum_LFC3D_pos = 0.0
        header_LFC3D = screenid + "_LFC3D"
        
        LFC3D = df_LFC_LFC3D.at[i, header_LFC3D]
        if (LFC3D != '-') and (float(LFC3D) < 0.0):
            this_res_sum_LFC3D_neg = this_res_sum_LFC3D_neg + float(LFC3D)
        if (LFC3D != '-') and (float(LFC3D) > 0.0):
            this_res_sum_LFC3D_pos = this_res_sum_LFC3D_pos + float(LFC3D)
        
        list_sum_LFC3D_neg.append(this_res_sum_LFC3D_neg)
        list_sum_LFC3D_pos.append(this_res_sum_LFC3D_pos)

    df_METAggregation[sum_colname] = list_sum_LFC3D_neg 
    df_METAggregation['SUM_LFC3D_pos'] = list_sum_LFC3D_pos

    for r in range(0, nRandom):
        list_sum_LFC3D_neg = []
        list_sum_LFC3D_pos = []
        for i in range(0, len(df_LFC_LFC3D)):

            this_res_sum_LFC3D_neg = 0.0
            this_res_sum_LFC3D_pos = 0.0
            header_LFC3D = screenid + "_LFC3Dr" + str(r+1);

            LFC3D = df_LFC_LFC3D.at[i, header_LFC3D]
            if (LFC3D != '-') and (float(LFC3D) < 0.0):
                this_res_sum_LFC3D_neg = this_res_sum_LFC3D_neg + float(LFC3D)
            if (LFC3D != '-') and (float(LFC3D) > 0.0):
                this_res_sum_LFC3D_pos = this_res_sum_LFC3D_pos + float(LFC3D)

            list_sum_LFC3D_neg.append(this_res_sum_LFC3D_neg)
            list_sum_LFC3D_pos.append(this_res_sum_LFC3D_pos)

        col_head_neg = 'SUM_LFC3Dr' + str(r+1) + '_neg'
        col_head_pos = 'SUM_LFC3Dr' + str(r+1) + '_pos'
        df_METAggregation[col_head_neg] = list_sum_LFC3D_neg 
        df_METAggregation[col_head_pos] = list_sum_LFC3D_pos

    # Average of SUM of randomized signal, this is the background noise #

    list_avg_LFC3D_neg = []
    list_avg_LFC3D_pos = []

    for i in range(0, len(df_METAggregation)):
        this_res_sum_LFC3D_neg, this_res_sum_LFC3D_pos = 0.0, 0.0
        this_res_avg_LFC3D_neg, this_res_avg_LFC3D_pos = 0.0, 0.0

        for r in range(0,nRandom):
            col_head_neg = 'SUM_LFC3Dr' + str(r+1) + '_neg'
            LFC3D = df_METAggregation.at[i, col_head_neg]
            this_res_sum_LFC3D_neg = this_res_sum_LFC3D_neg + float(LFC3D)
            if this_res_sum_LFC3D_neg != 0.0:
                this_res_avg_LFC3D_neg = this_res_sum_LFC3D_neg/nRandom
            else:
                this_res_avg_LFC3D_neg = this_res_sum_LFC3D_neg
                
            col_head_pos = 'SUM_LFC3Dr' + str(r+1) + '_pos'
            LFC3D = df_METAggregation.at[i, col_head_pos]
            this_res_sum_LFC3D_pos = this_res_sum_LFC3D_pos + float(LFC3D)
            if this_res_sum_LFC3D_pos != 0.0:
                this_res_avg_LFC3D_pos = this_res_sum_LFC3D_pos/nRandom
            else:
                this_res_avg_LFC3D_pos = this_res_sum_LFC3D_pos
        
        list_avg_LFC3D_neg.append(this_res_avg_LFC3D_neg)
        list_avg_LFC3D_pos.append(this_res_avg_LFC3D_pos)

    df_METAggregation[avg_colname] = list_avg_LFC3D_neg 
    df_METAggregation['AVG_LFC3Dr_pos'] = list_avg_LFC3D_pos

    # CONVERT SIGNAL TO Z SCORE #

    list_LFC3D_neg_zscore, list_LFC3D_neg_pval, list_LFC3D_neg_pval_label = [], [], []
    list_LFC3D_pos_zscore, list_LFC3D_pos_pval, list_LFC3D_pos_pval_label = [], [], []
    lists_LFC3D = [[list_LFC3D_neg_zscore, list_LFC3D_neg_pval, list_LFC3D_neg_pval_label], 
                   [list_LFC3D_pos_zscore, list_LFC3D_pos_pval, list_LFC3D_pos_pval_label]]

    mu_neg = df_METAggregation[avg_colname].mean()
    sigma_neg = df_METAggregation[avg_colname].std()
    mu_pos = df_METAggregation['AVG_LFC3Dr_pos'].mean()
    sigma_pos = df_METAggregation['AVG_LFC3Dr_pos'].std()
    params = [(mu_neg, sigma_neg), (mu_pos, sigma_pos)]

    pthr = 0.001
    pthr_sig_text = 'p<'+str(pthr)
    pthr_insig_text = 'p>='+str(pthr)

    for i in range(0, len(df_METAggregation)):
        # sensitizing and resistant
        for colname, params, lists in zip([sum_colname, 'SUM_LFC3D_pos'], 
                                          params, lists_LFC3D):
            signal = float(df_METAggregation.at[i, colname])
            signal_z = statistics.NormalDist(mu=params[0], sigma=params[1]).zscore(signal)
            signal_p = stats.norm.sf(abs(signal_z))
            if signal_p < pthr:
                signal_plabel = pthr_sig_text
            else:
                signal_plabel = pthr_insig_text

            lists[0].append(signal_z)
            lists[1].append(signal_p)
            lists[2].append(signal_plabel)
        
    df_METAggregation['SUM_LFC3D_neg_z'] = list_LFC3D_neg_zscore 
    df_METAggregation['SUM_LFC3D_neg_p'] = list_LFC3D_neg_pval
    df_METAggregation['SUM_LFC3D_neg_psig'] = list_LFC3D_neg_pval_label
    df_METAggregation['SUM_LFC3D_pos_z'] = list_LFC3D_pos_zscore 
    df_METAggregation['SUM_LFC3D_pos_p'] = list_LFC3D_pos_pval
    df_METAggregation['SUM_LFC3D_pos_psig'] = list_LFC3D_pos_pval_label

    filedir = Path(workdir + input_gene)
    filename = filedir / f"metaaggregation/{input_gene}_{input_uniid}_{structureid}_MetaAggr_LFC3D_and_randomized_background.tsv"
    df_METAggregation.to_csv(filename, "\t", index=False)

    # CALL PLOTS #
    results_sensitizing = metaaggregation_histogram(
        df_METAggregation, 
        filedir, input_gene, 
        avg_colname='AVG_LFC3Dr_neg', sum_colname='SUM_LFC3Dr_neg', out_keyword='sensitizing', 
    )
    results_resistant = metaaggregation_histogram(
        df_METAggregation, 
        filedir, input_gene, 
        avg_colname='AVG_LFC3Dr_pos', sum_colname='SUM_LFC3Dr_pos', out_keyword='resistant', 
    )
    df_METAggregation = binning_lfc3d(df_METAggregation)

    # DISPLOTS #
    metaaggregation_displot(
        df_METAggregation, 
        filter_col='SUM_LFC3D_neg_dis', x_col='SUM_LFC3D_neg', hue_col='SUM_LFC3D_neg_psig', 
        filedir=filedir, input_gene=input_gene, out_keyword='neg_pvalue'
    )
    metaaggregation_displot(
        df_METAggregation, 
        filter_col='SUM_LFC3D_neg_dis', x_col='SUM_LFC3D_neg', hue_col='SUM_LFC3D_neg_dis', 
        filedir=filedir, input_gene=input_gene, out_keyword='neg_dis'
    )
    metaaggregation_displot(
        df_METAggregation, 
        filter_col='SUM_LFC3D_pos_dis', x_col='SUM_LFC3D_pos', hue_col='SUM_LFC3D_pos_psig', 
        filedir=filedir, input_gene=input_gene, out_keyword='pos_pvalue'
    )
    metaaggregation_displot(
        df_METAggregation, 
        filter_col='SUM_LFC3D_pos_dis', x_col='SUM_LFC3D_pos', hue_col='SUM_LFC3D_pos_dis', 
        filedir=filedir, input_gene=input_gene, out_keyword='pos_dis'
    )

    # SCATTERPLOT #
    metaaggregation_scatterplot(
            df_METAggregation, 
            dis_col='SUM_LFC3D_neg_dis', pval_col='SUM_LFC3D_neg_psig', y_col='SUM_LFC3D_neg',
            filedir=filedir, input_gene=input_gene, out_keyword='neg_dis',
    )
    metaaggregation_scatterplot(
            df_METAggregation, 
            dis_col='SUM_LFC3D_pos_dis', pval_col='SUM_LFC3D_pos_psig', y_col='SUM_LFC3D_pos',
            filedir=filedir, input_gene=input_gene, out_keyword='pos_dis',
    )
    
    return df_METAggregation


def metaaggregation_histogram(
        df_METAggregation, 
        filedir, input_gene, 
        avg_colname, sum_colname, out_keyword, 
): 
    
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
    results['sum med'], results['sum std'] = temp.med(), sigma

    z = statistics.NormalDist(mu=mu, sigma=sigma).zscore(-4.6)
    p_cdf = stats.norm.cdf(z)
    p_sf = stats.norm.sf(abs(z))
    results['z'], results['p cdf'], results['p sf'] = z, p_cdf, p_sf

    temp = df_METAggregation_plot[avg_colname]
    results['sum min'], results['sum mean'] = temp.min(), temp.mean()
    results['sum med'], results['sum std'] = temp.med(), temp.std()

    df_METAggregation_plot.plot.area(x='unipos', alpha=0.55, stacked = False)
    plt.legend(loc='lower left', borderaxespad=0)
    plt.xticks(np.arange(0,len(df_METAggregation), 100))
    plt.axhline(y = mu, color = 'r', linestyle = '-')
    plt.axhline(y = mu-sigma, color = 'r', linestyle = '--')

    out_filename = filedir / f"plots/{input_gene}_signal_vs_randomizedbackground_{out_keyword}.png"
    plt.savefig(out_filename, dpi = 300)

    return results


def binning_lfc3d(
        df_METAggregation, 
        colnames = ['SUM_LFC3D_neg', 'SUM_LFC3D_pos'],
): 
    
    df_3daggr_list = [pd.DataFrame() for _ in colnames] # df_3daggr_neg, df_3daggr_pos
    quantile_numbers = {'SUM_LFC3D_neg': (0.1, 0.05), 
                        'SUM_LFC3D_pos': (0.9, 0.95),
                        }
    results = []

    for colname, df in zip(colnames, df_3daggr_list): 
        res = []
        df_3daggr_clean = df_METAggregation.loc[df_METAggregation[colname] != 0.0, ]
        df_3daggr_clean = df_3daggr_clean.reset_index(drop=True)
        df_3daggr_clean[colname] = df_3daggr_clean[colname].astype(float)
        df = df_3daggr_clean.loc[df_3daggr_clean[colname] < 0.0, ]
        df = df.reset_index(drop=True)
        print("length of neg: " + str(len(df)))

        df_stats = df[colname].describe()
        res.append(df_stats)
        p_1 = round(df.SUM_LFC3D_neg.quantile(quantile_numbers[colname][0]), 4) # (bottom 10th percentile)
        res.append(p_1)
        p_2 = round(df.SUM_LFC3D_neg.quantile(quantile_numbers[colname][1]), 4) # (bottom 5th percentile)
        res.append(p_2)
        results.append(res)

    df_3daggr_neg_stats, df_3daggr_pos_stats = res[0][0], res[1][0]
    NEG_10p_v, NEG_05p_v, POS_90p_v, POS_95p_v = res[0][1], res[0][2], res[1][1], res[1][2]

    ### how do i simplify this

    for colname, df in zip(colnames, df_3daggr_list): 
        arr_LFC3D_discrete = []

        for i in range(0, len(df_METAggregation)):
            LFC3D = df_METAggregation.at[i, colname]
            if LFC3D == 0.0:
                LFC3D_discrete = '-'
            else:
                LFC3Df = round(float(LFC3D), 3)
                
                if LFC3Df <= NEG_05p_v:
                    LFC3D_discrete = 'NEG_05p'
                elif (LFC3Df <= NEG_10p_v) and (LFC3Df > NEG_05p_v):
                    LFC3D_discrete = 'NEG_10p'
                elif (LFC3Df <= df_3daggr_neg_stats['25%']) and (LFC3Df > NEG_10p_v):
                    LFC3D_discrete = 'NEG_25p'
                elif (LFC3Df <= df_3daggr_neg_stats['50%']) and (LFC3Df > df_3daggr_neg_stats['25%']):
                    LFC3D_discrete = 'NEG_50p'
                elif (LFC3Df <= df_3daggr_neg_stats['75%']) and (LFC3Df > df_3daggr_neg_stats['50%']):
                    LFC3D_discrete = 'NEG_75p'
                elif (LFC3Df <= df_3daggr_neg_stats['max']) and (LFC3Df > df_3daggr_neg_stats['75%']):
                    LFC3D_discrete = 'NEG_100p'
                elif (LFC3Df >= df_3daggr_pos_stats['min']) and (LFC3Df < df_3daggr_pos_stats['25%']):
                    LFC3D_discrete = 'POS_0p'
                elif (LFC3Df >= df_3daggr_pos_stats['25%']) and (LFC3Df < df_3daggr_pos_stats['50%']):
                    LFC3D_discrete = 'POS_25p'
                elif (LFC3Df >= df_3daggr_pos_stats['50%']) and (LFC3Df < df_3daggr_pos_stats['75%']):
                    LFC3D_discrete = 'POS_50p'
                elif (LFC3Df >= df_3daggr_pos_stats['75%']) and (LFC3Df < POS_90p_v):
                    LFC3D_discrete = 'POS_75p'
                elif (LFC3Df >= POS_90p_v) and (LFC3Df < POS_95p_v):
                    LFC3D_discrete = 'POS_90p'
                elif LFC3Df >= POS_95p_v:
                    LFC3D_discrete = 'POS_95p'
                    
            arr_LFC3D_discrete.append(LFC3D_discrete)
            
        df_METAggregation[colname+'_dis'] = arr_LFC3D_discrete

        return df_METAggregation


def metaaggregation_displot(
        df_METAggregation, 
        filter_col, x_col, hue_col, 
        filedir, input_gene, out_keyword, 
): 
    
    df_combined_clean = df_METAggregation.loc[df_METAggregation[filter_col] != '-', ]
    df_combined_clean = df_combined_clean.reset_index(drop=True)
    df_combined_clean[x_col] = df_combined_clean[x_col].astype(float)

    plt.figure(figsize=(20, 20), dpi=300)
    sns.displot(df_combined_clean, x=x_col, hue=hue_col, bins=50, palette='tab10')

    out_name = filedir / f"plots/{input_gene}_Aggr_LFC3D_{out_keyword}_histogram.png"
    plt.savefig(out_name, dpi = 300) 


def metaaggregation_scatterplot(
        df_METAggregation, 
        dis_col, pval_col, y_col, 
        filedir, input_gene, out_keyword,
): 

    df_combined_clean = df_METAggregation.loc[df_METAggregation[dis_col] != '-', ]
    df_combined_psig = df_METAggregation.loc[df_METAggregation[pval_col] == 'p<0.001', ]
    v_combined_psig_SUM_LFC3D_neg_max = max(df_combined_psig[y_col])

    plt.figure(figsize=(12, 8), dpi=300)
    sns.scatterplot(data=df_combined_clean, x="unipos", y=y_col, hue=pval_col, palette='tab10')
    plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    plt.axhline(y = v_combined_psig_SUM_LFC3D_neg_max, color = 'r', linestyle = '--')
    plt.xticks(np.arange(0,len(df_METAggregation), 100))

    outname = filedir / f"plots/{input_gene}Aggr_LFC3D_{out_keyword}_dot_per_residue.png"
    plt.savefig(outname, dpi = 300)
