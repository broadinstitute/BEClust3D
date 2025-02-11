"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Helper Functions for Prioritize by Average Split Bin LFC3D / Average Split Bin Metaaggregate

"""

import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
import statistics
import scipy.stats as stats
from scipy.stats import mannwhitneyu

# HELPER FUNCTIONS #

def calculate_stats(signal, param, pthr):
    """
    Helper function to calculate stats: z, p, plabel
    """

    if np.isnan(signal) or param['s'] == 0: ### does this mean we should nromalize across row?
        return 0, 0, 0
    signal_z = statistics.NormalDist(mu=param['mu'], sigma=param['s']).zscore(signal)
    signal_p = stats.norm.sf(abs(signal_z))
    signal_plabel = f'p<{str(pthr)}' if signal_p < pthr and abs(signal) > abs(param['mu']) else f'p>={str(pthr)}'
    return signal_z, signal_p, signal_plabel

def binning_neg_pos(
        df_LFC_LFC3D, df_neg_stats, df_pos_stats, 
        quantile_vals, LFC3D_header
): 
    NEG_10p_v, POS_90p_v, NEG_05p_v, POS_95p_v = quantile_vals
    # binning and weighting 
    arr_LFC3D_disc, arr_LFC3D_weight = [], []

    for i in range(0, len(df_LFC_LFC3D)): 
        LFC3D = df_LFC_LFC3D.at[i, LFC3D_header]
        if LFC3D == '-' or LFC3D == 0.0:
            LFC3D_disc, LFC3D_weight = '-', 0.0
        else: 
            LFC3Df = round(float(LFC3D), 3)
            # ALIGNED FOR BETTER READABILITY #
            if                         LFC3Df <= NEG_05p_v:           LFC3D_disc, LFC3D_weight = 'NEG_05p', -0.95
            elif           NEG_05p_v < LFC3Df <= NEG_10p_v:           LFC3D_disc, LFC3D_weight = 'NEG_10p', -0.9
            elif           NEG_10p_v < LFC3Df <= df_neg_stats['25%']: LFC3D_disc, LFC3D_weight = 'NEG_25p', -0.75
            elif df_neg_stats['25%'] < LFC3Df <= df_neg_stats['50%']: LFC3D_disc, LFC3D_weight = 'NEG_50p', -0.5
            elif df_neg_stats['50%'] < LFC3Df <= df_neg_stats['75%']: LFC3D_disc, LFC3D_weight = 'NEG_75p', -0.25
            elif df_neg_stats['75%'] < LFC3Df <= df_neg_stats['max']: LFC3D_disc, LFC3D_weight = 'NEG_100p', -0.05
            
            elif df_pos_stats['25%'] > LFC3Df >= df_pos_stats['min']: LFC3D_disc, LFC3D_weight = 'POS_0p', 0.05
            elif df_pos_stats['50%'] > LFC3Df >= df_pos_stats['25%']: LFC3D_disc, LFC3D_weight = 'POS_25p', 0.25
            elif df_pos_stats['75%'] > LFC3Df >= df_pos_stats['50%']: LFC3D_disc, LFC3D_weight = 'POS_50p', 0.50
            elif           POS_90p_v > LFC3Df >= df_pos_stats['75%']: LFC3D_disc, LFC3D_weight = 'POS_75p', 0.75
            elif           POS_95p_v > LFC3Df >= POS_90p_v:           LFC3D_disc, LFC3D_weight = 'POS_90p', 0.90
            elif                       LFC3Df >= POS_95p_v:           LFC3D_disc, LFC3D_weight = 'POS_95p', 0.95
            # else: print(LFC3Df, 'Binning Error')
            else: LFC3D_disc, LFC3D_weight = 'NA', 0.0

        arr_LFC3D_disc.append(LFC3D_disc)
        arr_LFC3D_weight.append(LFC3D_weight)

    return arr_LFC3D_disc, arr_LFC3D_weight

def binning_lfc3d(
        df_meta, col, 
): 
    """
    Description
        Helper function to bin the top 10 and bottom 10 % of points
    """
    
    df_3d_list = [pd.DataFrame() for _ in col]
    quantile_numbers = {col[0]: (0.1, 0.05), col[1]: (0.9, 0.95)}
    result = {}

    for colname, df in zip(col, df_3d_list): 
        res = {}
        df_3d_clean = df_meta.loc[df_meta[colname] != 0.0, ].reset_index(drop=True)
        df[colname] = df_3d_clean[colname].replace('-', np.nan).astype(float)
        # df = df_3d_clean.loc[df_3d_clean[colname] != 0.0, ].reset_index(drop=True)

        res['dfstats'] = df[colname].describe()
        res['p1'] = round(df[colname].quantile(quantile_numbers[colname][0]), 4) # (bottom 10th percentile)
        res['p2'] = round(df[colname].quantile(quantile_numbers[colname][1]), 4) # (bottom 5th percentile)
        result[colname] = res

    df_neg_stats, df_pos_stats = result[col[0]]['dfstats'], result[col[1]]['dfstats']
    bins = [result[col[0]]['p1'], result[col[1]]['p1'], result[col[0]]['p2'], result[col[1]]['p2']]

    for colname, df in zip(col, df_3d_list): 
        arr_LFC3D_disc, _ = binning_neg_pos(df_meta, df_neg_stats, df_pos_stats, bins, colname)
        df_meta[colname+'_dis'] = arr_LFC3D_disc

    return df_meta

def metaaggregation_histogram(
        df_meta, params, out_filename, 
): 
    """
    Description
        Helper function to plot histograms of the values along the length of the gene
    """
    
    fig, ax = plt.subplots(1, 2, figsize=(16, 6), dpi=300)
    results_list = []

    for i, (avg, sum, out) in enumerate(params): 
        res = {}
        df_meta_plot = pd.DataFrame()
        df_meta_plot['unipos'] = df_meta['unipos']
        df_meta_plot[sum] = df_meta[sum].replace('-', np.nan).astype(float) ### can we fix the default format 250120
        df_meta_plot[avg] = df_meta[avg].replace('-', np.nan).astype(float) ### can we fix the default format 250120
        # df_meta_plot_sum = df_meta_plot[sum].dropna().tolist()
        # df_meta_plot_avg = df_meta_plot[avg].dropna().tolist()
        df_meta_filtered = df_meta_plot.dropna(subset=[sum, avg])
        df_meta_plot_sum = df_meta_filtered[sum].tolist()
        df_meta_plot_avg = df_meta_filtered[avg].tolist()

        U1, p = mannwhitneyu(df_meta_plot[sum], df_meta_plot[avg], method="asymptotic" )
        res['mannwhitneyu U1'], res['mannwhitneyu p'] = U1, p
        r, p = stats.pearsonr(df_meta_plot_sum, df_meta_plot_avg )
        res['pearsonr r'], res['pearsonr p'] = r, p
        # SUM #
        res['sum min'], res['sum mean'] = df_meta_plot[sum].min(), df_meta_plot[sum].mean()
        res['sum med'], res['sum std'] = df_meta_plot[sum].median(), df_meta_plot[sum].std()

        if res['sum std'] == 0: 
            return None
        z = statistics.NormalDist(mu=res['sum mean'], sigma=res['sum std']).zscore(-4.6)
        res['z'], res['p cdf'], res['p sf'] = z, stats.norm.cdf(z), stats.norm.sf(abs(z))
        # AVG #
        res['avg min'], res['avg mean'] = df_meta_plot[avg].min(), df_meta_plot[avg].mean()
        res['avg med'], res['avg std'] = df_meta_plot[avg].median(), df_meta_plot[avg].std()

        # PLOT #
        df_meta_plot.plot.area(x='unipos', alpha=0.55, stacked = False, ax=ax[i])
        ax[i].axhline(y = res['sum mean'], color = 'r', linestyle = '-')
        ax[i].axhline(y = res['sum mean']-res['sum std'], color = 'r', linestyle = '--')

        ax[i].legend(loc='lower left', borderaxespad=0)
        ax[i].set_xticks(np.arange(0,len(df_meta), 100))
        ax[i].set_title(out)

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

        del df_meta_plot
        results_list.append(res)
    
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(out_filename, dpi=300)
    return results_list[0], results_list[1]

def metaaggregation_hisplot(
        df_meta, params, out_name
): 
    """
    Description
        Helper function to plot the distributions for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, 4, figsize=(36, 6), dpi=300)

    for i, (dis, x, hue, name) in enumerate(params): 
    
        df_combined_clean = df_meta.loc[df_meta[dis] != '-', ].reset_index(drop=True)
        df_combined_clean[x] = df_combined_clean[x].astype(float)
        sns.histplot(df_combined_clean, x=x, hue=hue, bins=50, palette='tab10', ax=ax[i])
        ax[i].set_title(name)

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

    plt.subplots_adjust(wspace=0.3)
    plt.savefig(out_name, dpi=300) 

def metaaggregation_scatterplot(
        df_meta, params, input_gene, pthr, 
        outname, colors=False, 
): 
    """
    Description
        Helper function to plot scatterplots for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, 2, figsize=(18, 6), dpi=300)

    for i, (dis, pval, y, out) in enumerate(params): 

        df_combined_clean = df_meta.loc[df_meta[dis] != '-', ]
        if colors: # MULTIPLE COLORS
            if 'pos' in dis: 
                ax[i].axhline(y = 1.65, color = 'r', linestyle = '--')
                ax[i].axhline(y = 1.96, color = 'r', linestyle = '--')
                ax[i].axhline(y = 2.58, color = 'r', linestyle = '--')
            if 'neg' in dis: 
                ax[i].axhline(y = -1.65, color = 'r', linestyle = '--')
                ax[i].axhline(y = -1.96, color = 'r', linestyle = '--')
                ax[i].axhline(y = -2.58, color = 'r', linestyle = '--')
            sns.scatterplot(data=df_combined_clean, x="unipos", y=y, hue=pval, palette='tab10', ax=ax[i])

        else: # 2 COLORS #
            df_combined_psig = df_meta.loc[df_meta[pval] == 'p>='+str(pthr), ]
            if 'pos' in dis: 
                line_val = df_combined_psig[y][df_combined_psig[y] != '-'].astype(float).max()
            if 'neg' in dis: 
                line_val = df_combined_psig[y][df_combined_psig[y] != '-'].astype(float).min()
            ax[i].axhline(y = line_val, color = 'r', linestyle = '--')
            sns.scatterplot(data=df_combined_clean, x="unipos", y=y, hue=pval, palette='tab10', ax=ax[i])

        ax[i].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
        ax[i].set_xticks(np.arange(0, len(df_meta), 100))
        ax[i].set_title(f"{input_gene} {out}")

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

    plt.subplots_adjust(wspace=0.3)
    plt.savefig(outname, dpi=300)
