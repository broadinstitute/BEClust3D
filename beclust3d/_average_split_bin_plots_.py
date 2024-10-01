import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
import statistics
import scipy.stats as stats
from scipy.stats import mannwhitneyu

def binning(
        df_LFC_LFC3D, df_LFC3D_neg_stats, df_LFC3D_pos_stats, 
        quantile_vals, LFC3D_header
):
    NEG_10p_v, POS_90p_v, NEG_05p_v, POS_95p_v = quantile_vals
    # binning and weighting 
    arr_LFC3D_discrete = []
    arr_LFC3D_weight = []

    for i in range(0, len(df_LFC_LFC3D)): 
        LFC3D = df_LFC_LFC3D.at[i, LFC3D_header]
        if LFC3D == '-' or LFC3D == 0.0:
            LFC3D_discrete = '-'
            LFC3D_weight = 0.0
        else: 
            LFC3Df = round(float(LFC3D), 3)
            # aligned for better readability
            if                               LFC3Df <= NEG_05p_v:                 LFC3D_discrete, LFC3D_weight = 'NEG_05p', -0.95
            elif                 NEG_05p_v < LFC3Df <= NEG_10p_v:                 LFC3D_discrete, LFC3D_weight = 'NEG_10p', -0.9
            elif                 NEG_10p_v < LFC3Df <= df_LFC3D_neg_stats['25%']: LFC3D_discrete, LFC3D_weight = 'NEG_25p', -0.75
            elif df_LFC3D_neg_stats['25%'] < LFC3Df <= df_LFC3D_neg_stats['50%']: LFC3D_discrete, LFC3D_weight = 'NEG_50p', -0.5
            elif df_LFC3D_neg_stats['50%'] < LFC3Df <= df_LFC3D_neg_stats['75%']: LFC3D_discrete, LFC3D_weight = 'NEG_75p', -0.25
            elif df_LFC3D_neg_stats['75%'] < LFC3Df <= df_LFC3D_neg_stats['max']: LFC3D_discrete, LFC3D_weight = 'NEG_100p', -0.05
            elif df_LFC3D_pos_stats['25%'] > LFC3Df >= df_LFC3D_pos_stats['min']: LFC3D_discrete, LFC3D_weight = 'POS_0p', 0.05
            elif df_LFC3D_pos_stats['50%'] > LFC3Df >= df_LFC3D_pos_stats['25%']: LFC3D_discrete, LFC3D_weight = 'POS_25p', 0.25
            elif df_LFC3D_pos_stats['75%'] > LFC3Df >= df_LFC3D_pos_stats['50%']: LFC3D_discrete, LFC3D_weight = 'POS_50p', 0.50
            elif                 POS_90p_v > LFC3Df >= df_LFC3D_pos_stats['75%']: LFC3D_discrete, LFC3D_weight = 'POS_75p', 0.75
            elif                 POS_95p_v > LFC3Df >= POS_90p_v:                 LFC3D_discrete, LFC3D_weight = 'POS_90p', 0.90
            elif                             LFC3Df >= POS_95p_v:                 LFC3D_discrete, LFC3D_weight = 'POS_95p', 0.95

        arr_LFC3D_discrete.append(LFC3D_discrete)
        arr_LFC3D_weight.append(LFC3D_weight)

    return arr_LFC3D_discrete, arr_LFC3D_weight

def binning_lfc3d(
        df_meta, colnames, 
): 
    """
    Description
        Helper function to bin the top 10 and bottom 10 % of points
    """
    
    df_3daggr_list = [pd.DataFrame() for _ in colnames]
    quantile_numbers = {colnames[0]: (0.1, 0.05), 
                        colnames[1]: (0.9, 0.95), }
    results = {}

    for colname, df in zip(colnames, df_3daggr_list): 
        res = {}
        df_3daggr_clean = df_meta.loc[df_meta[colname] != 0.0, ]
        df_3daggr_clean = df_3daggr_clean.reset_index(drop=True)
        df_3daggr_clean[colname] = df_3daggr_clean[colname].astype(float)
        df = df_3daggr_clean.loc[df_3daggr_clean[colname] != 0.0, ]
        df = df.reset_index(drop=True)
        # print(f"length of {colname}: " + str(len(df)))
        res['dfstats'] = df[colname].describe()
        res['p1'] = round(df[colname].quantile(quantile_numbers[colname][0]), 4) # (bottom 10th percentile)
        res['p2'] = round(df[colname].quantile(quantile_numbers[colname][1]), 4) # (bottom 5th percentile)
        results[colname] = res

    df_3daggr_neg_stats, df_3daggr_pos_stats = results[colnames[0]]['dfstats'], results[colnames[1]]['dfstats']
    NEG_10p_v, NEG_05p_v = results[colnames[0]]['p1'], results[colnames[0]]['p2']
    POS_90p_v, POS_95p_v = results[colnames[1]]['p1'], results[colnames[1]]['p2']

    for colname, df in zip(colnames, df_3daggr_list): 
        arr_LFC3D_discrete, _ = binning(df_meta, df_3daggr_neg_stats, df_3daggr_pos_stats, 
                                        [NEG_10p_v, POS_90p_v, NEG_05p_v, POS_95p_v], colname)
        df_meta[colname+'_dis'] = arr_LFC3D_discrete

    return df_meta

def LFC3D_plots(
        df_Z, edits_filedir, input_gene, pthr, name='', 
): 
    # HISTOGRAMS #
    histogram_params = [(name+'AVG_LFC3Dr_neg', name+'SUM_LFC3D_neg', 'sensitizing'), 
                        (name+'AVG_LFC3Dr_pos', name+'SUM_LFC3D_pos', 'resistant'), ]
    res_neg, res_pos = metaaggregation_histogram(
        df_Z, histogram_params, edits_filedir, input_gene, )
    if res_neg is None or res_pos is None: 
        return None

    df_Z = binning_lfc3d(df_Z, colnames=[name+'SUM_LFC3D_neg',name+'SUM_LFC3D_pos'],)

    # HISPLOTS #
    hisplots_params = [(name+'SUM_LFC3D_neg_dis', name+'SUM_LFC3D_neg', name+'SUM_LFC3D_neg_psig', 'neg_pvalue'), 
                       (name+'SUM_LFC3D_neg_dis', name+'SUM_LFC3D_neg', name+'SUM_LFC3D_neg_dis', 'neg_pvalue'), 
                       (name+'SUM_LFC3D_pos_dis', name+'SUM_LFC3D_pos', name+'SUM_LFC3D_pos_psig', 'pos_pvalue'), 
                       (name+'SUM_LFC3D_pos_dis', name+'SUM_LFC3D_pos', name+'SUM_LFC3D_pos_dis', 'pos_pvalue'), ]
    metaaggregation_hisplot(
        df_Z, hisplots_params, filedir=edits_filedir, input_gene=input_gene, )

    # SCATTERPLOT #
    scatterplot_params = [(name+'SUM_LFC3D_neg_dis', name+'SUM_LFC3D_neg_psig', name+'SUM_LFC3D_neg', 'neg_dis'), 
                          (name+'SUM_LFC3D_pos_dis', name+'SUM_LFC3D_pos_psig', name+'SUM_LFC3D_pos', 'pos_dis')]
    metaaggregation_scatterplot(
        df_Z, scatterplot_params, filedir=edits_filedir, input_gene=input_gene, pthr=pthr, )


def metaaggregation_histogram(
        df_meta, params, filedir, input_gene, 
): 
    """
    Description
        Helper function to plot histograms of the values along the length of the gene
    """
    
    fig, ax = plt.subplots(1, 2, figsize=(16, 6), dpi=300)
    results_list = []

    # avg_colname='AVG_LFC3Dr_neg', sum_colname='SUM_LFC3D_neg', out_keyword='sensitizing', 
    for i, (avg, sum, out) in enumerate(params): 

        results = {}
        df_meta_plot = pd.DataFrame()
        df_meta_plot['unipos'] = df_meta['unipos']
        df_meta_plot[sum] = df_meta[sum]
        df_meta_plot[avg] = df_meta[avg]

        U1, p = mannwhitneyu(df_meta_plot[sum], df_meta_plot[avg], method="asymptotic" )
        results['mannwhitneyu U1'], results['mannwhitneyu p'] = U1, p
        r, p = stats.pearsonr(df_meta_plot[sum], df_meta_plot[avg] )
        results['pearsonr r'], results['pearsonr p'] = r, p

        # SUM #
        temp = df_meta_plot[sum]
        results['sum min'], results['sum mean'] = temp.min(), temp.mean()
        results['sum med'], results['sum std'] = temp.median(), temp.std()

        if results['sum std'] == 0: 
            return None
        z = statistics.NormalDist(mu=results['sum mean'], sigma=results['sum std']).zscore(-4.6)
        results['z'], results['p cdf'], results['p sf'] = z, stats.norm.cdf(z), stats.norm.sf(abs(z))

        # AVG #
        temp = df_meta_plot[avg]
        results['avg min'], results['avg mean'] = temp.min(), temp.mean()
        results['avg med'], results['avg std'] = temp.median(), temp.std()

        # PLOT #
        df_meta_plot.plot.area(x='unipos', alpha=0.55, stacked = False, ax=ax[i])
        ax[i].axhline(y = results['sum mean'], color = 'r', linestyle = '-')
        ax[i].axhline(y = results['sum mean']-results['sum std'], color = 'r', linestyle = '--')

        ax[i].legend(loc='lower left', borderaxespad=0)
        ax[i].set_xticks(np.arange(0,len(df_meta), 100))
        ax[i].set_title(out)

        del temp, df_meta_plot
        results_list.append(results)

    out_filename = filedir / f"plots/{input_gene}_signal_vs_background.png"
    plt.savefig(out_filename, dpi=300)
    return results_list[0], results_list[1]

def metaaggregation_hisplot(
        df_meta, 
        params, filedir, input_gene, 
): 
    """
    Description
        Helper function to plot the distributions for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, 4, figsize=(36, 6), dpi=300)

    for i, (filter, x, hue, name) in enumerate(params): 
    
        df_combined_clean = df_meta.loc[df_meta[filter] != '-', ].reset_index(drop=True)
        df_combined_clean[x] = df_combined_clean[x].astype(float)
        sns.histplot(df_combined_clean, x=x, hue=hue, bins=50, palette='tab10', ax=ax[i])
        ax[i].set_title(name)

    out_name = filedir / f"plots/{input_gene}_Aggr_LFC3D_histogram.png"
    plt.savefig(out_name, dpi=300) 

def metaaggregation_scatterplot(
        df_meta, params, filedir, input_gene, pthr, 
): 
    """
    Description
        Helper function to plot scatterplots for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, 2, figsize=(18, 6), dpi=300)

    for i, (dis, pval, y, out) in enumerate(params): 

        df_combined_clean = df_meta.loc[df_meta[dis] != '-', ]
        df_combined_psig = df_meta.loc[df_meta[pval] == 'p<'+str(pthr), ]
        if 'pos' in dis: 
            v_combined_psig_SUM_LFC3D = min(df_combined_psig[y])
        if 'neg' in dis: 
            v_combined_psig_SUM_LFC3D = max(df_combined_psig[y])

        sns.scatterplot(data=df_combined_clean, x="unipos", y=y, hue=pval, palette='tab10', ax=ax[i])
        ax[i].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
        ax[i].axhline(y = v_combined_psig_SUM_LFC3D, color = 'r', linestyle = '--')
        ax[i].set_xticks(np.arange(0, len(df_meta), 100))
        ax[i].set_title(out)

    outname = filedir / f"plots/{input_gene}_Aggr_LFC3D_dot_per_residue.png"
    plt.savefig(outname, dpi=300)
