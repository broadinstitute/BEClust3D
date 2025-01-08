"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-01-01
Description: Plotting Function for Prioritize by Average Split Bin LFC3D / Average Split Bin Metaaggregate

"""

from pathlib import Path
import os
from _average_split_bin_helpers_ import *

def average_split_bin_plots(
        df_Z, workdir, input_gene, pthr=0.05, 
        name='', func='SUM', score_type='LFC3D', 
): 
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    if len(func) == 0: 
        neg = '_'.join([name, f'{score_type}_neg'])
        pos = '_'.join([name, f'{score_type}_pos'])
    elif len(name) == 0: 
        neg = '_'.join([func, f'{score_type}_neg'])
        pos = '_'.join([func, f'{score_type}_pos'])
    else: 
        neg = '_'.join([name, func, f'{score_type}_neg'])
        pos = '_'.join([name, func, f'{score_type}_pos'])

    # HISTOGRAMS #
    histogram_params = [(f'{name}_AVG_{score_type}r_neg'.strip('_'), neg, 'Negative'), 
                        (f'{name}_AVG_{score_type}r_pos'.strip('_'), pos, 'Positive'), ]
    res_neg, res_pos = metaaggregation_histogram(df_Z, histogram_params, edits_filedir, 
                                                 input_gene, name=name)
    if res_neg is None or res_pos is None: 
        return None
    df_Z = binning_lfc3d(df_Z, col=[neg, pos])

    # HISPLOTS #
    hisplots_params = [(f'{neg}_dis', neg, f'{neg}_psig', 'Negative P-Value'), 
                       (f'{neg}_dis', neg, f'{neg}_dis', 'Negative P-Value'), 
                       (f'{pos}_dis', pos, f'{pos}_psig', 'Positive P-Value'), 
                       (f'{pos}_dis', pos, f'{pos}_dis', 'Positive P-Value'), ]
    metaaggregation_hisplot(df_Z, hisplots_params, edits_filedir, input_gene, 
                            name, score_type)

    # SCATTERPLOT #
    scatterplot_params = [(f'{neg}_dis', f'{neg}_psig', neg, 'Negative'), 
                          (f'{pos}_dis', f'{pos}_psig', pos, 'Positive')]
    metaaggregation_scatterplot(df_Z, scatterplot_params, edits_filedir, input_gene, 
                                pthr, name, score_type)
    
    # Z SCORE SCATTERPLOT #
    scatterplot_params = [(f'{neg}_dis', f'{neg}_dis', f'{neg}_z', 'Negative'), 
                          (f'{pos}_dis', f'{pos}_dis', f'{pos}_z', 'Positive')]
    metaaggregation_scatterplot(df_Z, scatterplot_params, edits_filedir, input_gene, 
                                pthr, name, score_type, colors=True)
    