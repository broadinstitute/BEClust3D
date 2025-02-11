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

    pthr_str = str(pthr).split('.')[1]
    neg = '_'.join([name, func, score_type, 'neg']).replace('__', '_').strip('_')
    pos = '_'.join([name, func, score_type, 'pos']).replace('__', '_').strip('_')

    # HISTOGRAMS #
    if name == '': 
        histogram_params = [(f'{func}_{score_type}r_neg', neg, 'Negative'), 
                            (f'{func}_{score_type}r_pos', pos, 'Positive'), ]
    else: 
        histogram_params = [(f'{name}_AVG_{score_type}r_neg', neg, 'Negative'), 
                            (f'{name}_AVG_{score_type}r_pos', pos, 'Positive'), ]
    res_neg, res_pos = metaaggregation_histogram(df_Z, histogram_params, 
                                                 edits_filedir / f"plots/{input_gene}_{name}_signal_vs_background.png" )

    if res_neg is None or res_pos is None: return None
    df_Z = binning_lfc3d(df_Z, col=[neg, pos])

    # HISPLOTS #
    hisplots_params = [(f'{neg}_dis', neg, f'{neg}_{pthr_str}_psig', 'Negative P-Value'), 
                       (f'{neg}_dis', neg, f'{neg}_dis', 'Negative P-Value'), 
                       (f'{pos}_dis', pos, f'{pos}_{pthr_str}_psig', 'Positive P-Value'), 
                       (f'{pos}_dis', pos, f'{pos}_dis', 'Positive P-Value'), ]
    metaaggregation_hisplot(df_Z, hisplots_params, 
                            edits_filedir / f"plots/{input_gene}_{name}_{score_type}_histogram.png" )

    # SCATTERPLOT #
    scatterplot_params = [(f'{neg}_dis', f'{neg}_{pthr_str}_psig', neg, 'Negative'), 
                          (f'{pos}_dis', f'{pos}_{pthr_str}_psig', pos, 'Positive')]
    metaaggregation_scatterplot(df_Z, scatterplot_params, input_gene, pthr, 
                                edits_filedir / f"plots/{input_gene}_{name}_{score_type}_scatter.png" )
    
    # Z SCORE SCATTERPLOT #
    scatterplot_params = [(f'{neg}_dis', f'{neg}_dis', f'{neg}_{pthr_str}_z', 'Negative'), 
                          (f'{pos}_dis', f'{pos}_dis', f'{pos}_{pthr_str}_z', 'Positive')]
    metaaggregation_scatterplot(df_Z, scatterplot_params, input_gene, pthr, 
                                edits_filedir / f"plots/{input_gene}_{name}_{score_type}_scatter_colored.png", colors=True )
    