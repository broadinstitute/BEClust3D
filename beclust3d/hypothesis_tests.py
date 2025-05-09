import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import os

from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp


mut_categories_spaced = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"]
mut_categories_unspaced = [mc.replace(' ','_') for mc in mut_categories_spaced]

def hypothesis_tests(
    input_dfs, workdir, 
    screen_names, 
    cases, controls, comp_name, 
    mut_col='Mutation category', val_col='logFC', gene_col='Target Gene Symbol', 
): 
    """
    Description
        Conduct hypothesis 1 (one screen vs control from same screens) and hypothesis 2
        (one screen vs control from all screens) on the set of input screens and genes. 

    Params
        screen_names: list of str, optional
            the names of the input screens

    Returns
        
    """

    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir): 
        os.mkdir(edits_filedir)

    unique_genes = []
    for df in input_dfs: 
        unique = df[gene_col].unique().tolist()
        unique_genes = list(set(unique_genes+unique))

    # AGGREGATE ACROSS SCREENS FOR HYPOTHESIS #
    # MW AND KS TESTS HYPOTHESIS 1 #
    df_MW1_input = hypothesis_one(input_dfs, unique_genes, edits_filedir, cases, controls, comp_name, 
                                  screen_names, gene_col, mut_col, val_col, testtype='MannWhitney')
    df_KS1_input = hypothesis_one(input_dfs, unique_genes, edits_filedir, cases, controls, comp_name, 
                                  screen_names, gene_col, mut_col, val_col, testtype='KolmogorovSmirnov')
    hypothesis_plot(edits_filedir, screen_names, 'screenid', 'gene_name',  # This gives hypothesis1_scatterplot_by_screenid.png (this only works for a single screen but others)
                    testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='1')
    if len(screen_names) > 1: # This gives hypothesis1_scatterplot_by_gene_name.png
        hypothesis_plot(edits_filedir, unique_genes, 'gene_name', 'screenid', 
                        testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='1')
    # MW AND KS TESTS HYPOTHESIS 2 #
    df_MW2_input = hypothesis_two(input_dfs, unique_genes, edits_filedir, cases, controls, comp_name, 
                                  screen_names, gene_col, mut_col, val_col, testtype='MannWhitney')
    df_KS2_input = hypothesis_two(input_dfs, unique_genes, edits_filedir, cases, controls, comp_name, 
                                  screen_names, gene_col, mut_col, val_col, testtype='KolmogorovSmirnov')
    hypothesis_plot(edits_filedir, screen_names, 'screenid', 'gene_name', # This gives hypothesis2_scatterplot_by_screenid.png
                    testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='2')

    if len(screen_names) > 1:  # This gives hypothesis2_scatterplot_by_gene_name.png
        hypothesis_plot(edits_filedir, unique_genes, 'gene_name', 'screenid', 
                        testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='2')
        
    return df_MW1_input, df_MW2_input, df_KS1_input, df_KS2_input


def add_to_row(
    df1, df2, val_col, function, gene_col=False, gene_name=False
): 
    if gene_col and gene_name:
        df1 = df1[df1[gene_col] == gene_name]
        df2 = df2[df2[gene_col] == gene_name]
    else:
        pass

    if len(df1) > 0 and len(df2) > 0: 
        if function == 'KolmogorovSmirnov': 
            U1, p = ks_2samp(df1[val_col], df2[val_col])
            return [U1, p]
        if function == 'MannWhitney': 
            U1, p = mannwhitneyu(df1[val_col], df2[val_col], method="asymptotic")
            return [U1, p]
    return [-999, -999]

def negative_log_transformation(value):
    if pd.notna(value) and value > 0:
        return -np.log10(value)
    return value

def hypothesis_plot(
        edits_filedir, category_names, cat_colname, hue_colname, testtype1, testtype2, hypothesis, 
        partial_col_header='Splice_Site_Nonsense_vs_Silent_No_Mutation', 
): 

    # SETUP PLOT BY NAME (SCREEN or GENE) #
    plt.rcParams.update({'font.size': 8})
    fig, axes = plt.subplots(nrows=len(category_names), ncols=2, sharey=True, 
                             figsize=(12, 5*len(category_names)))

    # PREP DATAFRAME MW #
    qc_filename = f"qc_validation/{testtype1}_hypothesis{hypothesis}.tsv"
    df_MW_input = pd.read_csv(edits_filedir / qc_filename, sep='\t')
    df_MW_input.replace(-999, pd.NA, inplace=True)  # replace -999 with NaN
    df_MW_input[f"p_{partial_col_header}"] = df_MW_input[f"p_{partial_col_header}"].apply(negative_log_transformation)

    all_handles = []
    all_labels = []
    # PLOT TO EXTRACT LEGEND #
    fig_dummy, ax_dummy = plt.subplots()
    scatter_dummy = sns.scatterplot(
        ax=ax_dummy, data=df_MW_input, 
        x=f"U_{partial_col_header}", y=f"p_{partial_col_header}", 
        hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k'
    )
    all_handles, all_labels = ax_dummy.get_legend_handles_labels()
    plt.close(fig_dummy)

    if len(category_names) == 1: axes_list = [axes[0]] # FOR ONE SCREEN #
    else: axes_list = [axes[i,0] for i in range(len(category_names))] # FOR MULTIPLE SCREEN #

    # PLOT MW #
    handles, labels = None, None
    for ax, name in zip(axes_list, category_names):
        plot1 = sns.scatterplot(ax=ax, data=df_MW_input[df_MW_input[cat_colname]==name], 
                                x=f"U_{partial_col_header}", y=f"p_{partial_col_header}", 
                                hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False)
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05 (-log10 ≈ 1.3)')
        ax.axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p = 0.1 (-log10 ≈ 1.0)')

        # GET LEGEND #
        if handles is None and labels is None:
            handles, labels = plot1.get_legend_handles_labels()
        # LEGEND AND Y AXIS #
        ax.set_ylabel(f'-log10({f"p_{partial_col_header}"})')
        ax.set_title(f'Hypothesis {hypothesis}: Mann-Whitney {name}')

        # GRAY BACKGROUND #
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

        # LABELS #
        ax.set_xlabel('MW U-Value')
        ax.set_ylabel('MW -log(P-Value)')

    # PREP DATAFRAME KS #
    qc_filename = f"qc_validation/{testtype2}_hypothesis{hypothesis}.tsv"
    df_KS_input = pd.read_csv(edits_filedir / qc_filename, sep='\t')
    df_KS_input.replace(-999, pd.NA, inplace=True)  # replace -999 with NaN
    df_KS_input[f"p_{partial_col_header}"] = df_KS_input[f"p_{partial_col_header}"].apply(negative_log_transformation)

    if len(category_names) == 1: axes_list = [axes[1]] # FOR ONE SCREEN #
    else: axes_list = [axes[i,1] for i in range(len(category_names))] # FOR MULTIPLE SCREEN #

    for ax, name in zip(axes_list, category_names):
        plot1 = sns.scatterplot(ax=ax, data=df_KS_input[df_KS_input[cat_colname]==name], 
                                x=f"D_{partial_col_header}", y=f"p_{partial_col_header}", 
                                hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False)
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05 (-log10 ≈ 1.3)')
        ax.axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p = 0.1 (-log10 ≈ 1.0)')
        ax.set_title(f'Hypothesis {hypothesis}: Kolmogorov-Smirnov {name}')

        # GRAY BACKGROUND #
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)
        # LABELS #
        ax.set_xlabel('KS D-Value')
        ax.set_ylabel('KS -log(P-Value)')

    # SAVE PLOT #
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    plt.tight_layout()
    plot_filename = f"plots/hypothesis{hypothesis}_scatterplot_by_{cat_colname}.png"
    plt.savefig(edits_filedir / plot_filename, dpi=500)
    plt.close()
    
    # CREATE SEPARATE LEGEND PLOT #
    legend_fig, legend_ax = plt.subplots(figsize=(4, len(all_handles) * 0.3))
    legend_ax.axis('off')
    legend_ax.legend(all_handles, all_labels, title=hue_colname, loc='center', fontsize='small', frameon=False)

    # Save legend separately
    legend_filename = f"plots/hypothesis{hypothesis}_legend_by_{cat_colname}.png"
    legend_fig.savefig(edits_filedir / legend_filename, dpi=500)
    plt.close()

# HYPOTHESIS 1: There is a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations and none (silent/no mutations) per screen, per gene #

def hypothesis_one(
    df_inputs, unique_genes, edits_filedir, 
    cases, controls, comp_name, 
    screen_names, gene_col, mut_col, val_col, 
    testtype, 
): 
    col_names = ['screenid', 'gene_name']
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('D_', 'p_')])
    col_names.extend(['num_of_cases','num_of_controls'])
    df_output = pd.DataFrame(columns=col_names)
    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes: 
            df_edits = df_input[df_input[gene_col] == current_gene]
            new_row = [screen_name, current_gene]
            # PARSE DF FOR EACH MUT TYPE #
            df_case = pd.DataFrame()
            for case in cases: 
                df_case = pd.concat([df_case, df_edits.loc[df_edits[mut_col]==case].reset_index(drop=True)])
            df_control = pd.DataFrame()
            for control in controls: 
                df_control = pd.concat([df_control, df_edits.loc[df_edits[mut_col]==control].reset_index(drop=True)])
            # new_row.extend(add_to_row(df_case, df_control, val_col, testtype))
            new_row.extend(add_to_row(df_case, df_control, val_col, testtype, gene_col, current_gene))
            new_row.extend([len(df_case),len(df_control)])

            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            # df_output.loc[len(df_output),'num_of_cases'] = len(df_case)
            # df_output.loc[len(df_output),'num_of_controls'] = len(df_control)
            del new_row, df_case, df_control

    # SAVE FILE #
    qc_filename = f"qc_validation/{testtype}_hypothesis1.tsv"
    df_output.to_csv(edits_filedir / qc_filename, sep = '\t', index=False)

    return df_output

# HYPOTHESIS 2: There's a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations per gene and none (silent/no mutations) from entire screen #
# the control is not per gene level.
def hypothesis_two(
    df_inputs, unique_genes, edits_filedir, 
    cases, controls, comp_name, 
    screen_names, gene_col, mut_col, val_col, 
    testtype, 
): 
    col_names = ['screenid', 'gene_name']
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('D_', 'p_')])
    col_names.extend(['num_of_cases','num_of_controls'])

    df_output = pd.DataFrame(columns=col_names)
    df_control = pd.DataFrame()

    # GLOBAL SILENT AND NO MUTATION #
    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes:
            df_edits = df_input[df_input[gene_col] == current_gene]

            # PARSE DF FOR EACH MUT TYPE, CONCAT TO PREVIOUS GENE #
            for control in controls: 
                df_control = pd.concat([df_control, df_edits.loc[df_edits[mut_col]==control].reset_index(drop=True)])

    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes: 
            df_edits = df_input[df_input[gene_col] == current_gene]
            new_row = [screen_name, current_gene]

            # PARSE DF FOR EACH MUT TYPE #
            df_case = pd.DataFrame()
            for case in cases: 
                df_case = pd.concat([df_case, df_edits.loc[df_edits[mut_col]==case].reset_index(drop=True)])
            new_row.extend(add_to_row(df_case, df_control, val_col, testtype, gene_col, current_gene))
            new_row.extend([len(df_case),len(df_control)])
            
            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            # df_output.loc[len(df_output),'num_of_cases'] = len(df_case)
            # df_output.loc[len(df_output),'num_of_controls'] = len(df_control)
            del new_row, df_case
    del df_control

    # SAVE FILE #
    qc_filename = f"qc_validation/{testtype}_hypothesis2.tsv"
    df_output.to_csv(edits_filedir / qc_filename, sep = '\t', index=False)

    return df_output
