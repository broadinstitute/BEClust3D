import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np

from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp


mut_categories_spaced = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"]
mut_categories_unspaced = [mc.replace(' ' , '_') for mc in mut_categories_spaced]
comparisons1 = [
    f'Nonsense_vs_Missense', f'Nonsense_vs_No_Mutation', f'Nonsense_vs_Silent', 
    f'Splice_Site_vs_Missense', f'Splice_Site_vs_No_Mutation', f'Splice_Site_vs_Silent', 
    f'Missense_vs_Splice_Site_Nonsense', f'No_Mutation_vs_Splice_Site_Nonsense', f'Silent_vs_Splice_Site_Nonsense', 
    f'Missense_vs_Silent_No_Mutation', f'Nonsense_vs_Silent_No_Mutation', f'Splice_Site_vs_Silent_No_Mutation', 
    f'No_Mutation_vs_Silent', 
    f'Nonsense_vs_Splice_Site', 
    f'Splice_Site_Nonsense_vs_Silent_No_Mutation', 
]
comparisons2 = [
    f'Nonsense_vs_No_Mutation', f'Nonsense_vs_Silent', 
    f'Splice_Site_vs_No_Mutation', f'Splice_Site_vs_Silent', 
    f'No_Mutation_vs_Splice_Site_Nonsense', f'Silent_vs_Splice_Site_Nonsense', 
    f'Missense_vs_Silent_No_Mutation', f'Nonsense_vs_Silent_No_Mutation', f'Splice_Site_vs_Silent_No_Mutation', 
    f'No_Mutation_vs_Silent', 
    f'Splice_Site_Nonsense_vs_Silent_No_Mutation', 
]

def add_to_row(
    df1, df2, val_col, function
): 
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
        edits_filedir, screen_names, testtype1, testtype2, 
        hypothesis, partial_col_header='_Splice_Site_Nonsense_vs_Silent_No_Mutation', 
): 
    ### if by gene replace "i, gene in enumerate(unique_genes)"

    fig, axes = plt.subplots(nrows=len(screen_names), ncols=2, sharey=True, 
                             figsize=(12, 4*len(screen_names)))

    # PREP DATAFRAME MW #
    qc_filename = f"qc_validation/{testtype1}_hypothesis{hypothesis}.tsv"
    df_MW_input = pd.read_csv(edits_filedir / qc_filename, sep='\t')
    df_MW_input.replace(-999, pd.NA, inplace=True)  # replace -999 with NaN
    df_MW_input[f"p{partial_col_header}"] = df_MW_input[f"p{partial_col_header}"].apply(negative_log_transformation)

    # PLOT MW #
    if len(screen_names) == 1: 
    # FOR ONE SCREEN #
        plot1 = sns.scatterplot(ax=axes[0], data=df_MW_input[df_MW_input['screenid']==screen_names[0]], 
                                x=f"U{partial_col_header}", y=f"p{partial_col_header}", 
                                hue="gene_name", palette='tab20', s=100, alpha=0.7, edgecolor='k' )
        axes[0].axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05 (-log10 ≈ 1.3)')
        axes[0].axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p = 0.1 (-log10 ≈ 1.0)')

        # LEGEND AND Y AXIS #
        handles, labels = plot1.get_legend_handles_labels()
        axes[0].legend(handles, labels, title="Genes", bbox_to_anchor=(1.0, 1), loc='upper left')
        axes[0].set_ylabel(f'-log10({f"p{partial_col_header}"})')
        axes[0].set_title(f'Mann-Whiteney')
    else: 
    # FOR MULTIPLE SCREEN #
        for i, screen in enumerate(screen_names): 
            plot1 = sns.scatterplot(ax=axes[i,0], data=df_MW_input[df_MW_input['screenid']==screen], 
                                    x=f"U{partial_col_header}", y=f"p{partial_col_header}", 
                                    hue="gene_name", palette='tab20', s=100, alpha=0.7, edgecolor='k' )
            axes[i,0].axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05 (-log10 ≈ 1.3)')
            axes[i,0].axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p = 0.1 (-log10 ≈ 1.0)')

            # LEGEND AND Y AXIS #
            handles, labels = plot1.get_legend_handles_labels()
            axes[i,0].legend(handles, labels, title="Genes", bbox_to_anchor=(1.0, 1), loc='upper left')
            axes[i,0].set_ylabel(f'-log10({f"p{partial_col_header}"})')
        axes[0,0].set_title(f'Mann-Whiteney')

    # PREP DATAFRAME #
    qc_filename = f"qc_validation/{testtype2}_hypothesis{hypothesis}.tsv"
    df_KS_input = pd.read_csv(edits_filedir / qc_filename, sep='\t')
    df_KS_input.replace(-999, pd.NA, inplace=True)  # replace -999 with NaN
    df_KS_input[f"p{partial_col_header}"] = df_KS_input[f"p{partial_col_header}"].apply(negative_log_transformation)

    # PLOTE KS #
    if len(screen_names) == 1: 
    # FOR ONE SCREEN #
        plot2 = sns.scatterplot(ax=axes[1], data=df_KS_input[df_KS_input['screenid']==screen_names[0]], 
                                x=f"D{partial_col_header}", y=f"p{partial_col_header}", 
                                hue="gene_name", palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False )
        axes[1].axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05 (-log10 ≈ 1.3)')
        axes[1].axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p = 0.1 (-log10 ≈ 1.0)')
        axes[1].set_title(f'Kolmogorov-Smirnov')
    else: 
    # FOR MULTIPLE SCREEN #
        for i, screen in enumerate(screen_names): 
            plot2 = sns.scatterplot(ax=axes[i,1], data=df_KS_input[df_KS_input['screenid']==screen], 
                                    x=f"D{partial_col_header}", y=f"p{partial_col_header}", 
                                    hue="gene_name", palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False )
            axes[i,1].axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05 (-log10 ≈ 1.3)')
            axes[i,1].axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p = 0.1 (-log10 ≈ 1.0)')
        axes[0,1].set_title(f'Kolmogorov-Smirnov')

    plt.tight_layout()
    # SAVE PLOT #
    plot_filename = f"plots/hypothesis{hypothesis}_scatterplot.pdf"
    plt.savefig(edits_filedir / plot_filename, dpi=500)


# HYPOTHESIS 1: There is a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations and none (silent/no mutations) per screen, per gene #

def hypothesis_one(
    df_inputs, edits_filedir, 
    screen_names, gene_col, mut_col, val_col, 
    testtype, 
): 
    col_names = ['screenid', 'gene_name']
    col_names.extend([f'n_{mut}' for mut in mut_categories_unspaced])
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in comparisons1 for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in comparisons1 for pref in ('D_', 'p_')])

    # CREATE DF TO STORE TEST RESULTS #
    unique_genes = []
    for df_input in df_inputs: 
        unique = df_input[gene_col].unique().tolist()
        unique_genes = list(set(unique_genes+unique))

    df_output = pd.DataFrame(columns=col_names)
    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes: 

            df_edits = df_input[df_input[gene_col] == current_gene]
            new_row = [screen_name, current_gene]

            # PARSE DF FOR EACH MUT TYPE #
            df_nonsense = df_edits.loc[df_edits[mut_col]=='Nonsense'].reset_index(drop=True)
            df_splice = df_edits.loc[df_edits[mut_col]=='Splice'].reset_index(drop=True)
            df_missense = df_edits.loc[df_edits[mut_col]=='Missense'].reset_index(drop=True)
            df_nomutation = df_edits.loc[df_edits[mut_col]=='NoMutation'].reset_index(drop=True)
            df_silent = df_edits.loc[df_edits[mut_col]=='Silent'].reset_index(drop=True)
            df_splice_nonsense = pd.concat([df_nonsense, df_splice])
            df_silent_nomutation = pd.concat([df_silent, df_nomutation])
            new_row.extend([len(df_missense), len(df_silent), len(df_nonsense), len(df_nomutation), len(df_splice)])

            new_row.extend(add_to_row(df_nonsense, df_missense, val_col, testtype)) # nonsense vs. missense
            new_row.extend(add_to_row(df_nonsense, df_nomutation, val_col, testtype)) # nonsense vs. no mutation
            new_row.extend(add_to_row(df_nonsense, df_silent, val_col, testtype)) # nonsense vs. silent
            new_row.extend(add_to_row(df_splice, df_missense, val_col, testtype)) # splice vs. missense
            new_row.extend(add_to_row(df_splice, df_nomutation, val_col, testtype)) # splice vs. no mutation
            new_row.extend(add_to_row(df_splice, df_silent, val_col, testtype)) # splice vs. silent
            new_row.extend(add_to_row(df_missense, df_splice_nonsense, val_col, testtype)) # missense vs. splice + nonsense
            new_row.extend(add_to_row(df_nomutation, df_splice_nonsense, val_col, testtype)) # no mutation vs. splice + nonsense
            new_row.extend(add_to_row(df_silent, df_splice_nonsense, val_col, testtype)) # silent vs. splice + nonsense
            new_row.extend(add_to_row(df_missense, df_silent_nomutation, val_col, testtype)) # missense vs. silent + no mutation
            new_row.extend(add_to_row(df_nonsense, df_silent_nomutation, val_col, testtype)) # nonsense vs. silent + no mutation
            new_row.extend(add_to_row(df_splice, df_silent_nomutation, val_col, testtype)) # splice vs. silent + no mutation
            new_row.extend(add_to_row(df_nomutation, df_silent, val_col, testtype)) # no mutation vs. silent (negative control)
            new_row.extend(add_to_row(df_nonsense, df_splice, val_col, testtype)) # nonsense vs. splice (negative control)
            new_row.extend(add_to_row(df_splice_nonsense, df_silent_nomutation, val_col, testtype)) # splice + nonsense vs. silent + no mutation

            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            del new_row, df_nonsense, df_splice, df_missense, df_nomutation, df_silent, df_splice_nonsense, df_silent_nomutation

    # SAVE FILE #
    qc_filename = f"qc_validation/{testtype}_hypothesis1.tsv"
    df_output.to_csv(edits_filedir / qc_filename, sep = '\t', index=False)

    return df_output

# HYPOTHESIS 2: There's a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations per gene and none (silent/no mutations) from entire screen #

def hypothesis_two(
    df_inputs, edits_filedir, 
    screen_names, gene_col, mut_col, val_col, 
    testtype, 
): 
    col_names = ['screenid', 'gene_name']
    col_names.extend([f'n_{mut}' for mut in mut_categories_unspaced])
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in comparisons2 for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in comparisons2 for pref in ('D_', 'p_')])

    # CREATE DF TO STORE CONTROL RESULTS #
    unique_genes = []
    for df_input in df_inputs: 
        unique = df_input[gene_col].unique().tolist()
        unique_genes = list(set(unique_genes+unique))

    df_output = pd.DataFrame(columns=col_names)
    df_all_screen_nomutation = pd.DataFrame()
    df_all_screen_silent = pd.DataFrame()

    # GLOBAL SILENT AND NO MUTATION #
    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes:
            df_edits = df_input[df_input[gene_col] == current_gene]

            # PARSE DF FOR EACH MUT TYPE, CONCAT TO PREVIOUS GENE #
            df_nomutation = df_edits.loc[df_edits[mut_col]=='NoMutation'].reset_index(drop=True)
            df_all_screen_nomutation = pd.concat([df_all_screen_nomutation, df_nomutation]).reset_index(drop=True)
            df_silent = df_edits.loc[df_edits[mut_col]=='Silent'].reset_index(drop=True)
            df_all_screen_silent = pd.concat([df_all_screen_silent, df_silent]).reset_index(drop=True)
            del df_nomutation, df_silent

    df_all_silent_nomutation = pd.concat([df_all_screen_nomutation, df_all_screen_silent])

    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes: 

            df_edits = df_input[df_input[gene_col] == current_gene]
            new_row = [screen_name, current_gene]

            # PARSE DF FOR EACH MUT TYPE #
            df_nonsense = df_edits.loc[df_edits[mut_col]=='Nonsense'].reset_index(drop=True)
            df_splice = df_edits.loc[df_edits[mut_col]=='Splice'].reset_index(drop=True)
            df_missense = df_edits.loc[df_edits[mut_col]=='Missense'].reset_index(drop=True)
            df_nomutation = df_edits.loc[df_edits[mut_col]=='NoMutation'].reset_index(drop=True)
            df_silent = df_edits.loc[df_edits[mut_col]=='Silent'].reset_index(drop=True)
            df_splice_nonsense = pd.concat([df_nonsense, df_splice])
            new_row.extend([len(df_missense), len(df_silent), len(df_nonsense), len(df_nomutation), len(df_splice)])

            new_row.extend(add_to_row(df_nonsense, df_all_screen_nomutation, val_col, testtype)) # nonsense vs. no mutation
            new_row.extend(add_to_row(df_nonsense, df_all_screen_silent, val_col, testtype)) # nonsense vs. silent
            new_row.extend(add_to_row(df_splice, df_all_screen_nomutation, val_col, testtype)) # splice vs. no mutation
            new_row.extend(add_to_row(df_splice, df_all_screen_silent, val_col, testtype)) # splice vs. silent
            new_row.extend(add_to_row(df_nomutation, df_splice_nonsense, val_col, testtype)) # no mutation vs. splice + nonsense
            new_row.extend(add_to_row(df_silent, df_splice_nonsense, val_col, testtype)) # silent vs. splice + nonsense
            new_row.extend(add_to_row(df_missense, df_all_silent_nomutation, val_col, testtype)) # missense vs. silent + no mutation
            new_row.extend(add_to_row(df_nonsense, df_all_silent_nomutation, val_col, testtype)) # nonsense vs. silent + no mutation
            new_row.extend(add_to_row(df_splice, df_all_silent_nomutation, val_col, testtype)) # splice vs. silent + no mutation
            new_row.extend(add_to_row(df_nomutation, df_all_screen_silent, val_col, testtype)) # no mutation vs. silent (negative control)
            new_row.extend(add_to_row(df_splice_nonsense, df_all_silent_nomutation, val_col, testtype)) # splice + nonsense vs. silent + no mutation
            
            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            del new_row, df_nonsense, df_splice, df_missense, df_nomutation, df_silent, df_splice_nonsense

    # SAVE FILE #
    qc_filename = f"qc_validation/{testtype}_hypothesis2.tsv"
    df_output.to_csv(edits_filedir / qc_filename, sep = '\t', index=False)

    return df_output
