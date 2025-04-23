import pandas as pd
import os, sys
import warnings

def main(**kwargs):
    workdir=kwargs['workdir'] 
    screen_dir=kwargs['screen_dir']
    screens=kwargs['screens']
    input_gene=kwargs['input_gene']
    input_uniprot=kwargs['input_uniprot']
    output_dir=kwargs['output_dir']
    nRandom=kwargs['nRandom']
    pthr_LFC=kwargs['pthr_LFC']
    pthr_LFC3D=kwargs['pthr_LFC3D']
    mut_list_col=kwargs['mut_list_col']
    mut_col=kwargs['mut_col']
    val_col=kwargs['val_col']
    gene_col=kwargs['gene_col']
    edits_col=kwargs['edits_col']
    gRNA_col=kwargs['gRNA_col']
    input_pdb=kwargs['input_pdb']
    input_fasta=kwargs['input_fasta']
    assign_hierachy=kwargs['assign_hierachy']
    hierachy=kwargs['hierachy']
    
    structureid = f"AF-{input_uniprot}-F1-model_v4"   
    screens = screens.split(',')    
    os.makedirs(output_dir, exist_ok=True)

    hierachy_list = [x.strip() for x in hierachy.split(',')]
    mutation_hierarchy = {mutation: idx for idx, mutation in enumerate(hierachy_list)}

    # Function to determine the Mutation_Category based on Mutation_Type
    def prioritize_mutations(mutation_types):
        if not mutation_types or pd.isna(mutation_types):
            return "No Mutation"
        mutation_list = mutation_types.split(', ')

        # Replace "Splice-donor" and "Splice-acceptor" with "Splice Site"
        mutation_list = ["Splice Site" if "Splice" in mutation else mutation for mutation in mutation_list]
        # Filter out mutation types not in the hierarchy list
        mutation_list = [mutation for mutation in mutation_list if mutation in mutation_hierarchy]

        # If no recognized mutation types are found, return an empty string
        if not mutation_list:
            return ""
        # Sort the mutations based on the defined hierarchy
        mutation_list.sort(key=lambda x: mutation_hierarchy.get(x, float('inf')))

        return mutation_list[0]

    def count_commas(s):
        if pd.isna(s):
            return 0
        return s.count(',')

    unique_genes_out = []
    # CHECK INPUT #
    for i, screen in enumerate(screens):
        df_screen = pd.read_csv(f'{screen_dir}/{screen}', sep='\t')
        assert all(x in df_screen.columns for x in [mut_col, val_col, gene_col, edits_col]), "Input column(s) are not in the dataframe"

        # CHECK VAL_COL a list of values
        assert all(isinstance(x, float) for x in df_screen[val_col]), "val_col is not an float"

        # CHECK GENE_COL a list of gene names
        unique_genes = df_screen[gene_col].unique()
        unique_genes.sort()
        unique_genes_out.append(' '.join(unique_genes))
        # print(f'All genes present in screen {screen}', unique_genes)
        assert all(isinstance(x, str) for x in df_screen[gene_col]), "gene_col is not a string"

        # CHECK gRNA_COL
        # assert all(isinstance(x, str) for x in df_screen[gRNA_col]), "gRNA_col is not a string"
        # assert all(len(x)==20 for x in df_screen[gRNA_col]), "gRNA_col should be 20 bps"

        # CHECK MUT_LIST_COL and EDITS_COL MATCH
        if len(mut_list_col) > 0:
            mismatch_rows = df_screen[df_screen[mut_list_col].apply(lambda x: count_commas(x)) != df_screen[edits_col].apply(lambda x: count_commas(x))]
            assert mismatch_rows.empty, "MUT_LIST_COL and EDITS_COL do not match"
            if not mismatch_rows.empty:
                print('Check the following rows')
                print(df_screen[mismatch_rows])

        # PRIORITIZE MUT_COL BASED ON MUT_LIST_COL
        if assign_hierachy:
            df_screen[mut_col] = df_screen[mut_list_col].apply(lambda x: prioritize_mutations(x))
            new_screen = screen.split('.')[0]+'_prioritized.'+screen.split('.')[1]
            print(screen, 'renamed to', new_screen)
            df_screen.to_csv(f'{workdir}rawdata/{new_screen}', sep='\t', index=False)
            screens[i] = new_screen

    screen_names = [screen.split('.')[0] for screen in screens]

    ## Prepare AlphaFold 3D structure, FASTA sequence and DSSP result
    from BEClust3D.beclust3d.af_structural_features import af_structural_features

    af_structural_features(
        workdir=output_dir,
        input_gene=input_gene,
        input_uniprot=input_uniprot,
        structureid=structureid,
        # optional
        user_uniprot=input_fasta,
        user_pdb=input_pdb,
        # user_dssp=user_dssp,
        )

    ## Preprocess screen data and QC plots
    from BEClust3D.beclust3d.preprocess_be_results import parse_base_editing_results, plot_base_editing_results
    from BEClust3D.beclust3d.hypothesis_tests import hypothesis_tests
    warnings.filterwarnings('ignore')

    df_struc = pd.read_csv(f"{output_dir}/{structureid}_coord_struc_features.tsv", sep = "\t")
    df_inputs = [pd.read_csv(os.path.join(screen_dir,f'{screen}'), sep='\t') for screen in screens]

    _ = parse_base_editing_results(
        input_dfs = df_inputs,
        workdir = output_dir, input_gene = input_gene, screen_names = screen_names,
        mut_col=mut_col, val_col=val_col, gene_col=gene_col, edits_col=edits_col,
        split_char=', ',
        )

    _ = plot_base_editing_results(
        input_dfs = df_inputs,
        workdir = output_dir, input_gene = input_gene, screen_names = screen_names,
        mut_col=mut_col, val_col=val_col, gene_col=gene_col,
        )

    hypothesis_tests(
        input_dfs = df_inputs,
        workdir = output_dir,
        screen_names=screen_names,
        cases=['Splice','Nonsense'], controls=['No Mutation','Silent'], comp_name='Splice_Site_Nonsense_vs_Silent_No_Mutation',
        mut_col=mut_col, val_col=val_col, gene_col=gene_col,
    )
    
    ## Randomize screen data
    from BEClust3D.beclust3d.randomize_be_results import randomize_be_results

    df_missenses = [pd.read_csv(f'{output_dir}/screendata/{input_gene}_{screen_name}_Missense.tsv', sep='\t') for screen_name in screen_names]
    for df, screen in zip(df_missenses, screen_names):
        randomize_be_results(
            df_missense=df,
            workdir=output_dir,
            input_gene=input_gene,
            screen_name=screen,
            nRandom=nRandom,
            seed=True
            )

    ## Prioritize 
    from BEClust3D.beclust3d.prioritize_by_sequence import prioritize_by_sequence, plots_by_sequence
    from BEClust3D.beclust3d.randomize_by_sequence import randomize_by_sequence

    for screen in screen_names:
        file_dict = {
            'Missense' : pd.read_csv(f"{output_dir}/screendata/{input_gene}_{screen}_Missense.tsv", sep='\t'),
            'Silent' : pd.read_csv(f"{output_dir}/screendata/{input_gene}_{screen}_Silent.tsv", sep='\t'),
            'Nonsense' : pd.read_csv(f"{output_dir}/screendata/{input_gene}_{screen}_Nonsense.tsv", sep='\t'),
        }
        df_nomutation = pd.read_csv(f"{output_dir}/screendata/{input_gene}_{screen}_No_Mutation.tsv", sep = "\t")
        
        out_df = prioritize_by_sequence(
                df_struc     =df_struc,
                df_consrv    =None,
                df_nomutation=df_nomutation,
                file_dict=file_dict,
                workdir      =output_dir,
                input_gene   =input_gene,
                screen_name  =screen,
                function_name='mean',
            )

        plots_by_sequence(
            out_df, output_dir,
            input_gene, screen,
            function_name='mean',
        )

        df_missense = pd.read_csv(f"{output_dir}/screendata/{input_gene}_{screen}_proteinedits.tsv", sep='\t')
        df_rand = pd.read_csv(f"{output_dir}/randomized_screendata/{input_gene}_{screen}_Missense_rand.tsv", sep='\t')
        
        randomize_by_sequence(
            workdir      =output_dir,
            df_missense=df_missense,
            df_rand=df_rand,
            input_gene   =input_gene,
            screen_name  =screen,
            nRandom      =nRandom,
            conservation =False,
        )
    ## Calculate BEClust3D
    from BEClust3D.beclust3d.calculate_lfc3d import calculate_lfc3d

    df_str_cons = pd.read_csv(f"{output_dir}/{structureid}_coord_struc_features.tsv", sep = "\t")
    str_cons_dfs,str_cons_rand_dfs = [],[]
    for screen in screen_names:
        str_cons_dfs.append(pd.read_csv(f"{output_dir}/screendata/{input_gene}_{screen}_proteinedits.tsv", sep='\t'))
        str_cons_rand_dfs.append(pd.read_csv(f"{output_dir}/randomized_screendata/{input_gene}_{screen}_Missense_proteinedits_rand.tsv", sep='\t'))

    calculate_lfc3d(
        df_str_cons  =df_str_cons,
        workdir      =output_dir,
        input_gene   =input_gene,
        screen_names =screen_names,
        nRandom      =nRandom,
        str_cons_dfs=str_cons_dfs, str_cons_rand_dfs=str_cons_rand_dfs,
    )

    from BEClust3D.beclust3d.average_split_bin_lfc3d import average_split_bin
    from BEClust3D.beclust3d.annotate_spatial_clusters import clustering_union
    from BEClust3D.beclust3d.characterization_plots import lfc_vs_lfc3d_scatterplot
    from BEClust3D.beclust3d._average_split_bin_helpers_ import scatterplot_by_residue_LFC3D
    df_LFC_LFC3D_rand_tsv = os.path.join(workdir,'LFC3D',f'{input_gene}_LFC_LFC3D_LFC3Dr.tsv')
    df_LFC_LFC3D_rand = pd.read_csv(df_LFC_LFC3D_rand_tsv, sep='\t')

    # This is only needed in Single Screen LFC3D and Meta-aggregation LFC/LFC3D.
    average_split_bin( # only for LFC3D
        df_LFC_LFC3D_rand, 
        workdir=workdir, input_gene=input_gene, screen_names=screen_names, 
        pthr=pthr_LFC3D, score_type='LFC3D', # pthr_LFC3D will be used for identifying _psig 
    )

    ## Screen-level Clustering for LFC and LFC3D
    for each_screen in screens:
        ## LFC ##
        each_screen_name = each_screen.split('.')[0]
        df_LFC_psig = pd.read_csv(os.path.join(workdir,f'screendata/{input_gene}_{each_screen_name}_proteinedits.tsv'),sep='\t')
        df_LFC3D_psig_tsv = f"{output_dir}/LFC3D/{input_gene}_NonAggr_LFC3D.tsv" # for LFC3D
        df_LFC3D_psig = pd.read_csv(df_LFC3D_psig_tsv, sep = "\t")
        df_LFC3D_dis_wght_tsv = f"{output_dir}/LFC3D/{input_gene}_LFC3D_dis_wght.tsv" # for LFC3D
        df_LFC3D_dis_wght_pd = pd.read_csv(df_LFC3D_dis_wght_tsv, sep = '\t')
        clustering_union(
            df_str_cons[["x_coord", "y_coord", "z_coord"]],
            df_LFC_psig, df_LFC3D_psig,
            workdir=output_dir, input_gene=input_gene, 
            screen_name= each_screen_name, 
            score_type='LFC',
            max_distances=6,
            pthr_cutoff=float(pthr_LFC))

        clustering_union(
            df_str_cons[["x_coord", "y_coord", "z_coord"]],
            df_LFC_psig, df_LFC3D_psig,
            workdir=output_dir, input_gene=input_gene, 
            screen_name= each_screen_name, 
            score_type='LFC3D',
            max_distances=6,
            pthr_cutoff=float(pthr_LFC))
        
        clustering_union(
            df_str_cons[["x_coord", "y_coord", "z_coord"]],
            df_LFC_psig, df_LFC3D_psig,
            workdir=output_dir, input_gene=input_gene, 
            screen_name= each_screen_name, 
            score_type='union',
            max_distances=6,
            pthr_cutoff=float(pthr_LFC))
        
        # Plotting LFC3D scatterplot
        scatterplot_by_residue_LFC3D(df_LFC3D_dis_wght_pd,output_dir,input_gene,each_screen_name)
        lfc_vs_lfc3d_scatterplot(df_LFC3D_dis_wght_pd, df_LFC3D_psig, output_dir, input_gene, each_screen_name, f'{workdir}/plots/{each_screen_name}_lfc_vs_lfc3d.png')
    # Meta-aggregation Clustering for LFC and LFC3D
    from BEClust3D.beclust3d.average_split_bin_metaaggregation import average_split_meta, bin_meta, znorm_meta
    from BEClust3D.beclust3d.average_split_bin_plot import average_split_bin_plots
    from BEClust3D.beclust3d.annotate_spatial_clusters import clustering_union

    lfc_lfc3d_df = pd.read_csv(f"{output_dir}/LFC3D/{input_gene}_LFC_LFC3D_LFC3Dr.tsv", sep='\t')
    meta_dict = {'LFC':pthr_LFC,'LFC3D':pthr_LFC3D}
    
    for score_type, each_pthr in meta_dict.items():
        df_bidir_meta = average_split_meta(
            df_LFC_LFC3D = lfc_lfc3d_df,
            workdir = output_dir,
            input_gene = input_gene,
            screen_names = screen_names,
            score_type=score_type,
            nRandom=nRandom,
            )
        _, df_neg_stats, df_pos_stats = bin_meta(
            df_bidir_meta,
            workdir = output_dir,
            input_gene = input_gene,
            score_type=score_type,        
            )
        
        df_meta_Z = znorm_meta(
            df_bidir_meta,
            df_neg_stats, 
            df_pos_stats,
            workdir = output_dir,
            input_gene = input_gene,
            score_type=score_type,        
            )

        average_split_bin_plots(
            df_meta_Z,
            workdir = output_dir, input_gene = input_gene, name='Meta',
            func='SUM', 
            pthr=each_pthr,
            score_type=score_type,
            )
    
   
    LFC_meta_tsv = f"{output_dir}/metaaggregation/{input_gene}_MetaAggr_LFC.tsv"
    LFC3D_meta_tsv = f"{output_dir}/metaaggregation/{input_gene}_MetaAggr_LFC3D.tsv"
        
    df_LFC_meta = pd.read_csv(LFC_meta_tsv, sep='\t')
    df_LFC3D_meta = pd.read_csv(LFC3D_meta_tsv, sep='\t')
       
    clustering_union(
        df_str_cons[["x_coord", "y_coord", "z_coord"]],
        df_LFC_meta, df_LFC3D_meta,
        workdir=output_dir, input_gene=input_gene, 
        screen_name= each_screen_name, 
        score_type='LFC',
        max_distances=6,
        pthr_cutoff=float(pthr_LFC),
        meta_aggregation=True)

    clustering_union(
        df_str_cons[["x_coord", "y_coord", "z_coord"]],
        df_LFC_meta, df_LFC3D_meta,
        workdir=output_dir, input_gene=input_gene, 
        screen_name= each_screen_name, 
        score_type='LFC3D',
        max_distances=6,
        pthr_cutoff=float(pthr_LFC),
        meta_aggregation=True)
    
    clustering_union(
        df_str_cons[["x_coord", "y_coord", "z_coord"]],
        df_LFC_meta, df_LFC3D_meta,
        workdir=output_dir, input_gene=input_gene, 
        screen_name= each_screen_name, 
        score_type='union',
        max_distances=6,
        pthr_cutoff=float(pthr_LFC),
        meta_aggregation=True)

if __name__ == '__main__':
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
    
    workdir = './MEN1' 
    screen_dir = '.' 
    screens = 'molm13.tsv,mv411.tsv' 
    input_gene = 'MEN1' 
    input_uniprot = 'O00255'
    input_pdb = 'men1_AF3.pdb'
    input_fasta = 'men1.fasta'
    output_dir = './MEN1'

    nRandom = 300
    pthr_LFC = 0.001
    pthr_LFC3D = 0.001

    mut_list_col = ""
    mut_col = "mutation_category"
    val_col = "delta_beta_score"
    gene_col = "Gene Symbol"
    edits_col = "predicted_edit"
    gRNA_col = None
    
    assign_hierachy = False
    hierachy = "Nonsense, Splice, Intron, Missense, Silent, UTR,  No Mutation" # original, "Nonsense, Splice Site, Missense, Intron, Silent, UTR, Flank, No Mutation"
    main(workdir=workdir, screen_dir=screen_dir, screens=screens,input_gene=input_gene, input_uniprot=input_uniprot, output_dir=output_dir, nRandom=nRandom,\
        pthr_LFC=pthr_LFC,pthr_LFC3D=pthr_LFC3D,mut_list_col=mut_list_col,input_pdb=input_pdb, input_fasta=input_fasta,\
        mut_col=mut_col, val_col=val_col, gene_col=gene_col, edits_col=edits_col, gRNA_col=gRNA_col, assign_hierachy=assign_hierachy, hierachy=hierachy)
