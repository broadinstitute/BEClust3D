"""
File: uniprot_annot.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 6

"""

import pandas as pd
from pathlib import Path

def uniprot_annot(df_struc_consvr, 
        workdir, input_human_gene, input_human_uniprot, 
): 
    
    df_META = pd.DataFrame()
    df_META["unipos"] = df_struc_consvr["unipos"]
    df_META["unires"] = df_struc_consvr["unires"]

    edits_filedir = Path(workdir + input_human_gene)

    uniprot_source_filename = edits_filedir / f"{input_human_uniprot}.txt"
    uniprot_source_file = open(uniprot_source_filename, "r")
    uniprot_source_lines = uniprot_source_file.readlines()
    uniprot_source_file.close()

    uniprot_processed_lines = [idx for idx in uniprot_source_lines if idx[0:2] == "FT"]

    uniprot_processed_filename = edits_filedir / f"{input_human_gene}_{input_human_uniprot}_processed.uniprot";
    uniprot_processed_file = open(uniprot_processed_filename, "w")
    uniprot_processed_file.writelines(uniprot_processed_lines)
    uniprot_processed_file.close()

    arr_act_site, arr_binding, arr_site = [], [], []
    arr_dna_bind, arr_zn_fing = [], []
    arr_region, arr_motif, arr_coiled, arr_compbias, arr_repeat = [], [], [], [], []
    arr_domain, arr_topo_dom, arr_transmem, arr_intramem = [], [], [], []
    arr_mod_res, arr_lipid, arr_carbohyb, arr_disulfide, arr_crosslnk = [], [], [], [], []
    arr_signal, arr_transit, arr_propep, arr_peptide = [], [], [], []

    codes = ['ACT_SITE', 'BINDING', 'SITE', 
             'DNA_BIND', 'ZN_FING', 
             'REGION', 'MOTIF', 'COILED', 'COMPBIAS', 'REPEAT', 
             'DOMAIN', 'TOPO_DOM', 'TRANSMEM', 'INTRAMEM', 
             'MOD_RES', 'LIPID', 'CARBOHYD', 'DISULFID', 'CROSSLNK', 
             'SIGNAL', 'TRANSIT', 'PROPEP', 'PEPTIDE', ]
    arr_list = [arr_act_site, arr_binding, arr_site, 
                arr_dna_bind, arr_zn_fing,
                arr_region, arr_motif, arr_coiled, arr_compbias, arr_repeat,
                arr_domain, arr_topo_dom, arr_transmem, arr_intramem, 
                arr_mod_res, arr_lipid, arr_carbohyb, arr_disulfide, arr_crosslnk,
                arr_signal, arr_transit, arr_propep, arr_peptide,
                ]

    ### refactor as a dictionary instead

    for i in range(0, len(df_META)):
        unipos = df_META.at[i, 'unipos']
        
        # Sequence annotation: query fields: https://www.uniprot.org/help/sequence_annotation
                
        for arr, code in zip(arr_list, codes): 
            temp = '-'
            for j in range(0, len(uniprot_processed_lines)):
                featureline = uniprot_processed_lines[j].strip()
                if code in featureline:
                    positionline = featureline.strip()
                    positioninfo = positionline[21:len(positionline)]
                    positiontokens = positioninfo.split('..')
                    notesline = uniprot_processed_lines[j+1].strip()
                    notesinfo = notesline[21:len(notesline)]
                    notestoken = notesinfo.split('=')
                    notes = notestoken[1].strip()
                    notes = notes[1:len(notes)-1]
                    
                    if len(positiontokens) == 1:
                        position = int(positiontokens[0])
                        if unipos == position:
                            temp = str(positioninfo) + ":" + str(notes)
                    elif len(positiontokens) == 2:
                        start_position = int(positiontokens[0])
                        end_position = int(positiontokens[1])
                        if (unipos >= start_position) and (unipos <= end_position):
                            temp = str(positioninfo) + ":" + str(notes)
                    elif len(positiontokens) == 0:
                        print("feature position can't be empty!")
                        break
                    elif len(positiontokens) > 2:
                        print("feature positions can't be more than two!")
                        break
            arr.append(temp)
    
    df_META['active_site'] = arr_act_site
    df_META['binding_site'] = arr_binding
    df_META['site'] = arr_site
    df_META['dna_binding_region'] = arr_dna_bind
    df_META['zinc_finger'] = arr_zn_fing
    df_META['region_of_interest'] = arr_region
    df_META['motif'] = arr_motif
    df_META['coiled_coil'] = arr_coiled
    df_META['compositional_bias'] = arr_compbias
    df_META['repeat'] = arr_repeat
    df_META['domain'] = arr_domain
    df_META['topological_domain'] = arr_topo_dom
    df_META['transmembrane'] = arr_transmem
    df_META['intramembrane'] = arr_intramem
    df_META['modified_residue'] = arr_mod_res
    df_META['lipidation'] = arr_lipid
    df_META['glycosylation'] = arr_carbohyb
    df_META['disulfide_bind'] = arr_disulfide
    df_META['cross_links'] = arr_crosslnk
    df_META['signal_peptide'] = arr_signal
    df_META['trasit_peptide'] = arr_transit
    df_META['propeptide'] = arr_propep
    df_META['peptide'] = arr_peptide

    uniprot_feature_filename = edits_filedir / f"{input_human_uniprot}_{input_human_uniprot}_Uniprotf.txt"
    df_META.to_csv(uniprot_feature_filename, sep = '\t', index=False)

    return df_META
