"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 2

"""

import pandas as pd
from pathlib import Path

def conservation(
        workdir, input_human_gene, input_mouse_gene, 
): 
    
    edits_filedir = Path(workdir + input_human_gene)
    # prepare file name
    input_Alignment_file = edits_filedir / f"Human{input_mouse_gene}_Mouse{input_mouse_gene}.align"
    iAlign = open(input_Alignment_file, "r")
    i_lines = iAlign.readlines()
    n_lines = len(i_lines)
    iAlign.close()

    ### this needs to be fixed cuz parsing the file is hardcoded
    start = 3 # depends on from which line number the alignment starts, for MUSCLE it's the 4th line
    n_score_line = 2 # in a block of sequence and score, the number of score line
    block_size = 4 # number of sequence + score line + space line
    j = 0
    blocks = []
    total_length = 0
    score_line_full = ''
    sequence_line_full_1 = ''
    sequence_line_full_2 = ''
    for i in range(start, n_lines):    
        seq_line = i_lines[i].strip('\n')
        l_line = len(seq_line)
        seq_line = seq_line[27:l_line] # IMPORTANT: WILL BE GENE DEPENDENT (Atf7ip:27, setdb1:28, H2-T23:26)
        blocks.append(seq_line)
        j = j + 1
        if j % block_size == 0:
            score_line_full = score_line_full + blocks[n_score_line]
            sequence_line_full_2 = sequence_line_full_2 + blocks[n_score_line - 1]; # first sequence above the score
            sequence_line_full_1 = sequence_line_full_1 + blocks[n_score_line - 2]; # second sequence above the score
            total_length = total_length + len(blocks[n_score_line])
            blocks = []

    # Collect info into a CSV
    score_list_filename = edits_filedir / f"Human{input_human_gene}_Mouse{input_mouse_gene}_align_conservation.tsv"
    score_list_file = open(score_list_filename, "w")
    score_list_file.write("alignment_pos\thuman_aligned_res\tmouse_aligned_res\tscore\tdis_score\tv_score\n")

    for i in range(0, len(sequence_line_full_1)):
        this_pos = i + 1
        human_aligned_res = sequence_line_full_1[i]
        mouse_aligned_res = sequence_line_full_2[i]
        this_score = score_line_full[i]
        
        cons_dict = {
            '*': ('conserved', 3),
            ':': ('similar', 2),
            '.': ('weakly_similar', 1),
            ' ': ('not_conserved', -1),
        }
        this_dis_score, this_v_score = cons_dict[this_score]
            
        write_line = str(this_pos) + "\t" + str(human_aligned_res) + "\t" + str(mouse_aligned_res)
        write_line = write_line + "\t" + str(this_score) + "\t" + str(this_dis_score) + "\t" + str(this_v_score)
        score_list_file.write(write_line  + "\n")
        score_list_file.close()

        df_alignment = pd.read_csv(score_list_filename, sep = '\t')

        residuemap_list_filename = edits_filedir / f"Human{input_human_gene}_Mouse{input_mouse_gene}_residuemap_conservation.tsv"
        residuemap_list_file = open(residuemap_list_filename, "w")
        residuemap_list_file.write("alignment_pos\thuman_res_pos\thuman_res\tmouse_res_pos\tmouse_res\tconservation\n")

        human_res_pos = 0
        mouse_res_pos = 0
        for i in range(0,len(df_alignment)):
            alignment_pos = df_alignment.at[i, 'alignment_pos']
            
            human_aligned_res = df_alignment.at[i, 'human_aligned_res']
            mouse_aligned_res = df_alignment.at[i, 'mouse_aligned_res']
            dis_score = df_alignment.at[i, 'dis_score']
            
            if human_aligned_res != '-' and mouse_aligned_res != '-':
                human_res_pos = human_res_pos + 1
                mouse_res_pos = mouse_res_pos + 1
                output_line = '\t'.join([str(alignment_pos), str(human_res_pos), human_aligned_res, 
                                         str(mouse_res_pos), mouse_aligned_res, dis_score])
                residuemap_list_file.write(output_line + '\n')
            elif human_aligned_res == '-' and mouse_aligned_res != '-':
                mouse_res_pos = mouse_res_pos + 1
            elif mouse_aligned_res == '-' and human_aligned_res != '-': # nothing to map from mouse 
                human_res_pos = human_res_pos + 1
                output_line = '\t'.join([str(alignment_pos), str(human_res_pos), human_aligned_res, 
                                         str(mouse_res_pos), mouse_aligned_res, dis_score])
                residuemap_list_file.write(output_line + '\n')
                
        residuemap_list_file.close()

        residuemap_file = workdir + input_human_gene + "/" + "Human" + input_human_gene + "_Mouse" + input_mouse_gene + "_residuemap_conservation.tsv"
        df_residuemap = pd.read_csv(residuemap_file, sep = '\t')

        return df_residuemap
