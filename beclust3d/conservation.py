"""
File: conservation.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 2

"""

import pandas as pd
from pathlib import Path
import requests
import sys
# from Bio.Align.Applications import MuscleCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

def conservation(
        workdir, 
        input_human_gene, input_mouse_gene, 
        input_human_uniid, input_mouse_uniid, 
): 
    """
    Description
        Generate dataframes of sequence conservation for each residue. 

    Params
        workdir: str, required
            the working directory
        input_human_gene: str, required
            the name of the input human gene
        input_mouse_gene: str, required
            the name of the input mouse gene
        input_human_uniid: str, required
            the name of the input human uniprot id
        input_mouse_uniid: str, required
            the name of the input mouse uniprot id

    Returns
        df_alignconserv: pandas dataframe
            a dataframe listing sequence conservation and their categories and scores
            with column headers alignment_pos human_aligned_res mouse_aligned_res score dis_score v_score
        df_residuemap: pandas dataframe
            a dataframe listing sequence conservation and their indices for human and mouse
            with column headers alignment_pos human_res_pos human_res mouse_res_pos mouse_res conservation
    """

    # query .fasta protein sequences
    edits_filedir = Path(workdir + '/' + input_human_gene)
    human_seq = query_protein_fasta(input_human_uniid, edits_filedir)
    mouse_seq = query_protein_fasta(input_mouse_uniid, edits_filedir)
    with open(edits_filedir / f"sequences.fasta", "w") as text_file:
        text_file.write(human_seq)
        text_file.write(mouse_seq)

    # alignment
    alignment_file = edits_filedir / f"Human{input_human_gene}_Mouse{input_mouse_gene}.align"
    alignment_muscle(edits_filedir)
    
    # use MUSCLE to align human and mouse sequences
    iAlign = open(alignment_file, "r")
    i_lines = iAlign.readlines()[3:]
    iAlign.close()

    # parse .align file #
    ind = i_lines[0].rfind(' ')
    cleaned_ilines = [s[ind:].strip().strip('\n') for s in i_lines]
    cleaned_ilines = [s for s in cleaned_ilines if len(s) > 0]
    human_aligned_res = list(''.join([s for i, s in enumerate(cleaned_ilines) if i%3 == 0]))
    mouse_aligned_res = list(''.join([s for i, s in enumerate(cleaned_ilines) if (i-1)%3 == 0]))
    score             = list(''.join([s for i, s in enumerate(cleaned_ilines) if (i-2)%3 == 0]))
    # construct align_conservation dataframe #
    cons_dict = {
        '*': ('conserved', 3),
        ':': ('similar', 2),
        '.': ('weakly_similar', 1),
        ' ': ('not_conserved', -1),
    }
    index  = [i+1 for i in range(len(human_aligned_res))]
    dis, v = [cons_dict[s][0] for s in score], [cons_dict[s][1] for s in score]
    colnames = ['alignment_pos', 'human_aligned_res', 'mouse_aligned_res', 'score', 'dis_score', 'v_score']
    colvals  = [index, human_aligned_res, mouse_aligned_res, score, dis, v]
    df_alignconserv = pd.DataFrame()
    for name, col in zip(colnames, colvals): 
        df_alignconserv[name] = col
    # save align_conservation #
    alignconserv_filename = edits_filedir / f"Human{input_human_gene}_Mouse{input_mouse_gene}_align_conservation.tsv"
    df_alignconserv.to_csv(alignconserv_filename, sep='\t', index=False)

    # residue map conservation #
    i, j = 0, 0
    human_res_pos, mouse_res_pos = [], []
    for s in human_aligned_res: 
        if s != '-': i += 1
        human_res_pos.append(i)
    for s in mouse_aligned_res: 
        if s != '-': j += 1
        mouse_res_pos.append(j)
    # construct residuemap_conservation dataframe #
    df_residuemap = pd.DataFrame()
    colnames = ['alignment_pos', 'human_res_pos', 'human_res', 'mouse_res_pos', 'mouse_res', 'conservation']
    colvals  = [index, human_res_pos, human_aligned_res, mouse_res_pos, mouse_aligned_res, dis]
    for name, col in zip(colnames, colvals): 
        df_residuemap[name] = col
    df_residuemap = df_residuemap[df_residuemap['human_res'] != '-']
    # save residuemap_conservation #
    residuemap_list_filename = edits_filedir / f"Human{input_human_gene}_Mouse{input_mouse_gene}_residuemap_conservation.tsv"
    df_residuemap.to_csv(residuemap_list_filename, sep='\t', index=False)

    return df_alignconserv, df_residuemap

def query_protein_fasta(
        input_uniprot, edits_filedir
): 
    """
    Description
        Query a Uniprot ID for .fasta file

    Params
        input_uniprot: str, required
            The Uniprot ID for a particular protein

    Returns
        response_str: str
            The amino acid sequence corresponding to the Uniprot ID
    """

    url = "https://rest.uniprot.org/uniprotkb/"
    request = f"{input_uniprot}.fasta"
    requestURL = url+request
    response = requests.get(requestURL)

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    response_body = response.text
    response_list = response_body.split('\n')
    response_str = ''.join(response_list[1:])

    with open(edits_filedir / request, "w") as text_file:
        text_file.write(response_body)
    return response_body

# from Bio.Align.Applications import MuscleCommandline
# from Bio import AlignIO

def alignment_muscle(
        edits_filedir,
): 
    # unaligned_filepath = edits_filedir / "sequences.fasta"
    # out_filepath = edits_filedir / "aligned.fasta"
    # muscle_exe = edits_filedir / "muscle3.8.31_i86darwin64"
    # muscle_cline = MuscleCommandline(
    #     muscle_exe,
    #     input=unaligned_filepath, 
    #     out=out_filepath
    #     )
    # stdout, stderr = muscle_cline()
    # MultipleSeqAlignment = AlignIO.read(out_filepath, "fasta") 
    # print(MultipleSeqAlignment)
    return None

