"""
File: conservation.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 2

"""

import pandas as pd
from pathlib import Path
import requests
import sys, os
import time
import warnings

cons_dict = {
    '*': ('conserved', 3),
    ':': ('similar', 2),
    '.': ('weakly_similar', 1),
    ' ': ('not_conserved', -1),
}

def conservation(
        workdir, 
        input_human_gene, input_mouse_gene, 
        input_human_uniid, input_mouse_uniid, 
        email, title, 
        wait_time=30, 
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
    edits_filedir = Path(workdir + '/' + input_human_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)

    # QUERY UNIPROT AND WRITE TO sequences.fasta #
    request_filename_human, request_filename_mouse = f"{input_human_uniid}.fasta", f"{input_mouse_uniid}.fasta"
    human_seq = query_protein_fasta(edits_filedir, request_filename_human)
    mouse_seq = query_protein_fasta(edits_filedir, request_filename_mouse)
    seqs_filename = f"sequences.fasta"
    with open(edits_filedir / seqs_filename, "w") as text_file:
        text_file.write(human_seq)
        text_file.write(mouse_seq)

    # MUSCLE ALIGNMENT #
    align_filename = f"Human{input_human_gene}_Mouse{input_mouse_gene}.align"
    muscle_id = alignment_muscle(edits_filedir, seqs_filename, align_filename, 
                                 email, title, wait_time )
    
    # PARSE ALIGNMENT #
    alignconserv_filename = f"Human{input_human_gene}_Mouse{input_mouse_gene}_align_conservation.tsv"
    residuemap_filename =  f"Human{input_human_gene}_Mouse{input_mouse_gene}_residuemap_conservation.tsv"
    df_alignconserv, df_residuemap = parse_alignment(edits_filedir, align_filename, 
                                                     alignconserv_filename, residuemap_filename)

    return df_alignconserv, df_residuemap

def query_protein_fasta(
    edits_filedir, request_filename
): 
    """
    Description
        Query a Uniprot ID for .fasta file
    Params
        edits_filedir: str, required
            the Path to the main directory
    Returns
        response_body: str
            The amino acid sequence of the Uniprot ID
    """

    # QUERY UNIPROT #
    url = "https://rest.uniprot.org/uniprotkb/"
    requestURL = url+request_filename
    response = requests.get(requestURL)
    if not response.ok:
        response.raise_for_status()
        sys.exit()

    # RESPONSE TEXT #
    response_body = response.text
    with open(edits_filedir / request_filename, "w") as text_file:
        text_file.write(response_body)
    return response_body

def alignment_muscle(
    edits_filedir, 
    seqs_filename, align_filename, 
    email, title, wait_time, 
): 
    """
    Description
        Query ebi.ac.uk for MUSCLE alignment
    Params
        edits_filedir: str, required
            the Path to the main directory
        email: str, required
            email
        title: str, required
            title for your queried job
    Returns
        job_url_code: str
            a code for the MUSCLE job
    """
    # POST MUSCLE #
    url = 'https://www.ebi.ac.uk/Tools/services/rest/muscle/run'
    files = {'sequence': open(edits_filedir / seqs_filename, 'rb')}
    data = {
        'email': email, 'title': title, 'format': 'clw', 'tree': 'none',
    }
    response = requests.post(url, data=data, files=files)
    if response.status_code != 200: 
        warnings.warn('Error with MUSCLE query')
    time.sleep(wait_time)

    # GET MUSCLE #
    job_url_code = response.text.split(' ')[-1]
    print(f'Job ID: {job_url_code}')
    url = f'https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{job_url_code}/aln-clustalw'
    response = requests.get(url)

    with open(edits_filedir / align_filename, "wb") as f:
        f.write(response.content)
    return job_url_code

def parse_alignment(
    edits_filedir, 
    align_filename, alignconserv_filename, residuemap_filename, 
): 
    # PARSE .align FILE #
    iAlign = open(edits_filedir / align_filename, "r")
    i_lines = iAlign.readlines()[3:]
    iAlign.close()
    ind = i_lines[0].rfind(' ')

    # EXTRACT ALIGNMENT RESIDUES, ASSIGNS CONSERVATION SCORES, ADD TO DATAFRAME #
    # takes into account diff start positions and diff spacing
    cleaned_ilines = [s for s in [s[ind+1:].strip('\n') for s in i_lines] if len(s) > 0]
    human_align_res = list(''.join([s for i, s in enumerate(cleaned_ilines) if     i%3 == 0]))
    mouse_align_res = list(''.join([s for i, s in enumerate(cleaned_ilines) if (i-1)%3 == 0]))
    score           = list(''.join([s for i, s in enumerate(cleaned_ilines) if (i-2)%3 == 0]))

    index  = [i+1 for i in range(len(human_align_res))]
    dis, v = [cons_dict[s][0] for s in score], [cons_dict[s][1] for s in score]
    colnames = ['alignment_pos', 'human_aligned_res', 'mouse_aligned_res', 'score', 'dis_score', 'v_score']
    colvals  = [index, human_align_res, mouse_align_res, score, dis, v]

    df_alignconserv = pd.DataFrame()
    for name, col in zip(colnames, colvals): 
        df_alignconserv[name] = col
    df_alignconserv.to_csv(edits_filedir / alignconserv_filename, sep='\t', index=False)

    # COUNT CONSERVATION POSITIONS, ADD TO DATAFRAME #
    i, j = 0, 0
    human_res_pos, mouse_res_pos = [], []
    for s in human_align_res: 
        if s != '-': 
            i += 1
        human_res_pos.append(i)
    for s in mouse_align_res: 
        if s != '-': 
            j += 1
        mouse_res_pos.append(j)

    colnames = ['alignment_pos', 'human_res_pos', 'human_res', 'mouse_res_pos', 'mouse_res', 'conservation']
    colvals  = [index, human_res_pos, human_align_res, mouse_res_pos, mouse_align_res, dis]
    df_residuemap = pd.DataFrame()
    for name, col in zip(colnames, colvals): 
        df_residuemap[name] = col
    df_residuemap = df_residuemap[df_residuemap['human_res'] != '-']
    df_residuemap.to_csv(edits_filedir / residuemap_filename, sep='\t', index=False)

    del cleaned_ilines, index, human_align_res, mouse_align_res, score, dis, v, human_res_pos, mouse_res_pos, colnames, colvals

    return df_alignconserv, df_residuemap
