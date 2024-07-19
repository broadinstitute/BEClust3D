import pandas as pd
import filecmp
import pytest
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import warnings

from variables import *
from beclust3d.prioritize_by_conservation import prioritize_by_conservation

@pytest.mark.parametrize(("gene", "uniprot", "structid", "mouse_gene"), zip(all_genes, all_uniprots, all_structureids, all_mouse_genes))
@pytest.mark.parametrize("screen", all_human_screens)
def test_afstructuralfeatures_human(gene, uniprot, structid, mouse_gene, screen): 

    df_struc = pd.read_csv(f"{workdir}/{gene}/{structid}_coord_struc_features.tsv", sep = "\t")
    df_consrv = pd.read_csv(f"{workdir}/{gene}/Human{gene}_Mouse{mouse_gene}_residuemap_conservation.tsv", sep = '\t')

    prioritize_by_conservation(
        df_struc     =df_struc, 
        df_consrv    =df_consrv, 
        workdir      =workdir, 
        input_gene   =gene, 
        input_screen =screen, 
        structureid  =structid, 
    )

    screen_name = screen.split('.')[0]
    for f in os.listdir(f'{workdir}/{gene}/screendata'): 
        new_file = f'{workdir}/{gene}/screendata/{f}'
        ref_file = f'{refdir}/{gene}/screendata/{f}'
        if os.path.exists(new_file) and os.path.exists(ref_file): 
            if 'struc_consrv' in new_file and 'struc_consrv' in ref_file: 
                if screen_name in new_file and screen_name in ref_file: 
                    similarity_ratio = fuzzy_compare(new_file, ref_file)
                    if 0.9 > similarity_ratio > 0.8: 
                        warnings.warn(f"Similarity Ratio between {new_file} and {ref_file} are moderate")
                    assert similarity_ratio >= 0.8
                    # 0.9 is the arbitrary threshold of similarity

                    ### validate if this is a good threshold

# use a fuzzy compare function because there are some minor discrepancies with capitalization and ordering
# could also use difflib.SequenceMatcher, Levenshtein.ratio, ngrams word_tokenize (Jaccard similarity)
# but the TfidfVectorizer cosine_similarity is the fastest for large documents
def fuzzy_compare(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        text1 = f1.read()
        text2 = f2.read()

    vectorizer = TfidfVectorizer().fit_transform([text1, text2])
    vectors = vectorizer.toarray()
    similarity = cosine_similarity(vectors)
    score = similarity[0][1]
    return score
