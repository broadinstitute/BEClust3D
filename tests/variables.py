
workdir = 'tests/results'
refdir = 'tests/reference'

all_genes = [
    'SETDB2', 'PPHLN1', 'MPHOSPH8', 'MORC2', 
    'ATF7IP', 'EHMT1', 'TASOR', 'SUV39H2', 
    'SETDB1', 'EHMT2', 'ATF7IP2', 
    'TRIM28', 'SUV39H1', 
]
all_uniprots = [
    'Q96T68', 'Q8NEY8', 'Q99549', 'Q9Y6X9', 
    'Q6VMQ6', 'Q9H9B1', 'Q9UK61', 'Q9H5I1', 
    'Q15047', 'A2ABF9', 'Q5U623', 
    'Q13263', 'O43463', 
]
all_structureids = [f"AF-{uni}-F1-model_v4" for uni in all_uniprots]

all_human_screens = [
    'Human_InVitro_294T_Apobec_D14_Input.txt',
    'Human_InVitro_294T_TadA_D14_Input.txt',
    'Human_InVitro_SK-ES_Apobec_D14_Input.txt',
    'Human_InVitro_SK-ES_TadA_D17_Input.txt',
    'Human_InVitro_SW480_Apobec_D14_Input.txt',
    'Human_InVitro_SW480_TadA_D14_Input.txt',
]
all_mouse_screens = [
    'Mouse_InVivo_B16_Apobec_NSG_Input.txt', 
    'Mouse_InVivo_B16_Apobec_Tx_NSG.txt', 
    'Mouse_InVivo_B16_TadA_NSG_Input.txt', 
    'Mouse_InVivo_B16_TadA_Tx_NSG.txt', 
    'Mouse_InVivo_LLC_Apobec_NSG_Input.txt', 
    'Mouse_InVivo_LLC_Apobec_Tx_NSG.txt', 
]

all_mouse_genes = [
    'Setdb2', 'Pphln1', 'Mphosph8', 'Morc2a', 
    'Atf7ip', 'Ehmt1', 'Fam208a', 'Suv39h2', 
    'Setdb1', 'Ehmt2', 'Atf7ip2', 'Ptpn2', 
    'Trim28', 'Suv39h1'
]
all_mouse_uniprots = [
    'Q8C267', 'Q8K2H1', 'Q3TYA6', 'Q69ZX6', 
    'Q7TT18', 'Q5DW34', 'Q69ZR9', 'Q9EQQ0', 
    'D3YYC3', 'Q9Z148', 'Q3UL97', 'Q06180', 
    'Q62318', 'O54864', 
]

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

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
