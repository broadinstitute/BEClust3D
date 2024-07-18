
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
