#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 06/22/2022
# Last Modified Date: 06/22/2022
# Last Modified By  : Dian Li 
import os
__version__="0.01"


#Precompute data, to infer biomarker, infiltration....
# PrepareData_R = os.path.join('baseP', 'prepare_data', 'data_prepare.R')
# PrepareData_Py = os.path.join('baseP', 'prepare_data', 'data_prepare.py')

#Immune related CRISPR SCREEN
R_RMD = os.path.join('baseP', 'evaluator', 'InvokeR.R')

#################### FIGURE SETTING ####################
Figure_Style='paper'

#################### CCLE lineage SETTING to be plotted ####################
CCLE_show_lineages = {
    "Exprsn": [
        "blood",
        "bone",
        "breast",
        "central_nervous_system",
        "colorectal",
        "fibroblast",
        "gastric",
        "liver",
        "lung",
        "lymphocyte",
        "peripheral_nervous_system",
        "skin",
        "query_cell_lines"
    ],
    "CRISPR_Broad": [
        "blood",
        "bone",
        "breast",
        "central_nervous_system",
        "colorectal",
        "gastric",
        "kidney",
        "liver",
        "lung",
        "lymphocyte",
        "peripheral_nervous_system",
        "skin",
        "query_cell_lines"
    ],
    "Proteomics": [
        "Bone",
        "Breast",
        "Central Nervous System",
        "Haematopoietic and Lymphoid Tissue",
        "Large Intestine",
        "Liver",
        "Lung",
        "Lymphoma",
        "Skin",
        "Stomach",
        "query_cell_lines"
    ]
}


#################### GTEx lineage SETTING to be plotted ####################
GTEx_show_lineages = {
    "Exprsn": [
        "Blood",
        "Brain",
        "Breast",
        "Colon",
        "Heart",
        "Liver",
        "Lung",
        "Nerve",
        "Skin",
        "Small Intestine",
        "Spleen"
    ]
}


#################### default cell line in case user doesn't provide a list ####################
default_cell_line = [
    "A549",
    "AF22",
    "HaCaT",
    "HEK 293",
    "HeLa",
    "Jurkat",
    "K-562",
    "Raji",
    "SH-SY5Y",
    "SNU738",
    "THP-1",
    "U937"
]
    
    
HPA_show_lineages = {
    "Exprsn": [
        "all_blood_cells",
        "all_cell_lines",
        "query_cell_lines"
    ]
}
    


