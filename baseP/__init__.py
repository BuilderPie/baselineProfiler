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
        "gastric",
        "liver",
        "lung",
        "lymphocyte",
        "peripheral_nervous_system",
        "skin"
    ],
    "CRISPR_Broad": [
        "blood",
        "bone",
        "breast",
        "central_nervous_system",
        "colorectal",
        "gastric",
        "liver",
        "lung",
        "lymphocyte",
        "peripheral_nervous_system",
        "skin"
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
       "Stomach"
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




    
    

    


