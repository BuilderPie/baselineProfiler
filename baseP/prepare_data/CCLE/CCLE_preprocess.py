#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 07/07/2022
# Last Modified Date: 07/07/2022
# Last Modified By  : Dian Li
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from baseP.configs.data_configs import DATA_DIR, DATA_COMPUTE, HOM_MouseHumanSequence
from baseP.prepare_data.prepare_utils import txt_2_json, txt_2_datatables
from baseP.configs.data_configs import CCLE_Data

# def geneset(request,logger_P, name, threads, evaluators, resum, html_modules, weights):


def stat_CCLE_Exprsn():
    ######################### Part 1. preprocess CCLE expression file ###############################
    CCLE_expression_file = os.path.join(DATA_DIR, 'CCLE', 'CCLE_expression.csv')
    CCLE_sample_info_file = os.path.join(DATA_DIR, 'CCLE', 'sample_info.csv')

    CCLE_expression = pd.read_csv(CCLE_expression_file)
    CCLE_cell_line = pd.DataFrame(data=CCLE_expression.iloc[:,0].tolist(), columns=['DepMap_ID'])

    CCLE_sample_info = pd.read_csv(CCLE_sample_info_file)

    CCLE_sample_info_subset = pd.merge(
        left=CCLE_cell_line, 
        right=CCLE_sample_info.loc[:,['DepMap_ID', 'cell_line_name', 'stripped_cell_line_name', 'lineage', 'primary_disease']],
        how='left',
        on="DepMap_ID")

    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables'), exist_ok=True)

    output_file = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'CCLE_Expression_sample_info.txt')
    CCLE_sample_info_subset.to_csv(output_file, index = False, sep='\t')
    CCLE_sample_info_subset = pd.read_csv(output_file, sep='\t')
    
    ##### ======= convert CCLE_sample_info_subset to datatables format ============= #####
    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables'), exist_ok=True)
    txt_2_datatables(fn = 'CCLE_Expression_sample_info.txt',input_path = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables'),
    output_path = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables'),
    columnsSelect = ['DepMap_ID', 'cell_line_name', 'stripped_cell_line_name', 'lineage', 'primary_disease'])

    ##### ======= generate a table showing the distribution of lineage and primary disease
    disease_count = CCLE_sample_info_subset['primary_disease'].value_counts().to_frame()
    disease_count.columns = ['value']
    disease_count['individual'] = disease_count.index
    disease_count['group'] = 'Disease'

    lineage_count = CCLE_sample_info_subset['lineage'].value_counts().to_frame()
    lineage_count.columns = ['value']
    lineage_count['individual'] = lineage_count.index
    lineage_count['group'] = 'Lineage'

    df_out = pd.concat([disease_count, lineage_count], axis=0)
    output_file = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'stats_CCLE_Expression_sample_info.txt')
    df_out.to_csv(output_file, sep = '\t', index=False)

    ###### ======= generate a circular barplot from the count table ======== ######
    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'figs'), exist_ok=True)

    R_RMD = os.path.join('baseP', 'prepare_data', 'R_functions', 'src', 'plot_barplot.R')
    
    r_cmd = ' '.join(['Rscript', R_RMD, 
    '--file ', os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'stats_CCLE_Expression_sample_info.txt').replace(' ', '\ '), 
    '--output ', os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'figs', 'stats_CCLE_Expression_sample_info.png').replace(' ', '\ '),
    ])
    process = subprocess.Popen(r_cmd, shell=True).wait()

def stat_CCLE_CRISPR_Broad():
    ######################### Part 1. preprocess CCLE expression file ###############################
    CCLE_expression_file = os.path.join(DATA_DIR, 'CCLE', 'CRISPR_gene_effect.csv')
    CCLE_sample_info_file = os.path.join(DATA_DIR, 'CCLE', 'sample_info.csv')

    CCLE_expression = pd.read_csv(CCLE_expression_file)
    CCLE_cell_line = pd.DataFrame(data=CCLE_expression.iloc[:,0].tolist(), columns=['DepMap_ID'])

    CCLE_sample_info = pd.read_csv(CCLE_sample_info_file)

    CCLE_sample_info_subset = pd.merge(
        left=CCLE_cell_line, 
        right=CCLE_sample_info.loc[:,['DepMap_ID', 'cell_line_name', 'stripped_cell_line_name', 'lineage', 'primary_disease']],
        how='left',
        on="DepMap_ID")

    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables'), exist_ok=True)

    output_file = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'CCLE_CRISPR_Broad_sample_info.txt')
    CCLE_sample_info_subset.to_csv(output_file, index = False, sep='\t')
    CCLE_sample_info_subset = pd.read_csv(output_file, sep='\t')
    
    ##### ======= convert CCLE_sample_info_subset to datatables format ============= #####
    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables'), exist_ok=True)
    txt_2_datatables(fn = 'CCLE_CRISPR_Broad_sample_info.txt',input_path = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables'),
    output_path = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables'),
    columnsSelect = ['DepMap_ID', 'cell_line_name', 'stripped_cell_line_name', 'lineage', 'primary_disease'])

    ##### ======= generate a table showing the distribution of lineage and primary disease
    disease_count = CCLE_sample_info_subset['primary_disease'].value_counts().to_frame()
    disease_count.columns = ['value']
    disease_count['individual'] = disease_count.index
    disease_count['group'] = 'Disease'

    lineage_count = CCLE_sample_info_subset['lineage'].value_counts().to_frame()
    lineage_count.columns = ['value']
    lineage_count['individual'] = lineage_count.index
    lineage_count['group'] = 'Lineage'

    df_out = pd.concat([disease_count, lineage_count], axis=0)
    output_file = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'stats_CCLE_CRISPR_Broad_sample_info.txt')
    df_out.to_csv(output_file, sep = '\t', index=False)

    ###### ======= generate a circular barplot from the count table ======== ######
    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'figs'), exist_ok=True)

    R_RMD = os.path.join('baseP', 'prepare_data', 'R_functions', 'src', 'plot_barplot.R')
    
    r_cmd = ' '.join(['Rscript', R_RMD, 
    '--file ', os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'stats_CCLE_CRISPR_Broad_sample_info.txt').replace(' ', '\ '), 
    '--output ', os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'figs', 'stats_CCLE_CRISPR_Broad_sample_info.png').replace(' ', '\ '),
    ])
    process = subprocess.Popen(r_cmd, shell=True).wait()

def stat_CCLE_Proteomics():
    ######################### Part 1. preprocess CCLE expression file ###############################
    CCLE_sample_info_file = os.path.join(DATA_DIR, 'CCLE', 'sample_info.csv')
    CCLE_sample_info = pd.read_csv(CCLE_sample_info_file)

    CCLE_proteomics_meta_file = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'sample_info_proteomics.csv')
    CCLE_proteomics_meta = pd.read_csv(CCLE_proteomics_meta_file)
    CCLE_proteomics_meta = CCLE_proteomics_meta.loc[:, ['cell_line', 'cell_line_CCLE', 'lineage']]

    CCLE_sample_info_subset = pd.merge(
        left = CCLE_proteomics_meta, 
        right = CCLE_sample_info.loc[:,['DepMap_ID', 'stripped_cell_line_name', 'primary_disease']],
        how = 'left',
        left_on='cell_line_CCLE',
        right_on='stripped_cell_line_name')

    CCLE_sample_info_subset = CCLE_sample_info_subset.loc[:, ['DepMap_ID', 'cell_line', 'cell_line_CCLE', 'lineage', 'primary_disease']]

    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables'), exist_ok=True)
    output_file = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'CCLE_Proteomics_sample_info.txt')
    CCLE_sample_info_subset.to_csv(output_file, index = False, sep='\t')
    CCLE_sample_info_subset = pd.read_csv(output_file, sep='\t')
    ##### ======= convert CCLE_sample_info_subset to datatables format ============= #####
    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables'), exist_ok=True)
    txt_2_datatables(fn = 'CCLE_Proteomics_sample_info.txt',input_path = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables'),
    output_path = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables'),
    columnsSelect = ['DepMap_ID', 'cell_line', 'cell_line_CCLE', 'lineage', 'primary_disease'])

    ##### ======= generate a table showing the distribution of lineage and primary disease
    disease_count = CCLE_sample_info_subset['primary_disease'].value_counts().to_frame()
    disease_count.columns = ['value']
    disease_count['individual'] = disease_count.index
    disease_count['group'] = 'Disease'

    lineage_count = CCLE_sample_info_subset['lineage'].value_counts().to_frame()
    lineage_count.columns = ['value']
    lineage_count['individual'] = lineage_count.index
    lineage_count['group'] = 'Lineage'

    df_out = pd.concat([disease_count, lineage_count], axis=0)
    output_file = os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'stats_CCLE_Proteomics_sample_info.txt')
    df_out.to_csv(output_file, sep = '\t', index=False)

    ###### ======= generate a circular barplot from the count table ======== ######
    os.makedirs(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'figs'), exist_ok=True)

    R_RMD = os.path.join('baseP', 'prepare_data', 'R_functions', 'src', 'plot_barplot.R')
    
    r_cmd = ' '.join(['Rscript', R_RMD, 
    '--file ', os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'tables', 'stats_CCLE_Proteomics_sample_info.txt').replace(' ', '\ '), 
    '--output ', os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'figs', 'stats_CCLE_Proteomics_sample_info.png').replace(' ', '\ '),
    ])
    process = subprocess.Popen(r_cmd, shell=True).wait()

def analyze():
    stat_CCLE_Exprsn()
    stat_CCLE_CRISPR_Broad()
    stat_CCLE_Proteomics()