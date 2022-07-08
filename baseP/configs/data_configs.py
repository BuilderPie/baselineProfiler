#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 06/21/2022
# Last Modified Date: 06/21/2022
# Last Modified By  : Dian Li 
import os
import baseP.configs.dirpath
STATIC_DIR = baseP.configs.dirpath.static
DATA_DIR = os.path.join(STATIC_DIR,'data/')
DATA_COMPUTE = os.path.join(STATIC_DIR,'data_computed/')

#pathway annotations
RcmdPath=os.path.join('Invoke_Rscripts.R')
#################### ANNOTATION SETTING ####################
mm10_ANNO_GENE = os.path.join(DATA_DIR, 'Annotation','mm_ensembl.txt')
hg38_ANNO_GENE = os.path.join(DATA_DIR, 'Annotation','hg_ensembl.txt')

HOM_MouseHumanSequence = os.path.join(DATA_DIR, 'Annotation', 'symbols', 'HOM_MouseHumanSequence.rpt')
HOM_MouseHumanSequence_processed = os.path.join(DATA_COMPUTE, 'Annotation', 'HOM_MouseHumanSequence_processed.csv')
#################### CCLE SETTING####################
CCLE_DIR = os.path.join(DATA_DIR, 'CCLE')
CCLE_COMPUTE = os.path.join(DATA_COMPUTE, 'CCLE')
CCLE_Data = {
    'Exprsn':  os.path.join(CCLE_COMPUTE, 'expression', 'expression.csv'),
    'Exprsn_rownames':  os.path.join(CCLE_COMPUTE, 'expression', 'expression_rownames.csv'),
    'Exprsn_colnames':  os.path.join(CCLE_COMPUTE, 'expression', 'expression_colnames.csv'),
    'Exprsn_meta':  os.path.join(CCLE_COMPUTE, 'sample_info', 'sample_info.csv'),

    'CRISPR_Broad':os.path.join(CCLE_COMPUTE,'CRISPR', 'Broad', 'gene_effect.csv'),
    'CRISPR_Broad_rownames':os.path.join(CCLE_COMPUTE,'CRISPR', 'Broad', 'gene_effect_rownames.csv'),
    'CRISPR_Broad_colnames':os.path.join(CCLE_COMPUTE,'CRISPR', 'Broad', 'gene_effect_colnames.csv'),
    'CRISPR_Broad_meta':os.path.join(CCLE_COMPUTE, 'sample_info', 'sample_info.csv'),

    'Proteomics':os.path.join(CCLE_COMPUTE,'proteomics','protein_normalized.csv'),
    'Proteomics_rownames':os.path.join(CCLE_COMPUTE,'proteomics', 'protein_normalized_rownames.csv'),
    'Proteomics_colnames':os.path.join(CCLE_COMPUTE,'proteomics', 'protein_normalized_colnames.csv'),
    'Proteomics_meta':os.path.join(CCLE_COMPUTE, 'sample_info', 'sample_info_proteomics.csv'),

    
    # 'copy_number': os.path.join(CCLE_COMPUTE, 'copy_number'),
}
# if not os.path.isdir(CCLE_Data['Exprsn']):
#         os.makedirs(CCLE_Data['Exprsn'])   # create CCLE Exprsn if it doesn't exists at the beginning
# CCLE_Data['Datasets'] = [x for x in os.listdir(CCLE_Data['Exprsn']) if not x.startswith('.')]

#################### GTEx SETTING####################
GTEx_COMPUTE = os.path.join(DATA_COMPUTE, 'GTEx')
GTEx_Data = {
    'Exprsn':  os.path.join(GTEx_COMPUTE, 'data_subset'),
    'Exprsn_rownames':  os.path.join(GTEx_COMPUTE, 'sample_info', 'gene_name.csv'),
}

#################### HPA SETTING####################
HPA_COMPUTE = os.path.join(DATA_COMPUTE, 'HPA')
HPA_Data = {
    'Exprsn_cell_line':  os.path.join(HPA_COMPUTE, 'Exprsn_cell_line', "rna_celline_log2_nTPM.csv"),
    'Exprsn_cell_line_rownames':  os.path.join(HPA_COMPUTE, 'Exprsn_cell_line', 'rna_celline_rownames.csv'),
    'Exprsn_cell_line_colnames':  os.path.join(HPA_COMPUTE, 'Exprsn_cell_line', 'rna_celline_colnames.csv'),

    'Exprsn_blood_cell':  os.path.join(HPA_COMPUTE, 'Exprsn_blood_cell', "rna_blood_log2_nTPM.csv"),
    'Exprsn_blood_cell_rownames':  os.path.join(HPA_COMPUTE, 'Exprsn_blood_cell', 'rna_blood_rownames.csv'),
    'Exprsn_blood_cell_colnames':  os.path.join(HPA_COMPUTE, 'Exprsn_blood_cell', 'rna_blood_colnames.csv'),
    
    'Exprsn_cell_line_meta':  os.path.join(HPA_COMPUTE, 'sample_info', "Exprsn_cell_line_meta.csv"),
}

#################### others SETTING####################
others_COMPUTE = os.path.join(DATA_COMPUTE, 'HPA')
others_Data = {
    'Exprsn_Neuro2a': os.path.join(others_COMPUTE, '')
}



#################### post HTML SETTING####################
HTML_DIR = os.path.join(DATA_DIR, 'Annotation', 'html_template')
HTML_Data = {
    'CSS': os.path.join(HTML_DIR, 'style.css'),
    'site_libs': os.path.join(HTML_DIR, 'site_libs'),
}

EMAIL_DIR = os.path.join(DATA_DIR, 'Annotation', 'email_config')
EMAIL_Data = {
    'send_email_sh': os.path.join(EMAIL_DIR, 'send_email.sh'),
}
