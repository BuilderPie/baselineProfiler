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

from baseP.configs.data_configs import DATA_DIR, DATA_COMPUTE, HOM_MouseHumanSequence_processed
from baseP.prepare_data.prepare_utils import txt_2_json, txt_2_datatables

# def geneset(request,logger_P, name, threads, evaluators, resum, html_modules, weights):

def analyze():
    input_file = os.path.join(DATA_DIR, 'others', 'Neuro2a', 'raw_log2_TPM.csv')
    df = pd.read_csv(input_file)
    anno = pd.read_csv(HOM_MouseHumanSequence_processed)

    df_merge = pd.merge(
        left=anno.loc[:, ["mouse_gene_symbol", "human_gene_symbol"]],
        right=df.iloc[:,~df.columns.isin(['gene_id'])], 
        how='left',
        left_on='mouse_gene_symbol',right_on='gene_name')

    df_merge = df_merge.iloc[:, ~df_merge.columns.isin(['gene_name'])]

    os.makedirs(os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a'), exist_ok=True)
    output_file = os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'raw_log2_TPM_name_matched.csv')
    df_merge.to_csv(output_file, index = False)

    output_file = os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'raw_log2_TPM_gene_names.csv')
    df_merge.loc[:, ['mouse_gene_symbol', 'human_gene_symbol']].to_csv(output_file, index = False)

    output_file = os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'raw_log2_TPM_column_names.csv')
    pd.DataFrame(data=df_merge.columns.tolist(), columns=['colnames']).to_csv(output_file, index = False)

    # ============= extract meta table subset ================ #
    sample_info_file = os.path.join(DATA_DIR, 'others', 'Neuro2a', 'meta_sheet_extended.csv')
    Neuro2a_sample_info = pd.read_csv(sample_info_file)

    os.makedirs(os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'tables'), exist_ok=True)
    output_file = os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'tables', 'Neuro2a_sample_info.txt')
    Neuro2a_sample_info.to_csv(output_file, index=False, sep='\t')
    # ============= convert to datatables format ============= #
    os.makedirs(os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'datatables'), exist_ok=True)
    txt_2_datatables(fn = 'Neuro2a_sample_info.txt',input_path = os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'tables'),
    output_path = os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'datatables'),
    columnsSelect = ['Study', 'Name_part_1', 'Name_part_2', 'PMID', 'Condition', 'Journal', 'Date'])


    