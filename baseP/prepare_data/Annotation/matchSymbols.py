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

from baseP.configs.data_configs import DATA_DIR, DATA_COMPUTE, HOM_MouseHumanSequence

# def geneset(request,logger_P, name, threads, evaluators, resum, html_modules, weights):

def analyze():
    os.makedirs(os.path.join(DATA_COMPUTE, 'Annotation'), exist_ok=True)
    output_file = os.path.join(DATA_COMPUTE, 'Annotation', 'HOM_MouseHumanSequence_processed.csv')
    out_list = []
    flag = 0
    with open(HOM_MouseHumanSequence, 'r') as MHS:
        for line in MHS:
            line = line.split('\t')
            if line[1] == 'mouse, laboratory':
                flag = 1
                tmp_dict = {'mouse_gene_symbol': line[3], 'mouse_EntrezGene_ID': line[4]}
                continue
            if flag == 1:
                tmp_dict['human_gene_symbol'] = line[3]
                tmp_dict['human_EntrezGene_ID'] = line[4]
                flag = 0
                out_list.append(tmp_dict)
    df_out = pd.DataFrame(out_list)
    df_out.to_csv(output_file, index = False)
    