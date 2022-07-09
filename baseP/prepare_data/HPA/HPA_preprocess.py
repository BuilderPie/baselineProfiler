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


def stat_HPA_Cell_Line():
    ######################### Part 1. preprocess HPA meta file ###############################
    HPA_Cell_Line_sample_info_file = os.path.join(DATA_COMPUTE, 'HPA', 'sample_info', 'Exprsn_cell_line_meta_complete.csv')

    HPA_Cell_Line_sample_info = pd.read_csv(HPA_Cell_Line_sample_info_file)
    
    os.makedirs(os.path.join(DATA_COMPUTE, 'HPA', 'sample_info', 'tables'), exist_ok=True)
    output_file = os.path.join(DATA_COMPUTE, 'HPA', 'sample_info', 'tables', 'HPA_Cell_Line_sample_info.txt')
    HPA_Cell_Line_sample_info.to_csv(output_file, index=False, sep='\t')
    ##### ======= convert CCLE_sample_info_subset to datatables format ============= #####
    os.makedirs(os.path.join(DATA_COMPUTE, 'HPA', 'sample_info', 'datatables'), exist_ok=True)
    txt_2_datatables(fn = 'HPA_Cell_Line_sample_info.txt',input_path = os.path.join(DATA_COMPUTE, 'HPA', 'sample_info', 'tables'),
    output_path = os.path.join(DATA_COMPUTE, 'HPA', 'sample_info', 'datatables'),
    columnsSelect = ['HPA_cell_line', 'lineage', 'Description'])

    
    ##### ======= generate a table showing the distribution of lineage and primary disease
    



def analyze():
    stat_HPA_Cell_Line()