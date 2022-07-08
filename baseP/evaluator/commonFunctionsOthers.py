#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 07/07/2022
# Last Modified Date: 07/07/2022
# Last Modified By  : Dian Li
# -*- coding:utf-8 -*-

import os
from os import listdir, stat
import warnings
warnings.simplefilter('ignore')
# try:
#    import cPickle as pickle
# except:
#    import pickle

import numpy as np
import pandas as pd

from scipy import stats, linalg
import bisect
import re
import csv


# from baseP import _show_lineages

#from statannot import add_stat_annotation  #pip install git+https://github.com/webermarcolivier/statannot
from baseP.configs.data_configs import others_Data


"""
Part 1
functions for cancer specific regressin analysis
"""
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
def fetchSig_gene_index(gene_list, exprsn_type, logger, name, data_type):
    if data_type == 'others':
        gene_name = pd.read_csv(others_Data[exprsn_type+"_rownames"])
        gene_name_dict = {item: idx for idx, item in enumerate(gene_name['human_gene_symbol'])}
        
        gene_ind = [gene_name_dict.get(item) for item in gene_list.index.tolist()]
        gene_ind = [x for x in gene_ind if x is not None]
        
        # return {re.sub('_rownames', '', gene_name_file): gene_ind}
        return {exprsn_type: gene_ind}

def fetchSig_val_by_index(gene_ind, exprsn_type, output, logger, name, data_type):
    if data_type == 'others':
        if exprsn_type in gene_ind.keys():
            exprsn_file = os.path.join(others_Data[exprsn_type])
            # if the expression file exists, then extract gene list from it
            if os.path.isfile(exprsn_file):
                
                gene_ind = gene_ind[exprsn_type]
                # append 0 because we need sample name
                gene_ind = [idx+1 for idx in gene_ind] # add 1 because the first row is included in csv.reader results
                gene_ind.append(0)
                
                with open(exprsn_file) as fd:
                    reader = csv.reader(fd)
                    rows = [[np.nan if item == 'NA' else item for item in row] for idx, row in enumerate(reader) if idx in gene_ind]
                    # rows = rows.replace()
                    gene_name = [x[0] for x in rows[1:]]    # extract mouse gene name
                    tissue_name = rows[0][2:]           # the first two columns are mouse_gene_symbol and human_gene_symbol
                    rows = np.array(rows)               # convert list of list to numpy array
                    rows = rows[1:,2:].astype(float)     # remove gene name column and convert item to float type
                
                df_out = pd.DataFrame(data = rows).T
                df_out.columns = gene_name
                
                # convert to log2(1+value)
                df_out += 1
                df_out = df_out.applymap(np.log2)
                
                df_out.insert(0, "cell_line",  tissue_name, True)
                df_out.insert(1, "lineage",  'Neuro2a', True)

                df_out.to_csv(os.path.join(output, 'exprsn', 'tables', exprsn_type+'_query_cell_lines.csv'), index = False)
                # return {cohort: rows,
                #         cohort+'_gene_name': gene_name}









######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########



######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########