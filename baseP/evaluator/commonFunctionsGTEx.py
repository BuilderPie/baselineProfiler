#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 06/29/2022
# Last Modified Date: 06/29/2022
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


from baseP import GTEx_show_lineages

#from statannot import add_stat_annotation  #pip install git+https://github.com/webermarcolivier/statannot
from baseP.configs.data_configs import GTEx_Data


"""
Part 1
functions for cancer specific regressin analysis
"""
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
def fetchSig_gene_index(gene_list, exprsn_type, logger, name, data_type):
    if data_type == 'GTEx':
        gene_name = pd.read_csv(GTEx_Data[exprsn_type+"_rownames"])
        gene_name_dict = {item: idx for idx, item in enumerate(gene_name['gene'])}
        
        gene_ind = [gene_name_dict.get(item) for item in gene_list.index.tolist()]
        gene_ind = [x for x in gene_ind if x is not None]
        
        # return {re.sub('_rownames', '', gene_name_file): gene_ind}
        return {exprsn_type: gene_ind}

def fetchSig_val_by_index(gene_ind, exprsn_type, cohort, output, logger, name, data_type):
    if data_type == 'GTEx':
        if exprsn_type in gene_ind.keys():
            exprsn_file = os.path.join(GTEx_Data[exprsn_type], cohort+'.csv')
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
                    gene_name = [x[0] for x in rows[1:]]    # extract gene name
                    tissue_name = rows[0][1:]
                    rows = np.array(rows)               # convert list of list to numpy array
                    rows = rows[1:,1:].astype(float)     # remove gene name column and convert item to float type
                
                df_tissue = pd.DataFrame(data = rows).T
                df_tissue.columns = gene_name

                # convert to log2(1+value)
                df_tissue += 1
                df_tissue = df_tissue.applymap(np.log2)


                
                df_tissue.insert(0, "cell_line",  tissue_name, True)
                df_tissue.insert(1, "lineage",  cohort, True)
                
                df_tissue.to_csv(os.path.join(output, exprsn_type, 'tables', cohort+'.csv'), index = False)
                # return {cohort: rows,
                #         cohort+'_gene_name': gene_name}









######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########



######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########