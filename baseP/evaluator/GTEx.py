#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li <dianli@wustl.edu> Jingxin Fu
# Date              : 06/29/2022
# Last Modified Date: 06/29/2022
# Last Modified By  : Dian Li <dianli@wustl.edu>
# -*- coding:utf-8 -*-
import os
from os import listdir, stat
# import joblib
import warnings
warnings.simplefilter('ignore')
# try:
#    import cPickle as pickle
# except:
#    import pickle

import numpy as np
import pandas as pd

# from scipy import stats, linalg
# import lifelines
# from lifelines import CoxPHFitter
# import statsmodels.stats.multitest as multi
# from sklearn.preprocessing import LabelEncoder
# import bisect

import joblib

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import seaborn as sns

#from statannot import add_stat_annotation  #pip install git+https://github.com/webermarcolivier/statannot
from baseP.configs.data_configs import GTEx_Data
from baseP import GTEx_show_lineages

from .commonFunctionsGTEx import * # functions for each CCLE module

# from multiprocessing import Pool
import functools

def analyze(gene_list, cell_line, output,logger,name,threads):
    """
    Pipeline to analyze input signature on current GTEx cancer types and functional datasets
    
    ----------
    Parameters
    df : list / pd.Series
        Weighted or un-weighted gene list
    output : str
        output directory
    name : str
        Name of input signature [Default: 'SIG']
    threads : int, optional
        number of threads, by default 4
    
    Returns
    ----------
    """
    exprsn_type_list = ['Exprsn']
    

    #### task list 1. Perform correlation analysis on single cohort / cancer type in Compound screeing, ####
    #### Protein array, Proteomics, CRISPR screen ####
    
    tasks_expr_single = {
            'signature expr rna_seq':{
                'calculation': 'sigExprRNA',
                },
            }

    #### task list 2. Perform regression analysis on all the cohorts / cancer types in Compound screeing, ####
    #### Protein array, Proteomics, CRISPR screen ####
    
    # tasks_corr_all = {
    #         'signature corr Proteomics':{
    #             'process': 'sigExprRNA',
    #             },
            
    #         }
    #################### 1. Fetch signature gene index (row index) from GTEx datasets ####################
    ind_list = joblib.Parallel(n_jobs=threads, backend='threading')(joblib.delayed(fetchSig_gene_index)(
                 gene_list = gene_list,
                 exprsn_type = exprsn_type,
                 logger = logger,
                 name = name,
                 data_type = 'GTEx') for exprsn_type in exprsn_type_list )
                 
    ind_dict = {}
    for x in ind_list:
        ind_dict.update(x)
    
    #################### 2. Fetch signature expression values by gene index (row index) in GTEx datasets and write to files ########
    sig_list = joblib.Parallel(n_jobs=threads, backend='threading')(joblib.delayed(fetchSig_val_by_index)(
                 gene_ind = ind_dict,
                 exprsn_type = exprsn_type_list[0],
                 cohort = cohort,
                 output = output,
                 logger = logger,
                 name = name,
                 data_type = 'GTEx') for cohort in GTEx_show_lineages[exprsn_type_list[0]] )
    
    # sig_dict = {}
    # for x in sig_list:
    #     if x is not None:
    #         sig_dict.update(x)
    

    
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########


