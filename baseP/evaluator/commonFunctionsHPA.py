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


# from baseP import HPA_show_lineages

#from statannot import add_stat_annotation  #pip install git+https://github.com/webermarcolivier/statannot
from baseP.configs.data_configs import HPA_Data


"""
Part 1
functions for cancer specific regressin analysis
"""
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
def fetchSig_gene_index(gene_list, exprsn_type, logger, name, data_type):
    if data_type == 'HPA':
        gene_name = pd.read_csv(HPA_Data[exprsn_type+"_rownames"])
        gene_name_dict = {item: idx for idx, item in enumerate(gene_name['gene'])}
        
        gene_ind = [gene_name_dict.get(item) for item in gene_list.index.tolist()]
        gene_ind = [x for x in gene_ind if x is not None]
        
        # return {re.sub('_rownames', '', gene_name_file): gene_ind}
        return {exprsn_type: gene_ind}

def fetchSig_val_by_index(gene_ind, exprsn_type, output, logger, name, data_type):
    if data_type == 'HPA':
        if exprsn_type in gene_ind.keys():
            exprsn_file = os.path.join(HPA_Data[exprsn_type])
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

                df_tissue.insert(0, "cell_line",  tissue_name, True)
                # df_tissue.insert(1, "lineage",  cohort, True)
                return {exprsn_type: df_tissue,
                        exprsn_type+'_gene_name': gene_name}
                # df_tissue.to_csv(os.path.join(output, exprsn_type, 'tables', cohort+'.csv'), index = False)
                # return {cohort: rows,
                #         cohort+'_gene_name': gene_name}



def match_lineage(sig_dict, cell_line, exprsn_type, output, logger, name, data_type):
    if data_type == 'HPA':
        # colnames = pd.read_csv(HPA_Data[exprsn_type+"_colnames"]) # column names / cell line names of the expression matrix
        # df = pd.DataFrame(data=sig_dict[exprsn_type], index = sig_dict[exprsn_type+'_gene_name']) # expression matrix to be splitted
        df = sig_dict[exprsn_type]
        ##############======================================#################
        ##############======================================#################
        # part 1. cell line specific
        # if there is cell_line specified, create a subset of expression matrix based on query cell lines
        if isinstance(cell_line, list):
            cell_line = [x.upper() for x in cell_line]
            ##############======================================#################
            if exprsn_type in ['Exprsn_cell_line']:
                metaFile = pd.read_csv(HPA_Data[exprsn_type+"_meta"])  # meta information
                # index of query cell line in the meta table
                idx_0 = np.where(np.isin(metaFile['HPA_cell_line'].tolist(), cell_line))[0]
                idx_1 = np.where(np.isin(metaFile['CCLE_cell_line_name'].tolist(), cell_line))[0]
                idx_2 = np.where(np.isin(metaFile['CCLE_stripped_cell_line_name'].tolist(), cell_line))[0]

                idx_union = np.union1d(idx_0, np.union1d(idx_1, idx_2))
                
                # index in column stripped_cell_line_name without the overlapped ones
                # idx_2_unique = idx_2[np.in1d(idx_2, idx_1, invert=True)]
                # concatanate two dataframes in the row axis
                # meta_subset = pd.concat([metaFile.loc[idx_1,['HPA_cell_line', 'lineage']],
                # metaFile.loc[idx_2_unique,['DepMap_ID', 'stripped_cell_line_name', 'lineage']].rename(columns={'stripped_cell_line_name':'cell_line_name'})], ignore_index=True)
                meta_subset = metaFile.loc[idx_union,['HPA_cell_line', 'lineage']]
                
                # index of query cell line in the expression file
                # idx_cell_line = np.where(np.isin(colnames['cell_line'].tolist(), meta_subset['DepMap_ID'].tolist()))[0]
                idx_cell_line = np.where(np.isin(df['cell_line'].tolist(), meta_subset['HPA_cell_line'].tolist()))[0]
                
                # create a dataframe with matched DepMap_ID and corresponding column index in the expression file 
                meta_matched = pd.DataFrame(data={'HPA_cell_line': df.loc[idx_cell_line, 'cell_line'], 
                'col_index': idx_cell_line})

                # left join two dataframes meta_matched and meta_subset, to get the final meta_matched that 
                # provides the index of query cell lines in the expression file 
                meta_matched = meta_matched.merge(meta_subset, on='HPA_cell_line', how='left')
                
                df_query_cell_line = df.iloc[meta_matched['col_index'], :]
                # add cell_line name and lineage columns
                # df_query_cell_line.insert(0, "cell_line",  meta_matched['cell_line_name'].tolist(), True)
                
                df_query_cell_line.insert(1, "lineage",  meta_matched['lineage'].tolist(), True)
                # sort dataframe by cell_line name
                df_query_cell_line = df_query_cell_line.sort_values(by=['cell_line'])
                
                df_query_cell_line.to_csv(os.path.join(output, 'Exprsn', 'tables', 'query_cell_lines.csv'), index = False)
            ##############======================================#################

        ##############======================================#################
        ##############======================================#################
        # part 2. export combined-lineage files
        if exprsn_type in ['Exprsn_cell_line']:
            metaFile = pd.read_csv(HPA_Data[exprsn_type+"_meta"])  # meta information
            meta_subset = metaFile.loc[:, ['HPA_cell_line', 'lineage']]
            idx_cell_line = np.where(np.isin(df['cell_line'].tolist(), meta_subset['HPA_cell_line'].tolist()))[0]

            meta_matched = pd.DataFrame(data={'HPA_cell_line': df.loc[idx_cell_line, 'cell_line'], 
                'col_index': idx_cell_line})
            meta_matched = meta_matched.merge(meta_subset, on='HPA_cell_line', how='left')
            
            df_combined = df.iloc[meta_matched['col_index'], :]
            # add cell_line name and lineage columns
            # df_query_cell_line.insert(0, "cell_line",  meta_matched['cell_line_name'].tolist(), True)
            
            df_combined.insert(1, "lineage",  meta_matched['lineage'].tolist(), True)
            # sort dataframe by cell_line name
            df_combined = df_combined.sort_values(by=['cell_line'])
            
            df_combined.to_csv(os.path.join(output, 'Exprsn', 'tables', 'all_cell_lines.csv'), index = False)
        elif exprsn_type in ['Exprsn_blood_cell']:
            
            df.to_csv(os.path.join(output, 'Exprsn', 'tables', 'all_blood_cells.csv'), index = False)





######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########



######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########