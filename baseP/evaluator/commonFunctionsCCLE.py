#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 06/27/2022
# Last Modified Date: 06/27/2022
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


# from baseP import CCLE_show_lineages

#from statannot import add_stat_annotation  #pip install git+https://github.com/webermarcolivier/statannot
from baseP.configs.data_configs import CCLE_Data


"""
Part 1
functions for cancer specific regressin analysis
"""
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
def fetchSig_gene_index(gene_list, exprsn_type, logger, name, data_type):
    if data_type == 'CCLE':
        gene_name = pd.read_csv(CCLE_Data[exprsn_type+"_rownames"])
        gene_name_dict = {item: idx for idx, item in enumerate(gene_name['gene'])}
        
        gene_ind = [gene_name_dict.get(item) for item in gene_list.index.tolist()]
        gene_ind = [x for x in gene_ind if x is not None]
        
        # return {re.sub('_rownames', '', gene_name_file): gene_ind}
        return {exprsn_type: gene_ind}

def fetchSig_val_by_index(gene_ind, exprsn_type, logger, name, data_type):
    if data_type == 'CCLE':
        if exprsn_type in gene_ind.keys():
            exprsn_file = CCLE_Data[exprsn_type]
            gene_ind = gene_ind[exprsn_type]
            
            gene_ind = [idx+1 for idx in gene_ind] # add 1 because the first row is included in csv.reader results

            with open(exprsn_file) as fd:
                reader = csv.reader(fd)
                rows = [[np.nan if item == 'NA' else item for item in row] for idx, row in enumerate(reader) if idx in gene_ind]
                # rows = rows.replace()
                gene_name = [x[0] for x in rows]    # extract gene name
                rows = np.array(rows)               # convert list of list to numpy array
                rows = rows[:,1:].astype(float)     # remove gene name column and convert item to float type
                
            return {exprsn_type: rows,
                    exprsn_type+'_gene_name': gene_name}

def splitSig_by_lineage(sig_dict, cell_line, exprsn_type, output, logger, name, data_type, show_lineages):
    if data_type == 'CCLE':
        colnames = pd.read_csv(CCLE_Data[exprsn_type+"_colnames"]) # column names / cell line names of the expression matrix
        df = pd.DataFrame(data=sig_dict[exprsn_type], index = sig_dict[exprsn_type+'_gene_name']) # expression matrix to be splitted
        
        metaFile = pd.read_csv(CCLE_Data[exprsn_type+"_meta"])  # meta information
        ##############======================================#################
        ##############======================================#################
        # part 1. cell line specific
        # if there is cell_line specified, create a subset of expression matrix based on query cell lines
        if isinstance(cell_line, list):
            cell_line = [x.upper() for x in cell_line]
            ##############======================================#################
            if exprsn_type in ['Exprsn', 'CRISPR_Broad']:
                # index of query cell line in the meta table
                idx_1=np.where(np.isin(metaFile['cell_line_name'].tolist(), cell_line))[0]
                idx_2=np.where(np.isin(metaFile['stripped_cell_line_name'].tolist(), cell_line))[0]
                # index in column stripped_cell_line_name without the overlapped ones
                idx_2_unique = idx_2[np.in1d(idx_2, idx_1, invert=True)]
                # concatanate two dataframes in the row axis
                meta_subset = pd.concat([metaFile.loc[idx_1,['DepMap_ID', 'cell_line_name', 'lineage']],
                metaFile.loc[idx_2_unique,['DepMap_ID', 'stripped_cell_line_name', 'lineage']].rename(columns={'stripped_cell_line_name':'cell_line_name'})], ignore_index=True)
                
                # index of query cell line in the expression file
                idx_cell_line = np.where(np.isin(colnames['cell_line'].tolist(), meta_subset['DepMap_ID'].tolist()))[0]
                
                # create a dataframe with matched DepMap_ID and corresponding column index in the expression file 
                meta_matched = pd.DataFrame(data={'DepMap_ID': colnames.loc[idx_cell_line, 'cell_line'], 
                'col_index': idx_cell_line})

                # left join two dataframes meta_matched and meta_subset, to get the final meta_matched that 
                # provides the index of query cell lines in the expression file 
                meta_matched = meta_matched.merge(meta_subset, on='DepMap_ID', how='left')

                df_query_cell_line = df.iloc[:, meta_matched['col_index']].T
                # add cell_line name and lineage columns
                df_query_cell_line.insert(0, "cell_line",  meta_matched['cell_line_name'].tolist(), True)
                df_query_cell_line.insert(1, "lineage",  meta_matched['lineage'].tolist(), True)
                # sort dataframe by cell_line name
                df_query_cell_line = df_query_cell_line.sort_values(by=['cell_line'])
                
                df_query_cell_line.to_csv(os.path.join(output, exprsn_type, 'tables', 'query_cell_lines.csv'), index = False)
            ##############======================================#################
            
            elif exprsn_type in ['Proteomics']:
                # index of query cell line in the meta table
                idx_1=np.where(np.isin(metaFile['cell_line'].tolist(), cell_line))[0]
                idx_2=np.where(np.isin(metaFile['cell_line_CCLE'].tolist(), cell_line))[0]
                # index in column stripped_cell_line_name without the overlapped ones
                idx_2_unique = idx_2[np.in1d(idx_2, idx_1, invert=True)]
                # concatanate two dataframes in the row axis
                meta_subset = pd.concat([metaFile.loc[idx_1,['cell_line_CCLE', 'lineage']],
                metaFile.loc[idx_2_unique,['cell_line_CCLE', 'lineage']]], ignore_index=True)
                
                # index of query cell line in the expression file
                idx_cell_line = np.where(np.isin(colnames['cell_line'].tolist(), meta_subset['cell_line_CCLE'].tolist()))[0]
                # create a dataframe with matched DepMap_ID and corresponding column index in the expression file 
                meta_matched = pd.DataFrame(data={'cell_line_CCLE': colnames.loc[idx_cell_line, 'cell_line'], 
                'col_index': idx_cell_line})
                # left join two dataframes meta_matched and meta_subset, to get the final meta_matched that 
                # provides the index of query cell lines in the expression file 
                meta_matched = meta_matched.merge(meta_subset, on='cell_line_CCLE', how='left')
                # create the final output for the query cell line and gene list dataframe
                df_query_cell_line = df.iloc[:, meta_matched['col_index']].T
                # add cell_line name and lineage columns
                df_query_cell_line.insert(0, "cell_line",  meta_matched['cell_line_CCLE'].tolist(), True)
                df_query_cell_line.insert(1, "lineage",  meta_matched['lineage'].tolist(), True)
                # sort dataframe by cell_line name
                df_query_cell_line = df_query_cell_line.sort_values(by=['cell_line'])
                # write df_query_cell_line to csv table
                df_query_cell_line.to_csv(os.path.join(output, exprsn_type, 'tables', 'query_cell_lines.csv'), index = False)
        ##############======================================#################
        ##############======================================#################
        # part 2. split by lineages
        if exprsn_type in ['Exprsn', 'CRISPR_Broad']:
            for item in show_lineages[exprsn_type]:
                idx_1=np.where(np.isin(metaFile['lineage'].tolist(), item))[0]
                meta_subset = metaFile.loc[idx_1,['DepMap_ID', 'cell_line_name', 'lineage']]
                
                # find the indices  of query lineage in the expression file
                idx_lineage = np.where(np.isin(colnames['cell_line'].tolist(), meta_subset['DepMap_ID'].tolist()))[0]
                if len(idx_lineage)>0:
                    # create a dataframe with matched DepMap_ID and corresponding column index in the expression file 
                    meta_matched = pd.DataFrame(data={'DepMap_ID': colnames.loc[idx_lineage, 'cell_line'], 'col_index': idx_lineage})
                    # left join two dataframes meta_matched and meta_subset, to get the final meta_matched that 
                    # provide the index of query cell lines in the expression file 
                    meta_matched = meta_matched.merge(meta_subset, on='DepMap_ID', how='left')
                    
                    df_query_lineage = df.iloc[:, meta_matched['col_index']].T
                    # add cell_line name and lineage columns
                    df_query_lineage.insert(0, "cell_line",  meta_matched['cell_line_name'].tolist(), True)
                    df_query_lineage.insert(1, "lineage",  meta_matched['lineage'].tolist(), True)
                    # write df_query_lineage to csv table
                    df_query_lineage.to_csv(os.path.join(output, exprsn_type, 'tables', item+'.csv'), index = False)
        #############======================================#################
        elif exprsn_type in ['Proteomics']:
            for item in show_lineages[exprsn_type]:
                idx_1=np.where(np.isin(metaFile['lineage'].tolist(), item))[0]
                meta_subset = metaFile.loc[idx_1,['cell_line_CCLE', 'lineage']]

                # find the indices of query lineage in the expression file
                idx_lineage = np.where(np.isin(colnames['cell_line'].tolist(), meta_subset['cell_line_CCLE'].tolist()))[0]
                if len(idx_lineage)>0:
                    # create a dataframe with matched DepMap_ID and corresponding column index in the expression file 
                    meta_matched = pd.DataFrame(data={'cell_line_CCLE': colnames.loc[idx_lineage, 'cell_line'], 'col_index': idx_lineage})
                    # left join two dataframes meta_matched and meta_subset, to get the final meta_matched that 
                    # provide the index of query cell lines in the expression file 
                    meta_matched = meta_matched.merge(meta_subset, on='cell_line_CCLE', how='left')

                    df_query_lineage = df.iloc[:, meta_matched['col_index']].T
                    # add cell_line name and lineage columns
                    df_query_lineage.insert(0, "cell_line",  meta_matched['cell_line_CCLE'].tolist(), True)
                    df_query_lineage.insert(1, "lineage",  meta_matched['lineage'].tolist(), True)
                    # write df_query_lineage to csv table
                    df_query_lineage.to_csv(os.path.join(output, exprsn_type, 'tables', item+'.csv'), index = False)







######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########



######## =========================================================== ########
######## =========================================================== ########
######## =========================================================== ########