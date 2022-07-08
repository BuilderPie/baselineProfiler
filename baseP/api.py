#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 06/23/2022
# Last Modified Date: 07/05/2022
# Last Modified By  : Dian Li
# -*- coding: utf-8 -*-

import os
import glob
import subprocess
import pandas as pd
# from Biopyutils import Comm
# import baseP.utils.pyUtils as utils
# from baseP.evaluator import CCLE 
import baseP.evaluator.CCLE as CCLE 
import baseP.evaluator.GTEx as GTEx
import baseP.evaluator.HPA as HPA
import baseP.evaluator.others as others
# from baseP import R_RMD
from baseP.utils import postProcess


from baseP.configs.data_configs import DATA_DIR

from baseP.configs.folder_configs import formOutput

# def geneset(request,logger_P, name, threads, evaluators, resum, html_modules, weights):
def geneset(request,logger_P, name, threads, evaluators, html_modules, resum, weights):
    ''' Evaluate gene signature composite of a list of gene '''
    
    formOutput(request['output'])
    
    #create all folders
    #based on the input gene list, assign the species
    species = "hg38"
    df_in=request['data']
    if isinstance(df_in,list):
        df_in = pd.Series(index=df_in,name='weight')
    if df_in.index[0].startswith('ENSMUSG') or  (df_in.index[0][0].isupper() and df_in.index[0][1:].islower() ):
        species='mm10'
    else:
        data = df_in.copy()
        for k,v in {'hg':data}.items():
            outfile = os.path.join(request['output'],'input_geneset_%s.txt'% k)
            if isinstance(v,list):
                with open(outfile,'w') as f:
                    for item in v:
                        f.write("%s\n" % item)
            else:
                v.to_csv(outfile,sep='\t',header=False)
    
    # write cellLine to the output folder 
    cell_line = request['cell_line']
    if len(cell_line) > 0:
        outfile = os.path.join(request['output'],'input_cell_line.txt')
        if isinstance(cell_line,list):
            with open(outfile,'w') as f:
                for item in cell_line:
                    f.write("%s\n" % item)
        else:
            cell_line.to_csv(outfile,sep='\t',header=False)
    #start evaluating the gene signature in different data modules
    #################### R CMD for CRISPR and scRNAICB and scRNANonICB####################
    # inputGeneList_mouse = os.path.join(
    #     request['output'], 'input_geneset_mm.txt')   
    inputGeneList_human = os.path.join(
        request['output'], 'input_geneset_hg.txt')
    
    # if 'scICB' in evaluators:
    #     logger = logger_P.getChild('[ICBscRNA]')
    #     logger.info(
    #         'Start evaluating association of input geneset with single-cell ICB data.')
    #     r_cmd = ' '.join(['Rscript', R_RMD, '--name', ('\''+name+'\''), '--file_mouse ', inputGeneList_mouse.replace(' ', '\ '), '--file_human ', inputGeneList_human.replace(' ', '\ '), '-m', 'scICB',
    #                       '-d', DATA_DIR, '-o', os.path.join(request['output'], 'scICB_evaluator').replace(' ', '\ '), '-t', str(threads)])
    #     print(r_cmd)
    #     process = subprocess.Popen(r_cmd, shell=True).wait()
    #     print(process)  # [0] is stdout
    #     logger_P.info('Successfully finish single-cell ICB evaluator.')

    if 'CCLE' in evaluators:
        logger = logger_P.getChild('[CCLE]')
        logger.info('Start evaluating association of input gene list in CCLE data.')
        CCLE.analyze(data, cell_line, output=os.path.join(request['output'],'CCLE_evaluator'), logger=logger, name=name, threads= threads)
        logger_P.info('Successfully finish CCLE evaluator.')
    
    if 'GTEx' in evaluators:
        logger = logger_P.getChild('[GTEx]')
        logger.info('Start evaluating association of input gene list in GTEx data.')
        GTEx.analyze(data, cell_line, output=os.path.join(request['output'],'GTEx_evaluator'), logger=logger, name=name, threads= threads)
        logger_P.info('Successfully finish GTEx evaluator.')
    
    if 'HPA' in evaluators:
        logger = logger_P.getChild('[HPA]')
        logger.info('Start evaluating association of input gene list in HPA data.')
        HPA.analyze(data, cell_line, output=os.path.join(request['output'],'HPA_evaluator'), logger=logger, name=name, threads= threads)
        logger_P.info('Successfully finish HPA evaluator.')
    
    if 'others' in evaluators:
        logger = logger_P.getChild('[others]')
        logger.info('Start evaluating association of input gene list in others data.')
        others.analyze(data, cell_line, output=os.path.join(request['output'],'others_evaluator'), logger=logger, name=name, threads= threads)
        logger_P.info('Successfully finish others evaluator.')


    # if 'Crispr' in evaluators:
    #     os.environ['MKL_THREADING_LAYER'] = 'GNU'  # test- Dian
    #     logger = logger_P.getChild('[Crispr]')
    #     logger.info(
    #         'Start evaluating association of input geneset with Immune-related CRISPR evaluator')
    #     # outputDire=request['output']
    #     r_cmd = ' '.join(['Rscript', R_RMD, '--name', ('\''+name+'\''), '--file_mouse ', inputGeneList_mouse.replace(' ', '\ '), '--file_human ', inputGeneList_human.replace(' ', '\ '), '-m', 'Crispr', '-d', CRISPR_Data,
    #                       '-o', os.path.join(
    #                           request['output'], 'Crispr_evaluator').replace(' ', '\ '),
    #                       '--conda_env ', os.environ['CONDA_PREFIX']])
    #     process = subprocess.Popen(r_cmd, shell=True).wait()
    #     print(process)

         

    if ('post' in evaluators) or ('Post' in evaluators) :
        logger = logger_P.getChild('[Post Process]')
        logger.info('Start output data post processing')
        postProcess.analyze(data, request=request, logger=logger, name=name, threads=threads, modules_selected=html_modules, weights = weights)
        logger_P.info('Successfully finish data post processing')
