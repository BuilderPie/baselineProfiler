#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 06/20/2022
# Last Modified Date: 07/05/2022
# Last Modified By  : Dian Li
# function to process user input variables

import argparse
import logging
import os
import pandas as pd
import baseP.configs.dirpath
import importlib

from baseP import default_cell_line

__doc__="""
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genelist',help="A gene list")
    parser.add_argument('--cell_line',help="A cell line list")
    parser.add_argument('-e', '--evaluators', nargs='+', default=[
                        'CCLE', 'GTEx', 'HPA', 'post'], help="Evaluators to run")
    parser.add_argument('-d', '--data_folder', default='static',help="preprocessed data folder")
    parser.add_argument('-n', '--sig_name', default='SIG',help="Name of input gene list")
    parser.add_argument('-r', '--resum', default=True,
                        help="resume the process")
    parser.add_argument('-t','--threads',default=4,type=int,help='Number of threads')
    parser.add_argument('-o', '--output', default='inputed sig_name',type=str, help="Output directory")
    parser.add_argument('--html_modules', nargs='+', default=[
                        'CCLE', 'GTEx', 'HPA'], help = "modules user selects to view in html")
    parser.add_argument('-p', '--precompute', default=False,
                         help='pre-compute data, input could be (only) ONE of the following modules: symbols, Neuro2a')
    

    args = parser.parse_args()

    #config output directories
    if args.output == 'inputed sig_name':
        args.output = os.path.abspath(os.path.join('analysis',args.sig_name))
    if not os.path.isdir(args.output): #generate output folder if not there
        os.makedirs(args.output)
    request = {
        'output':args.output,
        'input_dir':os.path.abspath(os.path.dirname(args.genelist)),
    }
    print('output directory: '+args.output)

    baseP.configs.dirpath.initialize_dr(args.data_folder)
    
    from baseP import api 
    
    # precompute data
    from baseP.prepare_data import prepare_data
    
    #infer secondary data
    if args.precompute:
        if "symbols".casefold() in args.precompute.casefold():
            # config_logger.info('run pyhton to preprocess CCLE data to generate biomarker, infiltration and so on')
            print('Annotation symbols precompute started')
            prepare_data.symbols_precompute()
            # importlib.reload(baseP.configs.data_configs)
            # from baseP.configs.data_configs import CCLE_Data
            print('Annotation symbols precompute finished')
        if "Neuro2a".casefold() in args.precompute.casefold():
            print('Annotation symbols precompute started')
            prepare_data.Neuro2a_precompute()  
            print('Annotation symbols precompute finished')
        
    # set up logging to file
    
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='%s/baseP.log' % args.output,
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')     # set a format which is simpler for console use
    console.setFormatter(formatter)     # tell the handler to use this format
    logging.getLogger('').addHandler(console)  # add the handler to the root logger

    ########################################################################
    #parse the input gene list, append to data
    data = []
    args.genelist = os.path.abspath(os.path.join(args.genelist))
    # weighted boolean variable indicates if figures should be split for pos and neg-weighted genes in the Crispr html module
    weighted = False

    # weight dictionary indicates presence of positively and negatively weighted genes
    weights = {"positive":False, "negative":False}

    with  open(args.genelist) as file:
        for line in file:
            tmp  = line.strip().split()
            data.append(tmp)
    if len(data[0]) == 1: #unweighted input gene list
        data = [x[0] for x in data]
        data.sort()
        request['data'] = data

     
    else: #weighted input gene list
        # set "weighted" variable to True to split Crispr figures for pos and neg-weighted genes
        # update the weight dictionary values according to presence of positive and negative weights
        weighted = True
        
        pos_weights = [x[1] for x in data if float(x[1]) >= 0]
        neg_weights = [x[1] for x in data if float(x[1]) < 0]
        if len(pos_weights) > 0:
            #print("positive weights found")
            weights["positive"] = True
        if len(neg_weights) > 0:
            #print("negative weights found")
            weights["negative"] = True      
        request['data'] = pd.Series([float(x[1]) for x in data],index=[x[0] for x in data],name='weight')
    
    ########################################################################
    # if args.cell_line is specified
    # parse the input cell line list, append to cellLine
    
    if args.cell_line is not None:
        
        cell_line = []
        args.cell_line = os.path.abspath(os.path.join(args.cell_line))

        with  open(args.cell_line) as file:
            for line in file:
                tmp  = line.strip()
                cell_line.append(tmp)
        request['cell_line'] = [x for x in cell_line]
        # if len(cell_line[0]) == 1: #unweighted input gene list
        #     request['cell_line'] = [x[0] for x in cell_line]
    else:
        # if args.cell_line is not specified, use the default_cell_line 
        request['cell_line'] = default_cell_line
    
    logging.info('Successfully Digested Input data and Initialize evaluation...')
    geneset_logger = logging.getLogger('[Geneset]')     #call api and pass the parameters
    # api.geneset(request=request, logger_P=geneset_logger,
    #             evaluators=args.evaluators, name=args.sig_name, resum=args.resum, threads=args.threads, html_modules = args.html_modules, weights = weights)
    api.geneset(request=request, logger_P=geneset_logger,
                evaluators=args.evaluators, name=args.sig_name, html_modules = args.html_modules,
                resum=args.resum, threads=args.threads, weights = weights)
    print('All results are at: ', os.path.abspath(args.output ))

