#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 06/30/2022
# Last Modified Date: 06/30/2022
# Last Modified By  : Dian Li  
# -*- coding:utf-8 -*-
import os
from os import listdir
import joblib
import warnings
warnings.simplefilter('ignore')
# import numpy as np
import pandas as pd
import json
import subprocess
from baseP.configs.data_configs import HTML_Data
from baseP.configs.data_configs import HTML_DIR
from baseP.utils.htmlReadTemplate import *
import shutil
from shutil import copyfile

def analyze(df,request,logger,name,threads, modules_selected, weights):
	"""
	Pipeline to do post analysis data processing
	Parameters
	df : list / pd.Series
		Weighted or un-weighted gene list
	"""
	# #====================#
	# logger.info('Start: Convert tables to JSON format')
	# CCLE_tables = os.path.join(output, 'CCLE_evaluator/tables')
	# CCLE_tables_JSON = os.path.join(output, 'CCLE_evaluator/tables/JSON')
	# meta_tables = os.path.join(HTML_Data['meta'], 'tables')
	# meta_tables_JSON = os.path.join(output, 'meta/tables/JSON')


	# #====================#
# 	file_names = [fn for fn in os.listdir(CCLE_tables)
# 				  if any(fn.endswith(ext) for ext in included_extensions)]
# 	for fn in file_names:
# 		txt_2_json(fn=fn, input_path=CCLE_tables, output_path=CCLE_tables_JSON)
# 	#====================#

# 	#====================#
# 	logger.info('Finish: Convert tables to JSON format')
# 	#====================#
	# copy GTEx figures to heatmp figs folder
	src = os.path.join(request['output'], 'GTEx_evaluator', 'Exprsn', 'plots')
	dst = os.path.join(request['output'], 'html', 'figs', 'GTEx_evaluator','Exprsn')
	shutil.rmtree(dst)
	shutil.copytree(src, dst) 
	#====================#
	# copy HPA figures to heatmp figs folder
	src = os.path.join(request['output'], 'HPA_evaluator', 'Exprsn', 'plots')
	dst = os.path.join(request['output'], 'html', 'figs', 'HPA_evaluator','Exprsn')
	shutil.rmtree(dst)
	shutil.copytree(src, dst) 

	#====================#
	# copy CCLE Exprsn figures to heatmp figs folder
	src = os.path.join(request['output'], 'CCLE_evaluator', 'Exprsn', 'plots')
	dst = os.path.join(request['output'], 'html', 'figs', 'CCLE_evaluator','Exprsn')
	shutil.rmtree(dst)
	shutil.copytree(src, dst) 

	#====================#
	# copy CCLE Proteomics figures to heatmp figs folder
	src = os.path.join(request['output'], 'CCLE_evaluator', 'Proteomics', 'plots')
	dst = os.path.join(request['output'], 'html', 'figs', 'CCLE_evaluator','Proteomics')
	shutil.rmtree(dst)
	shutil.copytree(src, dst) 

	#====================#
	# copy CCLE CRISPR_Broad figures to heatmp figs folder
	src = os.path.join(request['output'], 'CCLE_evaluator', 'CRISPR_Broad', 'plots')
	dst = os.path.join(request['output'], 'html', 'figs', 'CCLE_evaluator','CRISPR_Broad')
	shutil.rmtree(dst)
	shutil.copytree(src, dst) 

# 	# try: 
	logger.info('Start: Generate HTML report')
	html_build(analysis_path=request['output'],template_path=HTML_DIR,output_path=os.path.join(request['output'],'html'),name=name, modules_selected=modules_selected)
# 	cp_cmd = ''.join(['cp ', HTML_Data['CSS'], ' ', os.path.join(request['output'],'html')])
# 	print(cp_cmd)
# 	process=subprocess.Popen(cp_cmd, shell=True).wait()
# 	logger.info('Finish: Generate HTML report')
# 	# except Exception as e: print(e)	
# 	#========= copy html and css file to upper level folder ===========#
# 	cp_cmd = ''.join(['cp ', os.path.join(request['output'], 'html/* '), os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(request['output']))), '')])
# 	print(cp_cmd)
# 	process=subprocess.Popen(cp_cmd, shell=True).wait()
# 	logger.info('Finish: Sent analysis resutls to Karen TISET folder')
# 	#====================#
	

	


# 	logger.info('Finish: Post process')


# def txt_2_json(fn,input_path,output_path):
# 	df = pd.read_csv(os.path.join(input_path, fn), sep = '\t', na_values="")
# 	df = df.loc[:,~df.columns.str.contains('^Unnamed')]
# 	col_tmp = df.columns
# 	# col_tmp = df.columns[~df.columns.str.contains('^Unnamed')]
# 	df = df.fillna("")
# 	df = df.dropna(axis=1, how='all')
# 	df.columns = col_tmp

# 	for col in df.columns:
# 		if df[col].dtypes == "float64":
# 			df_temp = df[col][:]
# 			df_temp[abs(df[col])<0.001] = df[col][abs(df[col])<0.001].map(lambda x: '%2.1e' % x)
# 			df_temp[(abs(df[col])>=0.001)&(abs(df[col])<1)] = df[col][(abs(df[col])>=0.001)&(abs(df[col])<1)].map(lambda x: '%2.2g' % x)
# 			df_temp[abs(df[col])>=1] = df[col][abs(df[col])>=1].map(lambda x: '%2.1f' % x)

# 			df[col] = df_temp
# 		if 'HIGH' in col.upper():
# 			df = df.rename(columns={col: "Sig High"})
# 		if 'LOW' in col.upper():
# 			df = df.rename(columns={col: "Sig Low"})
# 	df_json = json.dumps(json.loads(df.to_json(orient='records')), indent=2)
# 	df_json = '{ \n' + '  "' + 'data' + '"' + ': '  + df_json + '\n}'
	
# 	df_json = df_json.replace('_', ' ')
# 	df_json = df_json.replace(';', '; ')

# 	with open(os.path.join(output_path, fn), 'w') as outfile:
# 		outfile.write(df_json)
