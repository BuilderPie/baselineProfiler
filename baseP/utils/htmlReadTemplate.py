#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li <dianli@wustl.edu> Jingxin Fu
# Date              : 06/30/2022
# Last Modified Date: 06/30/2022
# Last Modified By  : Dian Li <dianli@wustl.edu>
# -*- coding:utf-8 -*-


import os
from os import listdir
from os import path
import joblib
import warnings
warnings.simplefilter('ignore')
import numpy as np
import pandas as pd
import json
import subprocess
from baseP.configs.data_configs import HTML_Data
from baseP.utils.utilsPlotly import plotly_heatmap, plotly_heatmap_wo_dendrogram
import datetime as dt

from baseP import CCLE_show_lineages
from baseP import GTEx_show_lineages

import shutil

def html_build(analysis_path,template_path,output_path,name, modules_selected):
	
	# copy and paste library directory and css files to the destination folder
	dst = os.path.join(output_path, 'site_libs')
	src = os.path.join(template_path, 'site_libs')
	
	# try:
	# 	shutil.rmtree(dst)
	# except:
	# 	print(dst + " doesn't exist")
	if not os.path.isdir(dst):
		shutil.copytree(src, dst) 

	shutil.copy(os.path.join(template_path, 'style.css'), output_path)
	shutil.copy(os.path.join(template_path, '.nojekyll'), output_path)
	shutil.copy(os.path.join(template_path, 'index.html'), output_path)

	if "CCLE" in modules_selected:
		exprsn_type_list = ['Exprsn', 'Proteomics', 'CRISPR_Broad']
		# if the corresponding template exists
		for exprsn_type in exprsn_type_list:
			template_file = os.path.join(template_path, "CCLE_"+exprsn_type+'.html')
			if os.path.isfile(template_file):
				# create a replace dictionary that store the plotly html code
				# to be replaced into the template
				replace_dict = {}
				counter = 0
				for lineage in CCLE_show_lineages[exprsn_type]:
					lineage_csv = os.path.join(analysis_path, "CCLE_evaluator", exprsn_type, "tables", lineage+".csv")
					if os.path.isfile(lineage_csv):
						df = pd.read_csv(lineage_csv, index_col = 0)
						df = df.iloc[:,1:]

						df_index = df.index.tolist()
						nan_index = pd.isnull(df_index)
						if sum(nan_index) > 0:
							nan_index = [i for i, x in enumerate(nan_index) if x]
							for x in nan_index:
								df_index[x] = "NAN_"+str(x)
							df.index = df_index
						include_plotlyjs = True if counter == 0 else 'cdn'
						if exprsn_type == "Exprsn":
							fig_CCLE = plotly_heatmap(df = df)
							replace_dict['{{'+exprsn_type+'_'+lineage+'}}'] = fig_CCLE.to_html(full_html=False, include_plotlyjs=include_plotlyjs)
						elif exprsn_type == "CRISPR_Broad":
							fig_CCLE = plotly_heatmap(df = df, colorbar_title = "Essentiality")
							replace_dict['{{'+exprsn_type+'_'+lineage+'}}'] = fig_CCLE.to_html(full_html=False, include_plotlyjs=include_plotlyjs)
						elif exprsn_type == "Proteomics":
							fig_CCLE = plotly_heatmap_wo_dendrogram(df = df, colorbar_title = "Abundance")
							replace_dict['{{'+exprsn_type+'_'+lineage+'}}'] = fig_CCLE.to_html(full_html=False, include_plotlyjs=include_plotlyjs)
						counter = 1
				
				output_file = os.path.join(output_path, "CCLE_"+exprsn_type+'.html')
				with open(template_file, 'r') as template, open(output_file, 'w') as output:
				# with open(template_file, 'r') as template:
					for line in template:
						for key, val in replace_dict.items():
							line = line.replace(key, val)
						output.write(line)
	if "GTEx" in modules_selected:	
		exprsn_type_list = ['Exprsn']
		# if the corresponding template exists
		for exprsn_type in exprsn_type_list:
			# if template file exists, copy the template file to output folder
			template_file = os.path.join(template_path, "GTEx_"+exprsn_type+'.html')
			if os.path.isfile(template_file):
				output_file = os.path.join(output_path, "GTEx_"+exprsn_type+'.html')
				shutil.copyfile(template_file, output_file)
				








