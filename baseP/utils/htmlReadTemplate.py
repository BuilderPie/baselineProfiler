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
from os import path
import joblib
import warnings
warnings.simplefilter('ignore')
import numpy as np
import pandas as pd
import json
import subprocess
from baseP.configs.data_configs import HTML_Data, DATA_DIR, DATA_COMPUTE
from baseP.utils.utilsPlotly import plotly_heatmap, plotly_heatmap_wo_dendrogram
import datetime as dt

from baseP import CCLE_show_lineages, GTEx_show_lineages, HPA_show_lineages, others_show_lineages

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
	# shutil.copy(os.path.join(template_path, '.nojekyll'), output_path)
	# shutil.copy(os.path.join(template_path, 'index.html'), output_path)
	publish_date = dt.datetime.strftime(dt.datetime.now(),'%b %-d, %Y')
	
	### ============= update index.html ============= ###
	replace_dict = {}
	replace_dict['{{publish_date}}'] = publish_date
	
	replace_dict['{{htmlwidget_CCLE_Exprsn_table}}'] = open(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables', 'CCLE_Expression_sample_info.txt'),'r').read()
	replace_dict['{{htmlwidget_CCLE_Proteomics_table}}'] = open(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables', 'CCLE_Proteomics_sample_info.txt'),'r').read()
	replace_dict['{{htmlwidget_CCLE_CRISPR_Broad_table}}'] = open(os.path.join(DATA_COMPUTE, 'CCLE', 'sample_info', 'datatables', 'CCLE_CRISPR_Broad_sample_info.txt'),'r').read()
	replace_dict['{{htmlwidget_More_HPA_table}}'] = open(os.path.join(DATA_COMPUTE, 'HPA', 'sample_info', 'datatables', 'HPA_Cell_Line_sample_info.txt'),'r').read()
	replace_dict['{{htmlwidget_More_Neuro2a_table}}'] = open(os.path.join(DATA_COMPUTE, 'others', 'Exprsn_Neuro2a', 'datatables', 'Neuro2a_sample_info.txt'),'r').read()
	
	template_file = os.path.join(template_path, 'index.html')
	if os.path.isfile(template_file):
		output_file = os.path.join(output_path, 'index.html')
		with open(template_file, 'r') as template, open(output_file, 'w') as output:
			# with open(template_file, 'r') as template:
			for line in template:
				for key, val in replace_dict.items():
					line = line.replace(key, val)
				output.write(line)

	### ============= update CCLE_*.html ============= ###
	if "CCLE" in modules_selected:
		exprsn_type_list = ['Exprsn', 'Proteomics', 'CRISPR_Broad']
		# if the corresponding template exists
		for exprsn_type in exprsn_type_list:
			template_file = os.path.join(template_path, "CCLE_"+exprsn_type+'.html')
			if os.path.isfile(template_file):
				# create a replace dictionary that store the plotly html code
				# to be replaced into the template
				replace_dict = {}
				replace_dict['{{publish_date}}'] = publish_date
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
							if df.shape[1]<5:
								fig_CCLE = plotly_heatmap_wo_dendrogram(df = df)
							else:
								fig_CCLE = plotly_heatmap(df = df)
							replace_dict['{{'+exprsn_type+'_'+lineage+'}}'] = fig_CCLE.to_html(full_html=False, include_plotlyjs=include_plotlyjs)
						elif exprsn_type == "CRISPR_Broad":
							if df.shape[1]<5:
								fig_CCLE = plotly_heatmap_wo_dendrogram(df = df, colorbar_title = "Essentiality")
							else:
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
	### ============= update GTEx_*.html ============= ###
	if "GTEx" in modules_selected:	
		exprsn_type_list = ['Exprsn']
		# if the corresponding template exists
		for exprsn_type in exprsn_type_list:
			# if template file exists, copy the template file to output folder
			replace_dict = {}
			replace_dict['{{publish_date}}'] = publish_date

			template_file = os.path.join(template_path, "GTEx_"+exprsn_type+'.html')
			if os.path.isfile(template_file):
				output_file = os.path.join(output_path, "GTEx_"+exprsn_type+'.html')
				with open(template_file, 'r') as template, open(output_file, 'w') as output:
				# with open(template_file, 'r') as template:
					for line in template:
						for key, val in replace_dict.items():
							line = line.replace(key, val)
						output.write(line)
				# shutil.copyfile(template_file, output_file)
	### ============= update More_Exprsn.html ============= ###
	if "HPA" in modules_selected:	
		exprsn_type_list = ['Exprsn']
		# if the corresponding template exists
		for exprsn_type in exprsn_type_list:
			# if template file exists, copy the template file to output folder
			template_file = os.path.join(template_path, "More_"+exprsn_type+'.html')
			if os.path.isfile(template_file):
				replace_dict = {}
				replace_dict['{{publish_date}}'] = publish_date
				counter = 0
				for lineage in HPA_show_lineages[exprsn_type]:
					lineage_csv = os.path.join(analysis_path, "HPA_evaluator", exprsn_type, "tables", lineage+".csv")
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
							if df.shape[1]<5:
								fig_CCLE = plotly_heatmap_wo_dendrogram(df = df, colorbar_title = "log2(nTPM)")
							else:
								fig_CCLE = plotly_heatmap(df = df, colorbar_title = "log2(nTPM)")
							replace_dict['{{'+exprsn_type+'_'+lineage+'}}'] = fig_CCLE.to_html(full_html=False, include_plotlyjs=include_plotlyjs)
						counter = 1

				output_file = os.path.join(output_path, "More_"+exprsn_type+'.html')
				# shutil.copyfile(template_file, output_file)	
				with open(template_file, 'r') as template, open(output_file, 'w') as output:
				# with open(template_file, 'r') as template:
					for line in template:
						for key, val in replace_dict.items():
							line = line.replace(key, val)
						output.write(line)		

	if "others" in modules_selected:	
		exprsn_type_list = ['Exprsn_Neuro2a']
		# if the corresponding template exists
		for exprsn_type in exprsn_type_list:
			# if template file exists, copy the template file to output folder
			template_file = os.path.join(output_path, "More_Exprsn.html") # because HPA module already output
			# this template, we have to use the edited one as the template here
			if os.path.isfile(template_file):
				# rename the original template file because need to create an updated file with the same name
				template_file_new = os.path.join(output_path, "More_Exprsn.html.tmp")
				os.rename(template_file, template_file_new)

				replace_dict = {}
				replace_dict['{{publish_date}}'] = publish_date
				# counter = 0
				
				for lineage in others_show_lineages[exprsn_type]:
					lineage_csv = os.path.join(analysis_path, "others_evaluator", 'Exprsn', "tables", exprsn_type+'_'+lineage+".csv")
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
						include_plotlyjs = False # because it is already included when adding HPA figures

						if exprsn_type == "Exprsn_Neuro2a":
							if df.shape[1]<5:
								fig_CCLE = plotly_heatmap_wo_dendrogram(df = df, colorbar_title = "log2(TPM)")
							else:
								fig_CCLE = plotly_heatmap(df = df, colorbar_title = "log2(TPM)")
							replace_dict['{{'+exprsn_type+'_'+lineage+'}}'] = fig_CCLE.to_html(full_html=False, include_plotlyjs=include_plotlyjs)
						# counter = 1

				output_file = os.path.join(output_path, "More_"+'Exprsn'+'.html')
				# shutil.copyfile(template_file, output_file)	
				with open(template_file_new, 'r') as template, open(output_file, 'w') as output:
				# with open(template_file, 'r') as template:
					for line in template:
						for key, val in replace_dict.items():
							line = line.replace(key, val)
						output.write(line)	
				os.remove(template_file_new)





