#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 07/07/2022
# Last Modified Date: 07/07/2022
# Last Modified By  : Dian Li
import os
import pandas as pd
import json
import joblib
from baseP.configs.data_configs import CCLE_Data
import subprocess





def txt_2_json(fn,input_path,output_path):
	df = pd.read_csv(os.path.join(input_path, fn), sep = '\t', na_values="")
	df = df.loc[:,~df.columns.str.contains('^Unnamed')]
	col_tmp = df.columns
	# col_tmp = df.columns[~df.columns.str.contains('^Unnamed')]
	df = df.fillna("")
	df = df.dropna(axis=1, how='all')
	df.columns = col_tmp

	for col in df.columns:
		if df[col].dtypes == "float64":
			df_temp = df[col][:]
			df_temp[abs(df[col])<0.001] = df[col][abs(df[col])<0.001].map(lambda x: '%2.1e' % x)
			df_temp[(abs(df[col])>=0.001)&(abs(df[col])<1)] = df[col][(abs(df[col])>=0.001)&(abs(df[col])<1)].map(lambda x: '%2.2g' % x)
			df_temp[abs(df[col])>=1] = df[col][abs(df[col])>=1].map(lambda x: '%2.1f' % x)

			df[col] = df_temp
		if 'HIGH' in col.upper():
			df = df.rename(columns={col: "Sig High"})
		if 'LOW' in col.upper():
			df = df.rename(columns={col: "Sig Low"})
	df_json = json.dumps(json.loads(df.to_json(orient='records')), indent=2)
	df_json = '{ \n' + '  "' + 'data' + '"' + ': '  + df_json + '\n}'
	
	df_json = df_json.replace('_', ' ')
	df_json = df_json.replace(';', '; ')

	with open(os.path.join(output_path, fn), 'w') as outfile:
		outfile.write(df_json)

def txt_2_datatables(fn,input_path,output_path,columnsSelect = None, sep = '\t'):
	df = pd.read_csv(os.path.join(input_path, fn), sep = sep, na_values="")
	df = df.loc[:,~df.columns.str.contains('^Unnamed')]
	col_tmp = df.columns
	# col_tmp = df.columns[~df.columns.str.contains('^Unnamed')]
	df = df.fillna("")
	df = df.dropna(axis=1, how='all')
	df.columns = col_tmp

	for col in df.columns:
		if df[col].dtypes == "float64":
			df_temp = df[col][:]
			df_temp[abs(df[col])<0.001] = df[col][abs(df[col])<0.001].map(lambda x: '%2.1e' % x)
			df_temp[(abs(df[col])>=0.001)&(abs(df[col])<1)] = df[col][(abs(df[col])>=0.001)&(abs(df[col])<1)].map(lambda x: '%2.2g' % x)
			df_temp[abs(df[col])>=1] = df[col][abs(df[col])>=1].map(lambda x: '%2.1f' % x)

			df[col] = df_temp
		if 'HIGH' in col.upper():
			df = df.rename(columns={col: "Sig High"})
		if 'LOW' in col.upper():
			df = df.rename(columns={col: "Sig Low"})
	if columnsSelect is None:
		columnsSelect = df.columns.tolist()
	
	str_datatables = '<tbody> \n'
	for index, row in df.iterrows():
		tmp_row = '<tr>'
		for colname in columnsSelect:
			tmp_row = tmp_row + ' <td>' + str(row[colname]) + '</td> '
		tmp_row = tmp_row + '</tr> \n'
		str_datatables = str_datatables + tmp_row
	
	# print(str_datatables)

	with open(os.path.join(output_path, fn), 'w') as outfile:
		outfile.write(str_datatables)