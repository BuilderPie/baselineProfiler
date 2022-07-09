#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Dian Li Jingxin Fu
# Date              : 07/07/2022
# Last Modified Date: 07/07/2022
# Last Modified By  : Dian Li
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import shutil

from baseP.configs.data_configs import DATA_DIR, DATA_COMPUTE, HOM_MouseHumanSequence
from baseP.prepare_data.prepare_utils import txt_2_json, txt_2_datatables
from baseP.configs.data_configs import GTEx_Data

# def geneset(request,logger_P, name, threads, evaluators, resum, html_modules, weights):


def stat_GTEx_Exprsn():
    ######################### Part 1. preprocess GTEx meta file ###############################
    GTEx_sample_info_file = os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'sample_info.csv')

    GTEx_sample_info = pd.read_csv(GTEx_sample_info_file)

    ##### ======= generate a table showing the distribution of lineage and primary disease
    # disease_count = GTEx['primary_disease'].value_counts().to_frame()
    # disease_count.columns = ['value']
    # disease_count['individual'] = disease_count.index
    # disease_count['group'] = 'Disease'

    lineage_count = GTEx_sample_info['SMTS'].value_counts().to_frame()
    lineage_count.columns = ['value']
    lineage_count['individual'] = lineage_count.index
    lineage_count['group'] = 'Lineage'

    # df_out = pd.concat([disease_count, lineage_count], axis=0)
    os.makedirs(os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'tables'), exist_ok=True)
    df_out = lineage_count
    output_file = os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'tables', 'stats_GTEx_Expression_sample_info.txt')
    df_out.to_csv(output_file, sep = '\t', index=False)

    ###### ======= generate a circular barplot from the count table ======== ######
    os.makedirs(os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'figs'), exist_ok=True)

    # R_RMD = os.path.join('baseP', 'prepare_data', 'R_functions', 'src', 'plot_barplot.R')
    
    # r_cmd = ' '.join(['Rscript', R_RMD, 
    # '--file ', os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'tables', 'stats_GTEx_Expression_sample_info.txt').replace(' ', '\ '), 
    # '--output ', os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'figs', 'stats_GTEx_Expression_sample_info.png').replace(' ', '\ '),
    # ])
    # process = subprocess.Popen(r_cmd, shell=True).wait()
    
    # Figure Size
    fig, ax = plt.subplots(figsize =(10, 12))

    # Horizontal Bar Plot
    ax.barh(df_out['individual'], df_out['value'])

    # Remove axes splines
    for s in ['top', 'bottom', 'left', 'right']:
        ax.spines[s].set_visible(False)

    # Remove x, y Ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # Add padding between axes and labels
    ax.xaxis.set_tick_params(pad = 5)
    ax.yaxis.set_tick_params(pad = 10)

    # Add x, y gridlines
    ax.grid(b = True, color ='grey',
            linestyle ='-.', linewidth = 0.5,
            alpha = 0.2)

    # Show top values
    ax.invert_yaxis()

    # Add annotation to bars
    for i in ax.patches:
        plt.text(i.get_width()+0.2, i.get_y()+0.5,
                str(round((i.get_width()), 2)),
                fontsize = 12, fontweight ='bold',
                color ='grey')

    # Add Plot Title
    ax.set_title('22,951 GTEx Tissue Distribution from 979 Donors',
                loc ='left', )

    # Add Text watermark
    # fig.text(0.9, 0.15, 'watermark', fontsize = 12,
    #         color ='grey', ha ='right', va ='bottom',
    #         alpha = 0.7)

    # Show Plot
    # plt.show()
    output_file = os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'figs', 'stats_GTEx_Expression_sample_info.png').replace(' ', '\ ')
    fig.savefig(output_file, dpi = 300)

    dst = os.path.join(DATA_DIR, 'Annotation', 'html_template', 'figs', 'index_file')
    os.makedirs(dst, exist_ok=True)
    shutil.copyfile(src=os.path.join(DATA_COMPUTE, 'GTEx', 'sample_info', 'figs', 'stats_GTEx_Expression_sample_info.png').replace(' ', '\ '),
    dst=os.path.join(dst, 'stats_GTEx_Expression_sample_info.png'))



def analyze():
    stat_GTEx_Exprsn()
    