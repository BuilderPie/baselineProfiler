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
from baseP.prepare_data.Annotation import matchSymbols
from baseP.prepare_data.others.Neuro2a import Neuro2a_preprocess
from baseP.prepare_data.CCLE import CCLE_preprocess
from baseP.prepare_data.GTEx import GTEx_preprocess
from baseP.prepare_data.HPA import HPA_preprocess



def symbols_precompute():
    # print('CCLE_Data[Datasets] (Before CCLE_precompute):', CCLE_Data['Datasets'])
    # exec_path= os.path.join(os.path.dirname(os.path.abspath(__file__)), 'otherPreprocessing/CCLE')
    matchSymbols.analyze()

def Neuro2a_precompute():
    Neuro2a_preprocess.analyze()

def CCLE_precompute():
    CCLE_preprocess.analyze()

def GTEx_precompute():
    GTEx_preprocess.analyze()

def HPA_precompute():
    HPA_preprocess.analyze()
