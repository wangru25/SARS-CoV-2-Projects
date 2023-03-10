# -*- coding: utf-8 -*-
'''
@Author: Rui Wang
@Date: 2020-07-06 10:38:43
@LastModifiedBy: Rui Wang
LastEditTime: 2021-05-20 10:05:43
@Email: wangru25@msu.edu
FilePath: /38_Influ/analysis/Animation/seperateData.py
@Description: 
'''
# import plotly.graph_objects as go

import numpy as np
np.random.seed(1)
from scipy.signal import savgol_filter
import math
import pandas as pd
import sys

Date = sys.argv[1]
OldDate = sys.argv[2] 

a = pd.read_csv('./snpRecords_%s_new.csv'%Date)
b = pd.read_csv('./snpRecords_%s.csv'%OldDate)

new_csv = pd.concat([b,a],axis = 0)
new_csv.to_csv('./snpRecords_%s.csv'%Date, index=False)

# c = pd.read_csv('./snpRecords_%s.csv'%Date)

for i in range(math.ceil(a.shape[0]/10000)):
    save_data = a.iloc[i*10000:(i+1)*10000]
    save_data.to_csv('./snpRecords_%s_new_%d.csv'%(Date,i),index = False)
