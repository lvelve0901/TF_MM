
#!/usr/bin/python

import os
import sys
import numpy as np
import scipy.stats
import pandas as pd
import pylab as pl
from ensembtool.jsd import *


# mismatch code of database file
mismatches = ['AT','GC','TA','CG','GT','AC','CT','TT','CC','AA','GA2','GA1','GG2','GG1']
itemlist = ['WC','GT','AC','CT','TT','CC','AA','GA','GA*syn','G*Gsyn','GG*syn']

dflist = []
for mismatch in mismatches:
    dflist.append(pd.read_csv("../MM_free/Local/Data/myc_%s_free_all_step.csv"%mismatch).rename(index=str,columns={"minorgw1":"minorgw","majorgw1":"majorgw"}))


bin_mini = 15

plot_list = ['minorgw','majorgw']

output = [['mismatch']+plot_list]
std_df = pd.DataFrame([],columns=output[0])
std_df['mismatch'] = pd.Series(itemlist[1:])


for sigma in [0.5,1.0]:

    for step in range(3,7):
        flank = step - 5
        print(flank)
    
        for idx, param in enumerate(plot_list):
        
            print("Process the parameters: %s"%param)
            
            # append all the data to ys according to itemlist (combine AT, GC, TA, CG to be WC)
            ys = []
            for jdx, mismatch in enumerate(mismatches):
                ys.append(dflist[jdx][param].loc[dflist[jdx]['step'] == step].values)
        
            if param == 'minorgw' or param == 'majorgw':
                ys = np.array([y[~np.isnan(y)] for y in ys])
                ys = ys - 5.8
        
            # change the scale to be -180 to 180
            control = np.hstack((ys[0],ys[1],ys[2],ys[3]))
            
            ref_mean = np.mean(control)
            ref_std = sigma*np.std(control)
            
            results = []
            # calculate REsemble score for mismatch
            for i in range(4,len(ys)):
                mis_mean = np.mean(ys[i])
                if mis_mean >= ref_mean + ref_std:
                    result = 1
                elif mis_mean <= ref_mean - ref_std:
                    result = -1
                else:
                    result = 0
                results.append(result)
        
            std_df[param] = pd.Series(results)
        
        print(std_df)
        std_df.to_csv("./Data/cutoff_mm_value_groove_sigma%.1f_flank_%s.csv"%(sigma,flank),index=False,float_format='%.0f')
    
        std_df = std_df.replace(1.0,'Pos')
        std_df = std_df.replace(-1.0,'Neg')
    
        std_df.to_csv("./Data/cutoff_mm_string_groove_sigma%.1f_flank_%s.csv"%(sigma,flank),index=False)
