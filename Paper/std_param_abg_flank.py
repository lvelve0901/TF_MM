
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
    dflist.append(pd.read_table("../MM_free/ABG/myc_%s_free/abglist.txt"%mismatch,delim_whitespace=True))


bin_mini = 15

plot_list = ['beta','gamma','zeta']
range_list = [[0,90],[-200,200],[-60,60]]

std_df = pd.DataFrame([],columns=['mismatch'])
std_df['mismatch'] = pd.Series(itemlist[1:])

for sigma in [0.5,1.0]:

    for resi in range(7,10):
        flank = resi - 8
    
        for idx, param in enumerate(plot_list):
        
            print("Process the parameters: %s"%param)
            
            # append all the data to ys according to itemlist (combine AT, GC, TA, CG to be WC)
            ys = []
            for jdx, mismatch in enumerate(mismatches):
                df = dflist[jdx]
                ys.append(df.loc[(df.RMSD1 <= 2.0) & (df.RMSD2 <= 2.0) & (df.resi == resi)][param].values)
            ys = np.array([y[~np.isnan(y)] for y in ys])
        
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
        std_df.to_csv("./Data/cutoff_mm_value_abg_sigma%.1f_flank_%s.csv"%(sigma,flank),index=False,float_format='%.0f')
        
        std_df = std_df.replace(1.0,'Pos')
        std_df = std_df.replace(-1.0,'Neg')
    
        std_df.to_csv("./Data/cutoff_mm_string_abg_sigma%.1f_flank_%s.csv"%(sigma,flank),index=False)
