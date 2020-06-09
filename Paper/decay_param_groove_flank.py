
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
colorlist = ['gray','green','orange','cyan','blue','skyblue','maroon','red','purple','yellow','gold']

dflist = []
for mismatch in mismatches:
    dflist.append(pd.read_csv("../MM_free/Local/Data/myc_%s_free_all_step.csv"%mismatch).rename(index=str,columns={"minorgw1":"minorgw","majorgw1":"majorgw"}))


bin_mini = 15

plot_list = ['minorgw','majorgw']

output = [['mismatch']]
std_df = pd.DataFrame([],columns=output[0])
std_df['mismatch'] = pd.Series(itemlist[1:])


for idx, param in enumerate(plot_list):
    print("Process the parameters: %s"%param)

    steps = range(0,10)
    max_std = 0.
    for step in steps:
        flank = step - 5
        ys = []
        
        # append all the data to ys according to itemlist (combine AT, GC, TA, CG to be WC)
        for jdx, mismatch in enumerate(mismatches):
            ys.append(dflist[jdx][param].loc[dflist[jdx]['step'] == step].values)
    
        if param == 'minorgw' or param == 'majorgw':
            ys = np.array([y[~np.isnan(y)] for y in ys])
            ys = ys - 5.8
    
        # change the scale to be -180 to 180
        control = np.hstack((ys[0],ys[1],ys[2],ys[3]))
        
        ref_mean = np.mean(control)
        ref_std = np.std(control)
        if flank == 0:
            max_std = ref_std
        
        results = []
        # calculate REsemble score for mismatch
        for i in range(4,len(ys)):
            mis_mean = np.mean(ys[i])
            delta = mis_mean - ref_mean
            results.append(delta)
    
        std_df["step_"+str(flank)] = pd.Series(results)
    
    print(std_df)
    std_df.to_csv("./Data/decay_mismatch_%s_raw_value.csv"%param,index=False,float_format='%.3f')

    fig = pl.figure(1,figsize=(10,5))
    ax = fig.add_subplot(111)
    ax.axhspan(1.0*max_std,-1.0*max_std,color='gray',alpha=0.2)
    ax.axhspan(0.5*max_std,-0.5*max_std,color='gray',alpha=0.2)
    ax.axvspan(-0.1-0.5,0.1-0.5,color='red',alpha=0.5)
    ax.plot(range(-10,10),[0]*len(range(-10,10)),'-',color='black',linewidth=6)
    for i, mm in enumerate(std_df.mismatch):
        values = std_df.loc[std_df.mismatch == mm].values[0][1:]
        ax.plot(np.array(steps)-5,values,'--o',markersize=12,color=colorlist[i+1],label=mm)

    ax.set_xticks(np.array(steps)-5)
    ax.set_xticklabels(np.array(steps)-5)
    ax.set_xlim(-6,5)
    ax.set_ylim(-2.0*max_std,2.0*max_std)
    ax.set_xlabel("position index")
    ax.set_ylabel("delta %s"%param)
    ax.legend(loc=1,fontsize=10,ncol=2)

    pl.savefig("./Plot/decay_mismatch_%s.pdf"%param)
    pl.show()
