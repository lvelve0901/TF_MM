
#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
import pylab as pl

# mismatch code of database file
mismatches = ['AT','GC','TA','CG','GT','AC','CT','TT','CC','AA','GA2','GA1','GG2','GG1']
# mismatch code representation xticklabels
itemlist = ['WC','GT','AC','CT','TT','CC','AA','GA','GA*\nsyn','G*G\nsyn','GG*\nsyn']
# load TF annotation and create a list of TF-DNA PDBID
df_tf = pd.read_table('../../../Survey/TF_name_with_01_annotations.txt',sep='\t',usecols=(0,1))
tf_list = df_tf[df_tf.isTF == 1]['pdbid'].unique()

df_free = pd.read_csv('../../../Survey/Survey_goldWC/stem/free/Golden_WC_free.csv').rename(index=str,columns={"dist_c1pc1p":"C1C1_dist","shift_3p":"shift","slide_3p":"slide","rise_3p":"rise","tilt_3p":"tilt","roll_3p":"roll","twist_3p":"twist","minorgw_3p":"minorgw","majorgw_3p":"majorgw"})
df_bound = pd.read_csv('../../../Survey/Survey_goldWC/stem/bound/Golden_WC_bound.csv').rename(index=str,columns={"dist_c1pc1p":"C1C1_dist","shift_3p":"shift","slide_3p":"slide","rise_3p":"rise","tilt_3p":"tilt","roll_3p":"roll","twist_3p":"twist","minorgw_3p":"minorgw","majorgw_3p":"majorgw"})
#df_bound = df_bound.loc[df_bound.pdbid.isin(tf_list)]
df_hg = pd.read_csv('../../../Survey/Survey_HG/crystal/Golden_HG_bound_crystal.csv')

dflist = []
for mismatch in mismatches:
    dflist.append(pd.read_csv("../Data/myc_%s_free_step.csv"%mismatch).rename(index=str,columns={"minorgw1":"minorgw","majorgw1":"majorgw"}))


dflist.append(df_free)
dflist.append(df_bound)

bin_mini = 15

plot_list = ['shift','slide','rise',
              'tilt','roll','twist',
              'minorgw','majorgw']

range_list = [[-5,5],[-5,5],[1,8],
              [-40,40],[-80,80],[-40,80],
              [-2,20],[-5,20]]

fig = pl.figure(1,figsize=(25,10))

for idx, param in enumerate(plot_list):

    ax = fig.add_subplot(np.ceil(len(plot_list)/3.),3,idx+1)
    ax.patch.set_facecolor('none')

    # append all the data to ys according to itemlist (combine AT, GC, TA, CG to be WC)
    free = dflist[-2][param].values
    bound = dflist[-1][param].values
    if param == 'minorgw' or param == 'majorgw':
        free = free[~np.isnan(free)]
        free = free - 5.8
        bound = bound[~np.isnan(bound)]
        bound = bound - 5.8

    ys = []
    for jdx, item in enumerate(mismatches):
        y = dflist[jdx][param].loc[dflist[jdx]['step'] == 2].values
        y = y[~np.isnan(y)]
        if param == 'minorgw' or param == 'majorgw':
            y = y - 5.8
        ys.append(y)

    control = np.hstack((ys[0],ys[1],ys[2],ys[3]))

    parts = ax.violinplot([control]+ys[4:len(ys)],positions=range(len(ys)-4+1),
                  showmedians=True, showextrema=False)

    for i, pc in enumerate(parts['bodies']):
        
        pc.set_facecolor('orange')
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    parts['cmedians'].set_linewidth(3)
    parts['cmedians'].set_color('black')

    ax.set_xlim((-1,len(itemlist)))
    ax.set_ylim(range_list[idx][0],range_list[idx][1])
    ax.set_xticks(range(len(ys)-4+1))
    ax.set_xticklabels(itemlist,fontsize=16)
    ax.set_ylabel(param,fontsize=20)
    
    free_upper = np.max(free)
    free_lower = np.min(free)
    bound_upper = np.max(bound)
    bound_lower = np.min(bound)

    #ax.hlines(free_upper,-1,len(itemlist),colors='rosybrown',linestyles='dashed',linewidth=3)
    #ax.hlines(free_lower,-1,len(itemlist),colors='rosybrown',linestyles='dashed',linewidth=3)
    #ax.hlines(bound_upper,-1,len(itemlist),colors='slateblue',linestyles='dashed',linewidth=3)
    #ax.hlines(bound_lower,-1,len(itemlist),colors='slateblue',linestyles='dashed',linewidth=3)

    if param != "minorgw" and param != "majorgw":
        ax.axvspan(7.5,10.5,alpha=0.5,color='gray',edgecolor='None')

fig.tight_layout()
pl.savefig("./Plot/MM_MD_step_params_violin.pdf")
pl.show()

