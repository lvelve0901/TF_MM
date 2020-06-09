
#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
import pylab as pl

itemlist = ["Free\nDNA","Protein\nDNA","TF\nDNA"]


# load TF annotation and create a list of TF-DNA PDBID
df_tf = pd.read_table('../../../Survey/TF_name_with_01_annotations.txt',sep='\t',usecols=(0,1))
tf_list = df_tf[df_tf.isTF == 1]['pdbid'].unique()

df_free = pd.read_csv('../../../Survey/Survey_goldWC/stem/free/Golden_WC_free.csv').rename(index=str,columns={"dist_c1pc1p":"C1C1_dist","shift_3p":"shift","slide_3p":"slide","rise_3p":"rise","tilt_3p":"tilt","roll_3p":"roll","twist_3p":"twist","minorgw_3p":"minorgw","majorgw_3p":"majorgw"})
df_bound = pd.read_csv('../../../Survey/Survey_goldWC/stem/bound/Golden_WC_bound.csv').rename(index=str,columns={"dist_c1pc1p":"C1C1_dist","shift_3p":"shift","slide_3p":"slide","rise_3p":"rise","tilt_3p":"tilt","roll_3p":"roll","twist_3p":"twist","minorgw_3p":"minorgw","majorgw_3p":"majorgw"})
df_tf = df_bound.loc[df_bound.pdbid.isin(tf_list)]
df_hg = pd.read_csv('../../../Survey/Survey_HG/crystal/Golden_HG_bound_crystal.csv')

dflist = []

dflist.append(df_free)
dflist.append(df_bound)
dflist.append(df_tf)

bin_mini = 15

plot_list = ['shift','slide','rise',
              'tilt','roll','twist',
              'minorgw','majorgw']

range_list = [[-5,5],[-5,5],[1,8],
              [-40,40],[-80,80],[-40,80],
              [-2,20],[-5,20]]

fig = pl.figure(1,figsize=(8,10))

for idx, param in enumerate(plot_list):

    ax = fig.add_subplot(np.ceil(len(plot_list)/3.),3,idx+1)
    ax.patch.set_facecolor('none')

    # append all the data to ys according to dflist (combine AT, GC, TA, CG to be WC)
    free = dflist[0][param].values
    bound = dflist[1][param].values
    tf = dflist[2][param].values

    free = free[~np.isnan(free)]
    bound = bound[~np.isnan(bound)]
    tf = tf[~np.isnan(tf)]

    if param == "majorgw" or param == "minorgw":
        free = free - 5.8
        bound = bound - 5.8
        tf = tf - 5.8

    parts = ax.violinplot([free,bound,tf],positions=range(len(dflist)),
                  showmedians=True, showextrema=False)

    for i, pc in enumerate(parts['bodies']):
        
        if i == 0:
            pc.set_facecolor('rosybrown')
        elif i == 1:
            pc.set_facecolor('slateblue')
        elif i == 2:
            pc.set_facecolor('limegreen')
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    parts['cmedians'].set_linewidth(4)
    parts['cmedians'].set_color('black')

    ax.set_xlim((-1+0.5,len(dflist)-0.5))
    ax.set_ylim(range_list[idx][0],range_list[idx][1])
    ax.set_xticks(range(len(dflist)))
    ax.set_xticklabels(itemlist,fontsize=12)
    ax.set_ylabel(param,fontsize=20)
    
    free_upper = np.max(free)
    free_lower = np.min(free)
    bound_upper = np.max(bound)
    bound_lower = np.min(bound)
    tf_upper = np.max(tf)
    tf_lower = np.min(tf)

    ax.hlines(free_upper,-1,len(dflist),colors='rosybrown',linestyles='dashed',linewidth=3)
    ax.hlines(free_lower,-1,len(dflist),colors='rosybrown',linestyles='dashed',linewidth=3)
    ax.hlines(bound_upper,-1,len(dflist),colors='slateblue',linestyles='dashed',linewidth=3)
    ax.hlines(bound_lower,-1,len(dflist),colors='slateblue',linestyles='dashed',linewidth=3)
    ax.hlines(tf_upper,-1,len(dflist),colors='limegreen',linestyles='dashed',linewidth=3)
    ax.hlines(tf_lower,-1,len(dflist),colors='limegreen',linestyles='dashed',linewidth=3)


fig.tight_layout()
pl.savefig("./Plot/Survey_step_params_violin.pdf")
pl.show()

