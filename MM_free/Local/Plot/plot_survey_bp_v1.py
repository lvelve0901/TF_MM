
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

df_free = pd.read_csv('../../../Survey/Survey_goldWC/stem/free/Golden_WC_free.csv').rename(index=str,columns={"dist_c1pc1p":"C1C1_dist"})
df_bound = pd.read_csv('../../../Survey/Survey_goldWC/stem/bound/Golden_WC_bound.csv').rename(index=str,columns={"dist_c1pc1p":"C1C1_dist"})
df_tf = df_bound.loc[df_bound.pdbid.isin(tf_list)]
#df_hg = pd.read_csv('../../../Survey/Survey_HG/crystal/Golden_HG_bound_crystal.csv')

dflist = []

dflist.append(df_free)
#dflist.append(df_bound)
dflist.append(df_tf)


plot_list = ['shear','stretch','stagger',
              'buckle','propeller','opening',
              'C1C1_dist']

range_list = [[-8,8],[-4,6],[-3,3],
              [-60,80],[-80,80],[-60,120],
              [6,15]]

fig = pl.figure(1,figsize=(8,10))

for idx, param in enumerate(plot_list):

    ax = fig.add_subplot(np.ceil(len(plot_list)/3.),3,idx+1)
    ax.patch.set_facecolor('none')

    # append all the data to ys according to dflist (combine AT, GC, TA, CG to be WC)
    free = dflist[0][param].values
    #bound = dflist[1][param].values
    #tf = dflist[2][param].values
    tf = dflist[1][param].values

    # PDB bound include HG bp
    #if param == 'C1C1_dist':
    #    bound = np.hstack((bound,df_hg['C1C1_dist'].values))
    #    tf = np.hstack((tf,df_hg['C1C1_dist'].values))

    #parts = ax.violinplot([free,bound,tf],positions=range(len(dflist)),
    parts = ax.violinplot([free,tf],positions=range(len(dflist)),
                  showmedians=True, showextrema=False)

    for i, pc in enumerate(parts['bodies']):
        
        if i == 0:
            pc.set_facecolor('rosybrown')
        elif i == 1:
            #pc.set_facecolor('slateblue')
            pc.set_facecolor('limegreen')
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
    #bound_upper = np.max(bound)
    #bound_lower = np.min(bound)
    tf_upper = np.max(tf)
    tf_lower = np.min(tf)

    ax.hlines(free_upper,-1,len(dflist),colors='rosybrown',linestyles='dashed',linewidth=3)
    ax.hlines(free_lower,-1,len(dflist),colors='rosybrown',linestyles='dashed',linewidth=3)
    #ax.hlines(bound_upper,-1,len(dflist),colors='slateblue',linestyles='dashed',linewidth=3)
    #ax.hlines(bound_lower,-1,len(dflist),colors='slateblue',linestyles='dashed',linewidth=3)
    ax.hlines(tf_upper,-1,len(dflist),colors='limegreen',linestyles='dashed',linewidth=3)
    ax.hlines(tf_lower,-1,len(dflist),colors='limegreen',linestyles='dashed',linewidth=3)


fig.tight_layout()
pl.savefig("./Plot/Survey_bp_params_violin.pdf")
pl.show()

