

import pandas as pd
import numpy as np
import pylab as pl
import os, sys
from seq_gap import *


pdbinfo = pd.read_csv("../Data/refine_pdbinfo_final.csv")
bseqinfo = pd.read_csv("../Data/tf_bseq.csv")
pdbinfo['bseq1_full'] = bseqinfo['bseq1']
pdbinfo['bseq2_full'] = bseqinfo['bseq2']
pdbinfo['bidx1_full'] = bseqinfo['bidx1']
pdbinfo['bidx2_full'] = bseqinfo['bidx2']


pdblist = ['1an2','1dux','1hjb','1k79','1nkp','1nlw','1p47','1qne','1r4r','3gat','3kz8','5t00','5u01','5zk1','6qhd']
lblist =   [
            [[8.5,9.5],[11.5,12.5]],
            [[5.5,6.5]],
            [[13.5,14.5],[16.5,17.5]],
            [[4.5,7.5]],
            [[6.5,7.5],[9.5,10.5]],
            [[6.5,7.5],[9.5,10.5]],
            [[8.5,9.5]],
            [[1.5,3.5],[7.5,9.5]],
            [[2.5,3.5]],
            [[4.5,5.5]],
            [[3.5,5.5],[13.5,15.5]],
            [[4.5,6.5]],
            [[1.5,3.5],[21.5,23.5]],
            [[3.5,4.5],[6.5,7.5],[11.5,12.5],[14.5,15.5]],
            [[7.5,8.5],[9.5,10.5]]
            ]


pdbinfo = pdbinfo.loc[pdbinfo.pdbid.isin(pdblist)].reset_index(drop=True)


df = pd.read_csv('../Data/tf_local.csv')
df_range = pd.read_csv("../../../Survey/Survey_goldWC/stem/free/Canonical_all_bp_params.csv")

bp_params = ['shift','slide','rise','tilt','roll','twist']
ylim_list = [[-4,4],[-4,4],[0,8],[-60,60],[-60,60],[-20,60]]


# generate B-DNA envelope
envelopes = []
for i in range(9,15):
    
    ymean = df_range.ix[i]['param_mean']
    std = df_range.ix[i]['param_sigma']
    ymax = df_range.ix[i]['param_max']
    ymin = df_range.ix[i]['param_min']

    ylower1 = ymean - std
    yupper1 = ymean + std
    ylower2 = ymean - 2*std
    yupper2 = ymean + 2*std
    ylower3 = ymean - 3*std
    yupper3 = ymean + 3*std

    envelopes.append([ymean, ymin, ymax, ylower1, yupper1, ylower2, yupper2, ylower3, yupper3])


# find gap in helices
gappdb = pdbinfo.loc[~(pdbinfo.bseq1 == pdbinfo.bseq1_full)]['pdbid'].tolist()


# generate plot for each bp params
for idx, param in enumerate(bp_params):

    print("--- Working on %s [%d/%d] ---"%(param,idx,len(bp_params)))

    fig = pl.figure(1,figsize=(18,9))
    fig.patch.set_facecolor('none')

    ymean, ymin, ymax, ylower1, yupper1, ylower2, yupper2, ylower3, yupper3 = envelopes[idx]

    for i in range(len(pdbinfo)):

        pdbid = pdbinfo['pdbid'].ix[i]
        name = pdbinfo['name'].ix[i]


        bseq1 = pdbinfo['bseq1'].ix[i]
        bseq2 = pdbinfo['bseq2'].ix[i]

        bseq1_full = pdbinfo['bseq1_full'].ix[i]
        bseq2_full = pdbinfo['bseq2_full'].ix[i]

        bidx1 = list(map(int,pdbinfo['bidx1'].ix[i].split('|')))
        bidx2 = list(map(int,pdbinfo['bidx2'].ix[i].split('|')))

        bidx1_full = list(map(int,pdbinfo['bidx1_full'].ix[i].split('|')))
        bidx2_full = list(map(int,pdbinfo['bidx2_full'].ix[i].split('|')))
        
        if pdbid in gappdb:

            seqrange = np.array(range(len(bseq1_full)))
            subdf = df.loc[df['pdbid'] == pdbid]
            param_value = subdf[param].tolist()
            param_value = [np.nan] + param_value + [np.nan]
            pos = find_gap(bidx1,bidx1_full)
            param_value = fill_gap(param_value,pos)

        else:

            seqrange = np.array(range(len(bseq1)))
            subdf = df.loc[df['pdbid'] == pdbid]
            param_value = subdf[param].tolist()
            param_value = [np.nan] + param_value + [np.nan]


        ylim_min = np.floor(np.min([ylower1,np.nanmin(param_value[:-2])]))
        ylim_max = np.ceil(np.max([yupper1,np.nanmax(param_value[:-2])]))
        ylim_mid = (ylim_min + ylim_max) / 2.
        ylim_min = round(ylim_min,1)
        ylim_max = round(ylim_max,1)
        ylim_mid = round(ylim_mid,1)

        ax1 = fig.add_subplot(5,3,i+1)
        ax1.patch.set_facecolor('none')
        ax2 = ax1.twiny()

        ax1.set_xticks(seqrange)
        ax1.set_xticklabels(list(bseq1_full),fontsize=15)
        ax1.axhspan(ylower1,yupper1,color='red',alpha=0.4)
        ax1.axhspan(ylower2,yupper2,color='orange',alpha=0.2)
        ax1.axhspan(ylower3,yupper3,color='yellow',alpha=0.2)

        for lb in lblist[i]:
            ax1.axvspan(lb[0],lb[1],alpha=0.6,color='red',edgecolor='None')

        ax1.plot(seqrange[:-2]+0.5,param_value[:-2],'-o',color='blue',
                 markersize=10,linewidth=3,markeredgecolor='none')

        ax1.set_ylabel(name+'('+pdbid+')')
        ax1.set_xlim(0,len(bseq1_full)-1)
        ax1.set_ylim(ylim_min,ylim_max)
        ax1.set_yticks([ylim_min,ylim_mid,ylim_max])
        ax1.set_yticklabels([str(ylim_min),str(ylim_mid),str(ylim_max)])
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(seqrange)
        ax2.set_xticklabels(list(bseq2_full),fontsize=15)

        ax1.spines['left'].set_linewidth(3)
        ax1.spines['right'].set_linewidth(3)
        ax1.spines['top'].set_linewidth(3)
        ax1.spines['bottom'].set_linewidth(3)

        ax1.tick_params(length=6,width=3)
        ax2.tick_params(length=6,width=3)


    fig.tight_layout()


    pl.savefig('./tf_final_%s.pdf'%(param))

    pl.clf()


