

import pandas as pd
import numpy as np
import pylab as pl
import os, sys
from seq_gap import *
from matplotlib.font_manager import FontProperties

pname = sys.argv[1]

os.system("mkdir -p ./Plot/%s"%pname)

pdbinfo = pd.read_csv("./Data/refine_pdbinfo_final.csv")
bseqinfo = pd.read_csv("./Data/tf_bseq.csv")
pdbinfo['bseq1_full'] = bseqinfo['bseq1']
pdbinfo['bseq2_full'] = bseqinfo['bseq2']
pdbinfo['bidx1_full'] = bseqinfo['bidx1']
pdbinfo['bidx2_full'] = bseqinfo['bidx2']

pdbinfo = pdbinfo.loc[pdbinfo.name == pname].reset_index(drop=True)

bp_params = ['contact']

df = pd.read_csv('./Data/tf_contact.csv')

# find gap in helices
gappdb = pdbinfo.loc[~(pdbinfo.bseq1 == pdbinfo.bseq1_full)]['pdbid'].tolist()

# generate plot for each bp params
for idx, param in enumerate(bp_params):

    print("--- Working on %s [%d/%d] ---"%(param,idx,len(bp_params)))

    fig = pl.figure(1,figsize=(6,len(pdbinfo)*6.0))
    fig.patch.set_facecolor('none')

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
            pos = find_gap(bidx1,bidx1_full)
            hbonds_bb_1 = fill_gap([np.nan] + subdf['hbonds_bb_s1'].tolist() + [np.nan],pos) 
            hbonds_bb_2 = fill_gap([np.nan] + subdf['hbonds_bb_s2'].tolist() + [np.nan],pos)  
            hbonds_sc_1 = fill_gap([np.nan] + subdf['hbonds_sc_s1'].tolist() + [np.nan],pos) 
            hbonds_sc_2 = fill_gap([np.nan] + subdf['hbonds_sc_s2'].tolist() + [np.nan],pos) 

        else:

            seqrange = np.array(range(len(bseq1)))
            subdf = df.loc[df['pdbid'] == pdbid]
            hbonds_bb_1 = [np.nan] + subdf['hbonds_bb_s1'].tolist() + [np.nan]
            hbonds_bb_2 = [np.nan] + subdf['hbonds_bb_s2'].tolist() + [np.nan] 
            hbonds_sc_1 = [np.nan] + subdf['hbonds_sc_s1'].tolist() + [np.nan] 
            hbonds_sc_2 = [np.nan] + subdf['hbonds_sc_s2'].tolist() + [np.nan] 

        hbonds = [hbonds_sc_1,hbonds_bb_1,hbonds_bb_2,hbonds_sc_2]

        ax1 = fig.add_subplot(len(pdbinfo),1,i+1)
        ax1.patch.set_facecolor('none')
        ax2 = ax1.twiny()

        ax1.set_xticks(seqrange)
        ax1.set_xticklabels(list(bseq1_full),fontsize=15)
        ax1.axhspan(-3.0,-1.5,color='green',alpha=0.2)
        ax1.axhspan(-1.5,0,color='red',alpha=0.2)
        ax1.axhspan(0,1.5,color='red',alpha=0.2)
        ax1.axhspan(1.5,3.0,color='green',alpha=0.2)

        font = FontProperties()
        color = 'yellow'
        font.set_weight('bold')
        alpha_text=1.0
        alpha_box=0.5
        for l, hbond in enumerate(hbonds):
            for m in seqrange:
                if hbond[m] is np.nan:
                    hbond_string = ""
                else:
                    hbond_string = hbond[m]
                aa_list = hbond_string.split('&')
                for n, aa in enumerate(aa_list):
                    pl.text(m,-3.0+1.5*l+0.2*n+0.2,aa,
                        fontsize=5,
                        fontproperties=font,
                        horizontalalignment='center',
                        verticalalignment='center',
                        alpha=alpha_text,
                        bbox=dict(
                            boxstyle='round',
                            facecolor=color,
                            alpha=alpha_box)
                    )
            
        

        ax1.set_ylabel(name+'('+pdbid+')')
        ax1.set_xlim(0,len(bseq1_full)-1)
        ax1.set_ylim(-3.0,3.0)
        ax1.get_yaxis().set_ticks([])
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

    pl.savefig('./Plot/%s/tf_final_%s.pdf'%(pname,param))

    pl.clf()


