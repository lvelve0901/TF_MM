
import re
import os
import sys
import pandas as pd
import pylab as pl
import numpy as np
from commontool import read
from pprint import pprint

tf_name = sys.argv[1]
tf_id = sys.argv[2]
param = sys.argv[3]
mismatch_list = sys.argv[4].split('/')
color_list = sys.argv[5].split('/')
ylim_list = map(int,sys.argv[6].split(','))


inpdir = '../../TF_MD/Interaction/'
oupdir = './violin'
df_list = []
dfs = []
ys = []

for idx, mismatch in enumerate(mismatch_list):
    
    inpf = param + '_' + mismatch + '_net.txt'
    if param[0:2] == 'sa':
        df = pd.read_table(os.path.join(inpdir,tf_id,inpf),delim_whitespace=True,skiprows=1,names=['id','value'])
        df = df.loc[((df['value'] >= 0) & (df['value'] <= 3000))]
    else:
        df = pd.read_table(os.path.join(inpdir,tf_id,inpf),delim_whitespace=True,names=['id','value'])
        df = df
    dfs.append(df['value'])
    y = df['value'].values
    ys.append(y)

fig = pl.figure(1)
fig.patch.set_facecolor('none')
ax = fig.add_subplot(111)
ax.patch.set_facecolor('none')

if param[0:2] == 'sa':
    parts = ax.violinplot(ys,positions=np.arange(len(ys))+0.5,
                  showmeans=True, showextrema=False)
else:    
    parts = ax.violinplot(ys,positions=np.arange(len(ys))+0.5,
                  bw_method=0.5,
                  showmeans=True, showextrema=False)

for i, pc in enumerate(parts['bodies']):
   
    pc.set_facecolor(color_list[i])
    pc.set_edgecolor('black')
    pc.set_alpha(1.0)

parts['cmeans'].set_linewidth(3)
parts['cmeans'].set_color('black')

ax.set_xlim((0,len(ys)))
ax.set_xticks(np.array(range(len(ys)+1))+0.5)
ax.set_xticklabels(mismatch_list)
ax.set_ylim((ylim_list[0],ylim_list[1]))
ax.set_yticks([ylim_list[0],(ylim_list[0]+ylim_list[1])/2.,ylim_list[1]])

fig.tight_layout()
pl.savefig("%s/%s_%s_violin.pdf"%(oupdir,tf_id,param))
pl.clf()


