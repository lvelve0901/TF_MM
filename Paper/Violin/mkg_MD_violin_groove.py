
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
resi = int(sys.argv[4])
mismatch_list = sys.argv[5].split('/')
color_list = sys.argv[6].split('/')
ylim_list = map(int,sys.argv[7].split(','))


inpdir = '../../TF_MD/Local/'
oupdir = './violin'
df_list = []

for idx, mismatch in enumerate(mismatch_list):
    inpf = tf_name + "_" + mismatch + "_steps.csv"
    df = pd.read_csv(os.path.join(inpdir,tf_id,inpf))
    df_list.append(df)

fig = pl.figure(1)
fig.patch.set_facecolor('none')
ax = fig.add_subplot(111)
ax.patch.set_facecolor('none')

dfs = []
ys = []

for idx, df in enumerate(df_list):

    subdf = df.loc[df.nt1_i == resi].reset_index(drop=True)
    dfs.append(subdf[param])
    ys.append(subdf[param].values)
    print(mismatch_list[idx])
    print(np.nanmean(subdf[param].values))
    print(np.nanstd(subdf[param].values))

ys = [y[~np.isnan(y)]-5.8 for y in ys]

parts = ax.violinplot(ys,positions=np.arange(len(ys))+0.5,
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

fig.tight_layout()
pl.savefig("%s/%s_%s_violin.pdf"%(oupdir,tf_id,param))
pl.clf()


