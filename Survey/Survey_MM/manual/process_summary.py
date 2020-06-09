
import pandas as pd
import numpy as np
import os, sys

mmlist = os.listdir("./atul")
mmlist.sort(key=lambda x: (x.split('_')[2],x.split('_')[3]))

print mmlist

df0 = pd.read_csv("init_automatic.csv")
df0['method'] = pd.Series([])
dflist = []

for mmcsv in mmlist:
    mm = mmcsv.split('_')[2]
    method = mmcsv.split('_')[3]
    if method == 'crystal':
        method = 'X-ray'
    elif method == 'nmr':
        method = 'NMR'
    df = pd.read_csv("./atul/%s"%mmcsv)
    df['method'] = method
    df = df.loc[df['notes'] != 'Not relevant']
    dflist.append(df)

result = pd.concat(dflist,axis=0)
result = result[df0.columns]

result = result[['bp_name','pdbid','method','mmtype','reso','pairid','notes']]
result['pdbid'] = result['pdbid']
result['pairid'] = result['pairid'].apply(lambda x: x.replace("_^",""))

print result

result.to_csv("final_mm_survey.csv",index=False)
