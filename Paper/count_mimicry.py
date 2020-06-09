
import os
import sys
import pandas as pd
import numpy as np



input_file = "final_local_and_global_mimicry_sigma1.0.csv"
df = pd.read_csv(os.path.join("Data",input_file))
print(df)
df['mismatch_reduce'] = df['mismatch'].str.replace("*","")
df['mismatch_reduce'] = df['mismatch_reduce'].str[0:2]
df['Key'] = df['pdbid'] + "_" + df['resi'].astype(str) + "_" + df['mismatch_reduce']
df1 = df[['pdbid','resi','mismatch','mismatch_reduce','Key']]


mmdic = {
        '1an2':[
                [9,['GT']],
                #[12,['GT']],
                ],
        '1dux':[
                [4,['CT']],
                [6,['GA','GA*syn']],
                ],
        '1hjb':[
                #[14,['TT','CT']],
                [17,['AC','AA']],
                ],
        '1k79':[
                [5,['CT']],
                [6,['GG*syn','G*Gsyn']],
                [7,['GA','GA*syn','GT']],
                ],
        '1nkp':[
                [7,['GT']],
                #[10,['GT']],
                ],
        #'1nlw':[
        #        [7,['GT']],
        #        [10,['GT']],
        #        ],
        '1p47':[
                [1,['CT']],
                [9,['GG*syn','G*Gsyn']],
                ],
        '1qne':[
                #[2,['CC']],
                #[3,['CC']],
                [5,['TT']],
                [8,['TT']],
                [9,['CC','CT','AC']],
                ],
        '1r4r':[
                [3,['TT','GT']],
                [13,['GA','GA*syn']],
                ],
        #'3gat':[
        #        [5,['CT','AC']],
        #        ],
        '3kz8':[
                [4,['TT','CT']],
                #[5,['TT','CT']],
                [14,['TT','CT']],
                #[15,['TT','CT']],
                ],
        '5t00':[
                [2,['CC']],
                [5,['GT']],
                [6,['GG*syn','G*Gsyn']],
                ],
        '5u01':[
                #[2,['CT','AC','CC']],
                #[3,['GG*syn','G*Gsyn']],
                [22,['GG*syn','G*Gsyn']],
                [23,['CT','AC','CC']],
                ],
        '5zk1':[
                #[4,['CT','TT']],
                [7,['CC']],
                #[12,['CC']],
                #[15,['CT','TT']],
                ],
        #'6qhd':[
        #        [8,['AC']],
        #        [10,['GA']],
        #        ],
        }

pdblist = mmdic.keys()

output = [['pdbid','resi','mismatch']]
for pdbid in pdblist:
    sublist = mmdic[pdbid]
    for sub in sublist:
        for mm in sub[1]:
            output.append([pdbid,sub[0],mm])

df0 = pd.DataFrame(output[1:],columns=output[0])
df0['mismatch_reduce'] = df0['mismatch'].str.replace("*","")
df0['mismatch_reduce'] = df0['mismatch_reduce'].str[0:2]
df0['Key'] = df0['pdbid'] + "_" + df0['resi'].astype(str) + "_" + df0['mismatch_reduce']


print(df1)
print(df0)
explain = list(set(df0.Key.unique()).intersection(set(df1.Key.unique())))
unexplain = list(set(df0.Key.unique())-set(df1.Key.unique()))
print(len(explain))
print(len(unexplain))
print(len(df0.Key.unique()))
print(unexplain)

