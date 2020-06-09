
import pandas as pd

tf_list  = ['p53','Gata1','CTCF','Egr1','Elk1','Ets1','Creb','Runx1','Max-Myc','Mad-Max','Max-Max','RelA','GR','Stat3','TBP']
pdblist  = ['1an2','1dux','1hjb','1k79','1nkp','1nlw','1p47','1qne','1r4r','3gat','3kz8','5t00','5u01','5zk1','6qhd']
resolist = [ 2.9  , 2.1  , 3.0  , 2.4  , 1.8  , 2.0  , 2.2  , 1.9  , 3.0  , None , 1.91 , 2.19 , 2.5  , 3.05 , 2.85 ]

df0 = pd.DataFrame()
df0['pdbid'] = pd.Series(pdblist)
df0['reso'] = pd.Series(resolist)
df0 = df0.sort_values(['pdbid'])

df = pd.read_csv("../TF_PDB/Update/Data/refine_pdbinfo_final.csv")
subdf = df.loc[df.name.isin(tf_list)].reset_index()
count = subdf.groupby(['name']).count()['pdbid'].reset_index().rename(columns={'name':'name1','pdbid':'count'})
print(count)
df = df.loc[df.pdbid.isin(pdblist)][['name','pdbid','bseq1']].sort_values(['name']).reset_index()
a = pd.concat([df,count],axis=1)
a = a.sort_values(['pdbid']).reset_index().rename(columns={'pdbid':'pdbid1'})
b = pd.concat([df0,a],axis=1)
c = b[['name','count','pdbid','reso','bseq1']]
c.to_csv("./Data/final_pdbid.csv",index=False)

