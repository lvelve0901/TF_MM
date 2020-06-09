
import pylab as pl
import numpy as np
import pandas as pd


output_df = pd.read_csv("./Data/tf_DNAproDB.csv")

# remove vdw and hbond redundancy
#output_df = output_df.sort_values(['pdbid','int_type','nt_chain','nt_resi','aa_chain','aa_resi']).drop_duplicates(['pdbid','res_pair','nt_moiety'],keep='first')
output_df = output_df.reset_index(drop=True)

#a = output_df.loc[output_df.int_type == 'hbond']['distance'].values
#b = output_df.loc[output_df.int_type == 'vdw']['distance'].values

#print(a)
#print(b)

#pl.hist(a,bins=30,alpha=0.5)
#pl.hist(b,bins=30,alpha=0.5)
#pl.show()

#print(output_df)
l = output_df[output_df.pdbid == '1qne']
print(sorted(l.aa_resi.unique()))
#l = l.sort_values(['aa_resi'])
#l.to_csv("test_1qne_remo.csv",index=False,float_format='%.3f')
##print(len(l))
#l = l.groupby(['int_type','nt_resi','aa_resi','nt_moiety']).count()['distance'].reset_index()
#print(l)
#l.to_csv("test_1qne_remo_1qne.csv",float_format='%.3f')
#print(len(l))
#

m = output_df[output_df.pdbid == 'xac1']
m = m.loc[(m.nt_chain.isin(['B','C'])) & (m.aa_chain.isin(['A']))]
print(sorted(m.aa_resi.unique()))
#m = m.sort_values(['aa_resi'])
#m.to_csv("test_xac1_remo.csv",index=False,float_format='%.3f')
##print(len(m))
#m = m.groupby(['int_type','nt_resi','aa_resi','nt_moiety']).count()['distance'].reset_index()
#m.to_csv("test_1qne_remo_xac1.csv",float_format='%.3f')
#print(len(m))
#


n = output_df[output_df.pdbid == 'xcc2']
print(sorted(n.aa_resi.unique()))
#n = n.sort_values(['aa_resi'])
#n.to_csv("test_xcc2_remo.csv",index=False,float_format='%.3f')
##print(len(n))
#n = n.groupby(['int_type','nt_resi','aa_resi','nt_moiety']).count()['distance'].reset_index()
#n.to_csv("test_1qne_remo_xcc2.csv",float_format='%.3f')
#print(len(n))
#

#l['autokey'] = l.apply(lambda x: x.int_type + "_" + str(x.nt_resi) + '_' + str(x.aa_resi) + '_' + x.nt_moiety,axis=1)
#m['autokey'] = m.apply(lambda x: x.int_type + "_" + str(x.nt_resi) + '_' + str(x.aa_resi) + '_' + x.nt_moiety,axis=1)
#n['autokey'] = n.apply(lambda x: x.int_type + "_" + str(x.nt_resi) + '_' + str(x.aa_resi) + '_' + x.nt_moiety,axis=1)
#
#print(l)

#a = l.merge(n,how='outer',left_on='autokey',right_on='autokey',suffixes=('_l','_r'))
#print(a[['autokey','distance_l','distance_r']])
#a.to_csv('test.csv')

#a = pd.concat([l,m,n],keys='autokey',join='outer',axis=1)
#print(a)
#a.to_csv("test.csv")

l = output_df[output_df.pdbid == '6njq']
l = l.loc[(l.nt_chain.isin(['C','D'])) & (l.aa_chain.isin(['A']))]
print(sorted(l.aa_resi.unique()))
#l = l.groupby(['int_type','res_pair','nt_moiety']).count()
#print(len(l))
##
l = output_df[output_df.pdbid == 'xcca']
l = l.loc[(l.nt_chain.isin(['C','D'])) & (l.aa_chain.isin(['A']))]
print(sorted(l.aa_resi.unique()))
###print(len(l))
#l = l.groupby(['int_type','res_pair','nt_moiety']).count()
#print(len(l))
##
l = output_df[output_df.pdbid == 'xccb']
l = l.loc[(l.nt_chain.isin(['C','D'])) & (l.aa_chain.isin(['A']))]
print(sorted(l.aa_resi.unique()))
###print(len(l))
#l = l.groupby(['int_type','res_pair','nt_moiety']).count()
#print(len(l))


#l = output_df[output_df.pdbid == '1qne']
#l.to_csv("test_1qne_remo.csv",index=False,float_format='%.3f')

#l = output_df[((output_df.pdbid == '3kz8') & (output_df.int_type == 'hbond'))]
#print(l)

