import pandas as pd
import numpy as np


df = pd.read_csv("../../PairTable_nmr.csv")

screenlist = ['RNA','Protein#RNA']
canonical = ['DA','DG','DC','DT']

pair = df.loc[(~df.mmtype.isin(screenlist))]

print len(pair)

print pair.bp_name.unique()

# unmodified G-T mismatch

gt = pair.loc[(pair.bp_name == 'GT') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

gt.to_csv("DNA_mismatch_GT_nmr.csv",index=False)

# unmodified A-C mismatch

ac = pair.loc[(pair.bp_name == 'AC') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

ac.to_csv("DNA_mismatch_AC_nmr.csv",index=False)

# unmodified G-G mismatch

gg = pair.loc[(pair.bp_name == 'GG') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

gg.to_csv("DNA_mismatch_GG_nmr.csv",index=False)

# unmodified A-G mismatch

ag = pair.loc[(pair.bp_name == 'AG') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

ag.to_csv("DNA_mismatch_AG_nmr.csv",index=False)

# unmodified A-A mismatch

aa = pair.loc[(pair.bp_name == 'AA') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

aa.to_csv("DNA_mismatch_AA_nmr.csv",index=False)

# unmodified T-T mismatch

tt = pair.loc[(pair.bp_name == 'TT') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

tt.to_csv("DNA_mismatch_TT_nmr.csv",index=False)

# unmodified C-T mismatch

ct = pair.loc[(pair.bp_name == 'CT') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

ct.to_csv("DNA_mismatch_CT_nmr.csv",index=False)

# unmodified C-C mismatch

cc = pair.loc[(pair.bp_name == 'CC') & (pair.resn_1.isin(canonical)) & (pair.resn_2.isin(canonical))]

cc.to_csv("DNA_mismatch_CC_nmr.csv",index=False)

