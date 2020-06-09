
import sys
import pandas as pd
import numpy as np
import pylab as pl

delta = 1e-6

df_crystal = pd.read_csv("../../../StemTable_Crystal_5bp.csv",keep_default_na=False, na_values=[''])
df_nmr = pd.read_csv("../../../StemTable_Nmr_5bp.csv",keep_default_na=False, na_values=[''])

bp_params = ['shear','stretch','stagger','buckle','propeller','opening','dist_c1pc1p']

bdna_crystal = df_crystal.loc[(df_crystal['mmtype'].isin(['Protein#DNA','Protein#DNA#RNA','Protein#DNA#DNA/RNA Hybrid'])) & (df_crystal['reso'] <= 3.0) 
    & (df_crystal['ligand'].isnull()) & (df_crystal['sform'] == 'B') 
    & (df_crystal['hf'] != '....')
    & (df_crystal['mismatch'] == False) & (df_crystal['modify'] == False)
    & (df_crystal['ct_bpname'].isin(['AT','CG']))]

bdna_nmr = df_nmr.loc[(df_nmr['mmtype'].isin(['Protein#DNA','Protein#DNA#RNA','Protein#DNA#DNA/RNA Hybrid'])) 
    & (df_nmr['ligand'].isnull()) & (df_nmr['sform'] == 'B') 
    & (df_crystal['hf'] != '....')
    & (df_nmr['mismatch'] == False) & (df_nmr['modify'] == False)
    & (df_crystal['ct_bpname'].isin(['AT','CG']))]

# invert the sign of shear and buckle if pyrimidine comes first
s = bdna_crystal.ct_bpmotif.isin(['A-T','G-C']).map({True:1,False:-1})
bdna_crystal['shear'] = bdna_crystal['shear'].mul(s,axis=0)
bdna_crystal['buckle'] = bdna_crystal['buckle'].mul(s,axis=0)

s = bdna_nmr.ct_bpmotif.isin(['A-T','G-C']).map({True:1,False:-1})
bdna_nmr['shear'] = bdna_nmr['shear'].mul(s,axis=0)
bdna_nmr['buckle'] = bdna_nmr['buckle'].mul(s,axis=0)


# separate at and gc dataset
at_crystal = bdna_crystal.loc[bdna_crystal['ct_bpname'] == 'AT']
gc_crystal = bdna_crystal.loc[bdna_crystal['ct_bpname'] == 'CG']
cano_crystal = at_crystal.append(gc_crystal)

at_nmr = bdna_nmr.loc[bdna_nmr['ct_bpname'] == 'AT']
gc_nmr = bdna_nmr.loc[bdna_nmr['ct_bpname'] == 'CG']
cano_nmr = at_nmr.append(gc_nmr)

at_all = at_crystal.append(at_nmr)
gc_all = gc_crystal.append(gc_nmr)
cano_all = at_all.append(gc_all)

# include TF or not
df_tf = pd.read_table('../../../TF_name_with_01_annotations.txt',sep='\t',usecols=(0,1))
tf_list = df_tf[df_tf.isTF == 1]['pdbid'].unique()
cano_all = cano_all.loc[cano_all.pdbid.isin(tf_list)]

print len(cano_all)
print len(cano_all.pdbid.unique())

free = pd.read_csv("../free/Canonical_all_bp_params.csv")

free_env_1sd = {}
free_env_2sd = {}
free_env_3sd = {}
free_env_bound = {}
for param in free.param.tolist():
    mean = free['param_mean'][free.param == param].values[0]
    sigma = free['param_sigma'][free.param == param].values[0]
    max = free['param_max'][free.param == param].values[0]
    min = free['param_min'][free.param == param].values[0]
    free_env_1sd[param] = [mean - sigma - delta, mean + sigma + delta]
    free_env_2sd[param] = [mean - 2*sigma - delta, mean + 2*sigma + delta]
    free_env_3sd[param] = [mean - 3*sigma - delta, mean + 3*sigma + delta]
    free_env_bound[param] = [min-delta,max+delta]

print free_env_1sd
print free_env_2sd
print free_env_3sd
print free_env_bound


outside_1sd = cano_all.loc[(~cano_all.shear.between(free_env_1sd['shear'][0],free_env_1sd['shear'][1])) | (~cano_all.stretch.between(free_env_1sd['stretch'][0],free_env_1sd['stretch'][1])) | (~cano_all.stagger.between(free_env_1sd['stagger'][0],free_env_1sd['stagger'][1])) | (~cano_all.buckle.between(free_env_1sd['buckle'][0],free_env_1sd['buckle'][1])) | (~cano_all.propeller.between(free_env_1sd['propeller'][0],free_env_1sd['propeller'][1])) | (~cano_all.opening.between(free_env_1sd['opening'][0],free_env_1sd['opening'][1])) | (~cano_all.dist_c1pc1p.between(free_env_1sd['dist_c1pc1p'][0],free_env_1sd['dist_c1pc1p'][1]))]

outside_2sd = cano_all.loc[(~cano_all.shear.between(free_env_2sd['shear'][0],free_env_2sd['shear'][1])) | (~cano_all.stretch.between(free_env_2sd['stretch'][0],free_env_2sd['stretch'][1])) | (~cano_all.stagger.between(free_env_2sd['stagger'][0],free_env_2sd['stagger'][1])) | (~cano_all.buckle.between(free_env_2sd['buckle'][0],free_env_2sd['buckle'][1])) | (~cano_all.propeller.between(free_env_2sd['propeller'][0],free_env_2sd['propeller'][1])) | (~cano_all.opening.between(free_env_2sd['opening'][0],free_env_2sd['opening'][1])) | (~cano_all.dist_c1pc1p.between(free_env_2sd['dist_c1pc1p'][0],free_env_2sd['dist_c1pc1p'][1]))]

outside_3sd = cano_all.loc[(~cano_all.shear.between(free_env_3sd['shear'][0],free_env_3sd['shear'][1])) | (~cano_all.stretch.between(free_env_3sd['stretch'][0],free_env_3sd['stretch'][1])) | (~cano_all.stagger.between(free_env_3sd['stagger'][0],free_env_3sd['stagger'][1])) | (~cano_all.buckle.between(free_env_3sd['buckle'][0],free_env_3sd['buckle'][1])) | (~cano_all.propeller.between(free_env_3sd['propeller'][0],free_env_3sd['propeller'][1])) | (~cano_all.opening.between(free_env_3sd['opening'][0],free_env_3sd['opening'][1])) | (~cano_all.dist_c1pc1p.between(free_env_3sd['dist_c1pc1p'][0],free_env_3sd['dist_c1pc1p'][1]))]

outside_bound = cano_all.loc[(~cano_all.shear.between(free_env_bound['shear'][0],free_env_bound['shear'][1])) | (~cano_all.stretch.between(free_env_bound['stretch'][0],free_env_bound['stretch'][1])) | (~cano_all.stagger.between(free_env_bound['stagger'][0],free_env_bound['stagger'][1])) | (~cano_all.buckle.between(free_env_bound['buckle'][0],free_env_bound['buckle'][1])) | (~cano_all.propeller.between(free_env_bound['propeller'][0],free_env_bound['propeller'][1])) | (~cano_all.opening.between(free_env_bound['opening'][0],free_env_bound['opening'][1])) | (~cano_all.dist_c1pc1p.between(free_env_bound['dist_c1pc1p'][0],free_env_bound['dist_c1pc1p'][1]))]

print("outside")
print len(outside_1sd)
print len(outside_1sd.pdbid.unique())
print len(outside_2sd)
print len(outside_2sd.pdbid.unique())
print len(outside_3sd)
print len(outside_3sd.pdbid.unique())
print len(outside_bound)
print len(outside_bound.pdbid.unique())
