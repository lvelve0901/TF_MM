
import sys
import pandas as pd
import numpy as np
import pylab as pl

df_crystal = pd.read_csv("../../StemTable_Crystal_5bp.csv",keep_default_na=False, na_values=[''])
df_nmr = pd.read_csv("../../StemTable_Nmr_5bp.csv",keep_default_na=False, na_values=[''])

abg_params = ['beta_h_b','gamma_h_b','zeta_h_b','rmsd_pu_b','rmsd_py_b']

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
df_tf = pd.read_table('../../TF_name_with_01_annotations.txt',sep='\t',usecols=(0,1))
tf_list = df_tf[df_tf.isTF == 1]['pdbid'].unique()
cano_all = cano_all.loc[cano_all.pdbid.isin(tf_list)]

print len(cano_all.pdbid.unique())

table_list = [at_crystal,gc_crystal,cano_crystal,
              at_nmr,gc_nmr,cano_nmr,
              at_all,gc_all,cano_all]

title_list = ["Bend_AT_crystal","Bend_GC_crystal","Bend_Canonical_crystal",
              "Bend_AT_nmr","Bend_GC_nmr","Bend_Canonical_nmr",
              "Bend_AT_all","Bend_GC_all","Bend_Canonical_all"]

for i in range(len(table_list)):
    
    print title_list[i]
    print len(table_list[i])
    values = [['param','param_mean','param_max','param_min','param_sigma']]
    
    for param in abg_params:
        a = table_list[i][param][(table_list[i]['rmsd_pu_b'] < 2.0) & (table_list[i]['rmsd_py_b'] < 2.0)]
        param_mean = round(np.mean(a),2)
        param_max = round(np.max(a),2)
        param_min = round(np.min(a),2)
        param_sigma = round(np.std(a),2)
        values.append([param,param_mean,param_max,param_min,param_sigma])

    df = pd.DataFrame(values[1:],columns=values[0])
    df.to_csv('%s_abg_params.csv'%title_list[i],index=False)

cano_all.to_csv("Bend_WC_bound.csv",index=False)

