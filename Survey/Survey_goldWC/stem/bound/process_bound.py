
import sys
import pandas as pd
import numpy as np
import pylab as pl

df_crystal = pd.read_csv("../../../StemTable_Crystal_5bp.csv",keep_default_na=False, na_values=[''])
df_nmr = pd.read_csv("../../../StemTable_Nmr_5bp.csv",keep_default_na=False, na_values=[''])

bp_params = ['shear','stretch','stagger','buckle','propeller','opening','dist_c1pc1p','minorgw','majorgw','shift','slide','rise','tilt','roll','twist']

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

# generate a table of PDBID and its corresponding protein name

pname_df = bdna_crystal.drop_duplicates(subset='pdbid',keep='first',inplace=False)[['pdbid','uniprot']]
pname_df.to_csv("TF_name.csv",index=False)


# invert the sign of shear and buckle if pyrimidine comes first
s = bdna_crystal.ct_bpmotif.isin(['A-T','G-C']).map({True:1,False:-1})
bdna_crystal['shear'] = bdna_crystal['shear'].mul(s,axis=0)
bdna_crystal['buckle'] = bdna_crystal['buckle'].mul(s,axis=0)
bdna_crystal['shift_5p'] = bdna_crystal['shift_5p'].mul(s,axis=0)
bdna_crystal['shift_3p'] = bdna_crystal['shift_3p'].mul(s,axis=0)
bdna_crystal['tilt_5p'] = bdna_crystal['tilt_5p'].mul(s,axis=0)
bdna_crystal['tilt_3p'] = bdna_crystal['tilt_3p'].mul(s,axis=0)
bdna_crystal['shift'] = bdna_crystal['shift_3p']
bdna_crystal['slide'] = bdna_crystal['slide_3p']
bdna_crystal['rise'] = bdna_crystal['rise_3p']
bdna_crystal['tilt'] = bdna_crystal['tilt_3p']
bdna_crystal['roll'] = bdna_crystal['roll_3p']
bdna_crystal['twist'] = bdna_crystal['twist_3p']
bdna_crystal['minorgw'] = bdna_crystal['minorgw_3p']
bdna_crystal['majorgw'] = bdna_crystal['majorgw_3p']

s = bdna_nmr.ct_bpmotif.isin(['A-T','G-C']).map({True:1,False:-1})
bdna_nmr['shear'] = bdna_nmr['shear'].mul(s,axis=0)
bdna_nmr['buckle'] = bdna_nmr['buckle'].mul(s,axis=0)
bdna_nmr['shift_5p'] = bdna_nmr['shift_5p'].mul(s,axis=0)
bdna_nmr['shift_3p'] = bdna_nmr['shift_3p'].mul(s,axis=0)
bdna_nmr['tilt_5p'] = bdna_nmr['tilt_5p'].mul(s,axis=0)
bdna_nmr['tilt_3p'] = bdna_nmr['tilt_3p'].mul(s,axis=0)
bdna_nmr['shift'] = bdna_nmr['shift_3p']
bdna_nmr['slide'] = bdna_nmr['slide_3p']
bdna_nmr['rise'] = bdna_nmr['rise_3p']
bdna_nmr['tilt'] = bdna_nmr['tilt_3p']
bdna_nmr['roll'] = bdna_nmr['roll_3p']
bdna_nmr['twist'] = bdna_nmr['twist_3p']
bdna_nmr['minorgw'] = bdna_nmr['minorgw_3p']
bdna_nmr['majorgw'] = bdna_nmr['majorgw_3p']


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

print len(at_all)
print len(gc_all)
print len(cano_all.pdbid.unique())

table_list = [at_crystal,gc_crystal,cano_crystal,
              at_nmr,gc_nmr,cano_nmr,
              at_all,gc_all,cano_all]

title_list = ["AT_crystal","GC_crystal","Canonical_crystal",
              "AT_nmr","GC_nmr","Canonical_nmr",
              "AT_all","GC_all","Canonical_all"]

for i in range(len(table_list)):
    
    print title_list[i]
    print len(table_list[i])
    bp_params = ['shear','stretch','stagger','buckle','propeller','opening','dist_c1pc1p','minorgw','majorgw','shift','slide','rise','tilt','roll','twist']
    values = [['param','param_mean','param_max','param_min','param_sigma']]
    
    for param in bp_params:
        a = table_list[i][param]
        param_mean = round(np.nanmean(a),3)
        param_max = round(np.nanmax(a),3)
        param_min = round(np.nanmin(a),3)
        param_sigma = round(np.nanstd(a),3)
        values.append([param,param_mean,param_max,param_min,param_sigma])

    df = pd.DataFrame(values[1:],columns=values[0])
    df.to_csv('%s_bp_params.csv'%title_list[i],index=False,float_format="%.3f")

cano_all.to_csv("Golden_WC_bound.csv",index=False)
