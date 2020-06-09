
import sys
import numpy as np
import pandas as pd

tf_seqs = {
        'p53':'AGACATGCCCGGGCATGCCTC',
        'CTCF':'GCAGCGCCCTCTACTGGCAGC',
        'Stat3':'AAGTTCCTGGAATTT',
        'Gata1':'CGTCGGATATCCGGT',
        'Runx1':'TCTTATGTGGTTTTTG',
        'Creb':'CCGGTGACGTAAACG',
        'Max-Myc':'TGCCACCTGGTGGCCACGTGCC',
        'Mad-Max':'GTAGCACGCGTAAC',
        'Max-Max':'TAGATCACGTGAATG',
        'RelA':'CCTGGGGAATTTCCGGGA',
        'TBP':'ACTATAAAAAGTTC',
        'Ets1':'CCTGTTCCGGTAT',
        'Elk1':'CCCTGTTCCGGTAT',
        'Egr1':'GTGGCGTGGGCGATA',
        'GR':'CAGAACATGATGTTCTCA',
}


#for sigma in [0.5,1.0]:
for sigma in [1.0,2.0,3.0]:

    df1 = pd.read_csv("../TF_PDB/Update/Subset/distortion_all_bp_sigma%.1f.csv"%sigma)
    df2 = pd.read_csv("../TF_PDB/Update/Subset/distortion_all_step_sigma%.1f.csv"%sigma)
    df1['flank'] = np.nan
    df2['flank'] = np.nan
    df3 = pd.read_csv("../TF_PDB/Update/Subset/distortion_all_groove_flank_sigma%.1f.csv"%sigma)
    df4 = pd.read_csv("../TF_PDB/Update/Subset/distortion_all_abg_flank_sigma%.1f.csv"%sigma)
    
    df = pd.concat([df1,df2,df3,df4],axis=0)
    df = df.sort_values(['tf','pdbid','resi']).reset_index(drop=True)
    
    df = df[['tf','flank','pdbid','pairid','isHbond','isVdw','param','value']]
    print(df)
    df['param'] = df['param'].apply(lambda x: x.split('_')[0] if x in ['beta_h_b','gamma_h_b','zeta_h_b'] else x)
    df.to_csv("./Data/final_local_and_global_distortion_all_value_sigma%.1f.csv"%sigma,index=False,float_format='%.3f')
    
    df['param_flank'] = df['param'] + "(" + df['flank'].map(str) + ")"
    df['param_flank'] = df.apply(lambda x: x['param_flank'].split("(")[0] if np.isnan(x['flank']) else x['param_flank'], axis=1)
    
    b = df.groupby(['tf','pdbid','pairid','isHbond','isVdw'])['param_flank'].agg(lambda x: ' '.join(list(x)))
    c = df.groupby(['tf','pdbid','pairid','isHbond','isVdw'])['param_flank'].count()
    
    b = b.reset_index()
    c = c.reset_index()
    
    b = b.rename(columns={'param_flank':'param'})
    b['count'] = c['param_flank']
    
    print(b[['tf','pdbid','pairid','param','count']])
    
    b.to_csv("./Data/final_local_and_global_distortion_all_sigma%.1f.csv"%sigma,index=False)
