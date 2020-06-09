
import sys
import pandas as pd
import numpy as np

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
for sigma in [1.0]:

    df1 = pd.read_csv("../TF_PDB/Update/Subset/mimicry_bp_sigma%.1f.csv"%sigma)
    df2 = pd.read_csv("../TF_PDB/Update/Subset/mimicry_step_sigma%.1f.csv"%sigma)
    df1['flank'] = np.nan
    df2['flank'] = np.nan
    df3 = pd.read_csv("../TF_PDB/Update/Subset/mimicry_groove_flank_sigma%.1f.csv"%sigma)
    df4 = pd.read_csv("../TF_PDB/Update/Subset/mimicry_abg_flank_sigma%.1f.csv"%sigma)
    
    df = pd.concat([df1,df2,df3,df4],axis=0)
    df = df.sort_values(['tf','pdbid','resi']).reset_index(drop=True)
    
    df = df[['tf','sequence_star','strand1','strand2','flank','pdbid','resi','pairid','isHbond','isVdw','param','value','mismatch']]
    print(df)
    df['param'] = df['param'].apply(lambda x: x.split('_')[0] if x in ['beta_h_b','gamma_h_b','zeta_h_b'] else x)
    df.to_csv("./Data/final_local_and_global_mimicry_value_sigma%.1f.csv"%sigma,index=False,float_format='%.3f')
    
    df['param_flank'] = df['param'] + "(" + df['flank'].map(str) + ")"
    df['param_flank'] = df.apply(lambda x: x['param_flank'].split("(")[0] if np.isnan(x['flank']) else x['param_flank'], axis=1)
    
    b = df.groupby(['tf','sequence_star','strand1','strand2','pdbid','resi','pairid','isHbond','isVdw','mismatch'])['param_flank'].agg(lambda x: ' '.join(list(x)))
    
    b = b.reset_index()
    b['strand1'] = "5'-" + b['strand1'] + "-3'"
    b['strand2'] = "3'-" + b['strand2'] + "-5'"
    
    b = b.rename(columns={'param_flank':'param'})
    
    print(b[['tf','sequence_star','strand1','strand2','pdbid','resi','pairid','mismatch','param']])
    
    b.to_csv("./Data/final_local_and_global_mimicry_sigma%.1f.csv"%sigma,index=False)
