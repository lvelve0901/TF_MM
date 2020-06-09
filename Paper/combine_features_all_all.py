
import sys
import numpy as np
import pandas as pd

param_list = ['shear','stretch','stagger','buckle','propeller','opening','c1c1_dist','shift','slide','rise','tilt','roll','twist','minorgw(-2.0)','minorgw(-1.0)','minorgw(0.0)','minorgw(1.0)','majorgw(-2.0)','majorgw(-1.0)','majorgw(0.0)','majorgw(1.0)','beta(-1.0)','beta(0.0)','beta(1.0)','gamma(-1.0)','gamma(0.0)','gamma(1.0)','zeta(-1.0)','zeta(0.0)','zeta(1.0)']

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

tf_mm_pos = {
    '1an2':[9,12],
    '1dux':[4,6],
    '1hjb':[17],
    '1k79':[5,6,7],
    '1nkp':[7,10],
    '1p47':[1,9],
    '1qne':[5,8,9],
    '1r4r':[3,13],
    '3kz8':[4,5,14,15],
    '5t00':[2,5,6],
    '5u01':[2,3,22,23],
    '5zk1':[7,12],
}

#tf_seqs = {
#        '1an2':[
#                [9,'TAGATCA*CGTGAATG'],
#                [12,'TAGATCA*CGTGAATG'],
#                ],
#        '1dux':[
#                [4,'CCCTGTTCCG*GTAT'],
#                [6,'CCCTGTTC*CGGTAT'],
#                ],
#        '1hjb':[
#                [17,'TCTTATGT*GGTTTTTG'],
#                ],
#        '1k79':[
#                [5,'CCTGTTCCG*GTAT'],
#                [6,'CCTGTTCC*GGTAT'],
#                [7,'CCTGTTC*CGGTAT'],
#                ],
#        '1nkp':[
#                [7,'TGCCACCTGGTGGCCA*CGTGCC'],
#                [10,'TGCCACCTGGTGGCCA*CGTGCC'],
#                ],
#        '1p47':[
#                [1,'GTG*GCGTGGGCGATA'],
#                [9,'GTGGCGTGGGC*GATA'],
#                ],
#        '1qne':[
#                [5,'ACTATA*AAAAGTTC'],
#                [8,'ACTATAAAAA*GTTC'],
#                [9,'ACTATAAAAAG*TTC'],
#                ],
#        '1r4r':[
#                [3,'CAGA*ACATGATGTTCTCA'],
#                [13,'CAGAACATGATGTT*CTCA'],
#                ],
#        '3kz8':[
#                [4,'AGACA*TGCCCGGGCATGCCTC'],
#                [5,'AGACAT*GCCCGGGCATGCCTC'],
#                [14,'AGACATGCCCGGGCA*TGCCTC'],
#                [15,'AGACATGCCCGGGCAT*GCCTC'],
#                ],
#        '5t00':[
#                [2,'GCAG*CGCCCTCTACTGGCAGC'],
#                [5,'GCAGCGC*CCTCTACTGGCAGC'],
#                [6,'GCAGCGCC*CTCTACTGGCAGC'],
#                ],
#        '5u01':[
#                [2,'CCTGGGGAATTTCCG*GGA'],
#                [3,'CCTGGGGAATTTCC*GGGA'],
#                [22,'CCTGGGGAATTTCC*GGGA'],
#                [23,'CCTGGGGAATTTCCG*GGA'],
#                ],
#        '5zk1':[
#                [7,'CCGGTG*ACGTAAACG'],
#                [12,'CCGGTG*ACGTAAACG'],
#                ],
#}





df1 = pd.read_csv("../TF_PDB/Update/Subset/features_all_bp.csv")
df2 = pd.read_csv("../TF_PDB/Update/Subset/features_all_step.csv")
df1['flank'] = np.nan
df2['flank'] = np.nan
df3 = pd.read_csv("../TF_PDB/Update/Subset/features_all_groove_flank.csv")
df4 = pd.read_csv("../TF_PDB/Update/Subset/features_all_abg_flank.csv")

df = pd.concat([df1,df2,df3,df4],axis=0)


df = df.sort_values(['tf','pdbid','resi']).reset_index(drop=True)

df = df[['tf','flank','pdbid','resi','pairid','isHbond','isVdw','param','value','std_refine']]
df['param'] = df['param'].apply(lambda x: x.split('_')[0] if x in ['beta_h_b','gamma_h_b','zeta_h_b'] else x)

df['param_flank'] = df['param'] + "(" + df['flank'].map(str) + ")"
df['param_flank'] = df.apply(lambda x: x['param_flank'].split("(")[0] if np.isnan(x['flank']) else x['param_flank'], axis=1)

df['isMM'] = False

for i in range(len(df)):
    pdbid = df.ix[i]['pdbid']
    resi = df.ix[i]['resi']
    if resi in tf_mm_pos[pdbid]:
        df.at[i,'isMM'] = True

b = df.groupby(['tf','pdbid','resi','pairid','isMM','isHbond','isVdw'])['param_flank'].agg(lambda x: list(x))
b = b.reset_index()
c = df.groupby(['tf','pdbid','resi','pairid','isMM','isHbond','isVdw'])['std_refine'].agg(lambda x: list(x))
c = c.reset_index()
result = b[['tf','pdbid','resi','pairid','isMM','isHbond','isVdw']]

for param in param_list:
    result[param] = np.nan

for i in range(len(result)):
    
    param_flanks = b.ix[i]['param_flank']
    std_refines = c.ix[i]['std_refine']
    for param in param_list:
        
        if param in param_flanks:
            param_index = param_flanks.index(param)
            std_refine = std_refines[param_index]
            result.at[i,param] = std_refine 

result = result.sort_values(['tf','pdbid','resi']).reset_index(drop=True)
print(result)
result.to_csv("./Data/features_raluca.csv",index=False,float_format='%.2f')
