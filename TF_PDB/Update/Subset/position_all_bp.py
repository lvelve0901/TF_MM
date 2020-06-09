
from Bio.Seq import complement
import pandas as pd
import numpy as np
import pylab as pl
import os, sys
from seq_gap import *
from pprint import pprint


pdbinfo = pd.read_csv("../Data/refine_pdbinfo_final.csv")
bseqinfo = pd.read_csv("../Data/tf_bseq.csv")
contactinfo = pd.read_csv("../Data/tf_contact.csv")
DNAproDB = pd.read_csv("../Data/tf_DNAproDB.csv")

pdbinfo['bseq1_full'] = bseqinfo['bseq1']
pdbinfo['bseq2_full'] = bseqinfo['bseq2']
pdbinfo['bidx1_full'] = bseqinfo['bidx1']
pdbinfo['bidx2_full'] = bseqinfo['bidx2']


#pdblist = ['1an2','1dux','1hjb','1k79','1nkp','1nlw','1p47','1qne','1r4r','3gat','3kz8','5t00','5u01','5zk1','6qhd']
pdblist = ['1an2','1dux','1hjb','1k79','1nkp','1p47','1qne','1r4r','3kz8','5t00','5u01','5zk1']

tf_seqs = {
        '1an2':[[7,14],'TAGATCACGTGAATG'],
        '1dux':[[2,8],'CCCTGTTCCGGTAT'],
        '1hjb':[[15,21],'TCTTATGTGGTTTTTG'],
        '1k79':[[4,9],'CCTGTTCCGGTAT'],
        '1nkp':[[6,12],'TGCCACCTGGTGGCCACGTGCC'],
        '1p47':[[1,10],'GTGGCGTGGGCGATA'],
        '1qne':[[1,9],'ACTATAAAAAGTTC'],
        '1r4r':[[1,16],'CAGAACATGATGTTCTCA'],
        '3kz8':[[3,18],'AGACATGCCCGGGCATGCCTC'],
        '5t00':[[1,7],'GCAGCGCCCTCTACTGGCAGC'],
        '5u01':[[2,23],'CCTGGGGAATTTCCGGGA'],
        '5zk1':[[6,13],'CCGGTGACGTAAACG'],
}

pdbinfo = pdbinfo.loc[pdbinfo.pdbid.isin(pdblist)].reset_index(drop=True)

df = pd.read_csv('../Data/tf_local.csv')
df_range = pd.read_csv("../../../Survey/Survey_goldWC/stem/free/Canonical_all_bp_params.csv")



# find gap in helices
gappdb = pdbinfo.loc[~(pdbinfo.bseq1 == pdbinfo.bseq1_full)]['pdbid'].tolist()


output = [['tf','resi','pdbid','pairid','isHbond','isVdw']]


for i in range(len(pdbinfo)):

    pdbid = pdbinfo['pdbid'].ix[i]
    name = pdbinfo['name'].ix[i]


    bseq1 = pdbinfo['bseq1'].ix[i]
    bseq2 = pdbinfo['bseq2'].ix[i]

    bseq1_full = pdbinfo['bseq1_full'].ix[i]
    bseq2_full = pdbinfo['bseq2_full'].ix[i]

    bidx1 = list(map(int,pdbinfo['bidx1'].ix[i].split('|')))
    bidx2 = list(map(int,pdbinfo['bidx2'].ix[i].split('|')))

    bidx1_full = list(map(int,pdbinfo['bidx1_full'].ix[i].split('|')))
    bidx2_full = list(map(int,pdbinfo['bidx2_full'].ix[i].split('|')))
    
    if pdbid in gappdb:

        seqrange = np.array(range(len(bseq1_full)))
        subdf = df.loc[df['pdbid'] == pdbid].reset_index(drop=True)
        pos = find_gap(bidx1,bidx1_full)

    else:

        seqrange = np.array(range(len(bseq1)))
        subdf = df.loc[df['pdbid'] == pdbid].reset_index(drop=True)

    resi_list = np.array(range(tf_seqs[pdbid][0][0],tf_seqs[pdbid][0][1]+1)) - 1

    for resi in resi_list:
    #for resi in range(len(subdf)):

        pairid = subdf.ix[resi]['pair_id']

        pairid_list = pairid.split('_')
        nt1_code = pairid_list[2] + '.' + pairid_list[3] + pairid_list[4]
        nt2_code = pairid_list[5] + '.' + pairid_list[6] + pairid_list[7]
        contacts_hbond = DNAproDB.loc[((DNAproDB.pdbid == pdbid) & (DNAproDB.nt_code.isin([nt1_code,nt2_code])) & (DNAproDB.int_type == 'hbond') & (DNAproDB.nt_moiety.isin(['minorgw','majorgw','base'])))]
        contacts_vdw = DNAproDB.loc[((DNAproDB.pdbid == pdbid) & (DNAproDB.nt_code.isin([nt1_code,nt2_code])) & (DNAproDB.int_type == 'vdw') & (DNAproDB.nt_moiety.isin(['minorgw','majorgw','base'])))]
        num_hbond = len(contacts_hbond)
        num_vdw = len(contacts_vdw)

        #contact_num1 = contactinfo.loc[contactinfo.pair_id == pairid]['num_sc_1'].values
        #contact_num2 = contactinfo.loc[contactinfo.pair_id == pairid]['num_sc_2'].values
        #contact_num = contact_num1 + contact_num2
        if num_hbond == 0:
            isHbond = False
        else:
            isHbond = True

        if num_vdw == 0:
            isVdw = False
        else:
            isVdw = True

        output.append([name,resi+1,pdbid,pairid,isHbond,isVdw])


output_df = pd.DataFrame(output[1:],columns=output[0])
print(output_df)
output_df.to_csv("position_all_bp.csv",index=False,float_format='%.3f')
