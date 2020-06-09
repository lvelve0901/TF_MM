
from Bio.Seq import complement
import pandas as pd
import numpy as np
import pylab as pl
import os, sys
from seq_gap import *
from pprint import pprint

sigma = 1.0
#sigma = 0.5

pdbinfo = pd.read_csv("../Data/refine_pdbinfo_final.csv")
bseqinfo = pd.read_csv("../Data/tf_bseq.csv")
mminfos = []
for i in range(-1,2):
    mminfos.append(pd.read_csv("../../../Paper/Data/cutoff_mm_value_abg_sigma%.1f_flank_%s.csv"%(sigma,i)))
contactinfo = pd.read_csv("../Data/tf_contact.csv")
DNAproDB = pd.read_csv("../Data/tf_DNAproDB.csv")

pdbinfo['bseq1_full'] = bseqinfo['bseq1']
pdbinfo['bseq2_full'] = bseqinfo['bseq2']
pdbinfo['bidx1_full'] = bseqinfo['bidx1']
pdbinfo['bidx2_full'] = bseqinfo['bidx2']


#pdblist = ['1an2','1dux','1hjb','1k79','1nkp','1nlw','1p47','1qne','1r4r','3gat','3kz8','5t00','5u01','5zk1','6qhd']
pdblist = ['1an2','1dux','1hjb','1k79','1nkp','1p47','1qne','1r4r','3kz8','5t00','5u01','5zk1']

tf_seqs = {
        '1an2':[
                [9,'TAGATCA*CGTGAATG'],
                #[12,'TAGATCA*CGTGAATG'],
                ],
        '1dux':[
                [4,'CCCTGTTCCG*GTAT'],
                [6,'CCCTGTTC*CGGTAT'],
                ],
        '1hjb':[
                #[14,'TCTTA*TGTGGTTTTTG'],
                [17,'TCTTATGT*GGTTTTTG'],
                ],
        '1k79':[
                [5,'CCTGTTCCG*GTAT'],
                [6,'CCTGTTCC*GGTAT'],
                [7,'CCTGTTC*CGGTAT'],
                ],
        '1nkp':[
                [7,'TGCCACCTGGTGGCCA*CGTGCC'],
                #[10,'TGCCACCTGGTGGCCA*CGTGCC'],
                ],
        #'1nlw':[
        #        [7,'GTAGCA*CGCGTAAC'],
        #        [10,'GTAGCA*CGCGTAAC'],
        #        ],
        '1p47':[
                [1,'GTG*GCGTGGGCGATA'],
                [9,'GTGGCGTGGGC*GATA'],
                ],
        '1qne':[
                #[2,'ACT*ATAAAAAGTTC'],
                #[3,'ACTA*TAAAAAGTTC'],
                [5,'ACTATA*AAAAGTTC'],
                [8,'ACTATAAAAA*GTTC'],
                [9,'ACTATAAAAAG*TTC'],
                ],
        '1r4r':[
                [3,'CAGA*ACATGATGTTCTCA'],
                [13,'CAGAACATGATGTT*CTCA'],
                ],
        #'3gat':[
        #        [5,'CGTCG*GATATCCGGT'],
        #        ],
        '3kz8':[
                [4,'AGACA*TGCCCGGGCATGCCTC'],
                #[5,'AGACAT*GCCCGGGCATGCCTC'],
                [14,'AGACATGCCCGGGCA*TGCCTC'],
                #[15,'AGACATGCCCGGGCAT*GCCTC'],
                ],
        '5t00':[
                [2,'GCAG*CGCCCTCTACTGGCAGC'],
                [5,'GCAGCGC*CCTCTACTGGCAGC'],
                [6,'GCAGCGCC*CTCTACTGGCAGC'],
                ],
        '5u01':[
                #[2,'CCTGGGGAATTTCCG*GGA'],
                #[3,'CCTGGGGAATTTCC*GGGA'],
                [22,'CCTGGGGAATTTCC*GGGA'],
                [23,'CCTGGGGAATTTCCG*GGA'],
                ],
        '5zk1':[
                #[4,'CCG*GTGACGTAAACG'],
                [7,'CCGGTG*ACGTAAACG'],
                #[12,'CCGGTG*ACGTAAACG'],
                #[15,'CCG*GTGACGTAAACG'],
                ],
        #'6qhd':[
        #        [8,'AAGTTCCT*GGAATTT'],
        #        [10,'AAGTTC*CTGGAATTT'],
        #        ],
}

mmdf = {
        '1an2':[
                [9,['GT']],
                #[12,['GT']],
                ],
        '1dux':[
                [4,['CT']],
                [6,['GA','GA*syn']],
                ],
        '1hjb':[
                #[14,['TT','CT']],
                [17,['AC','AA']],
                ],
        '1k79':[
                [5,['CT']],
                [6,['GG*syn','G*Gsyn']],
                [7,['GA','GA*syn','GT']],
                ],
        '1nkp':[
                [7,['GT']],
                #[10,['GT']],
                ],
        #'1nlw':[
        #        [7,['GT']],
        #        [10,['GT']],
        #        ],
        '1p47':[
                [1,['CT']],
                [9,['GG*syn','G*Gsyn']],
                ],
        '1qne':[
                #[2,['CC']],
                #[3,['CC']],
                [5,['TT']],
                [8,['TT']],
                [9,['CC','CT','AC']],
                ],
        '1r4r':[
                [3,['TT','GT']],
                [13,['GA','GA*syn']],
                ],
        #'3gat':[
        #        [5,['CT','AC']],
        #        ],
        '3kz8':[
                [4,['TT','CT']],
                #[5,['TT','CT']],
                [14,['TT','CT']],
                #[15,['TT','CT']],
                ],
        '5t00':[
                [2,['CC']],
                [5,['GT']],
                [6,['GG*syn','G*Gsyn']],
                ],
        '5u01':[
                #[2,['CT','AC','CC']],
                #[3,['GG*syn','G*Gsyn']],
                [22,['GG*syn','G*Gsyn']],
                [23,['CT','AC','CC']],
                ],
        '5zk1':[
                #[4,['CT','TT']],
                [7,['CC']],
                #[12,['CC']],
                #[15,['CT','TT']],
                ],
        #'6qhd':[
        #        [8,['AC']],
        #        [10,['GA']],
        #        ],
        }

pprint(mmdf)

pdbinfo = pdbinfo.loc[pdbinfo.pdbid.isin(pdblist)].reset_index(drop=True)

df = pd.read_csv('../Data/tf_abg.csv')
df_range = pd.read_csv("../../../Survey/Survey_bend/free/Bend_Canonical_all_abg_params.csv")

bp_params = ['beta_h_b','gamma_h_b','zeta_h_b']


# generate B-DNA envelope
envelopes = []
for i in range(0,3):

    ymean = df_range.ix[i]['param_mean']
    std = df_range.ix[i]['param_sigma']
    ymax = df_range.ix[i]['param_max']
    ymin = df_range.ix[i]['param_min']
    
    if sigma == 0.5:
        ylower1 = ymean - 0.5*std
        yupper1 = ymean + 0.5*std
    elif sigma == 1.0:
        ylower1 = ymean - std
        yupper1 = ymean + std
    ylower2 = ymean - 2*std
    yupper2 = ymean + 2*std
    ylower3 = ymean - 3*std
    yupper3 = ymean + 3*std

    envelopes.append([ymean, ymin, ymax, ylower1, yupper1, ylower2, yupper2, ylower3, yupper3])

# find gap in helices
gappdb = pdbinfo.loc[~(pdbinfo.bseq1 == pdbinfo.bseq1_full)]['pdbid'].tolist()

output = [['tf','sequence_star','strand1','strand2','resi','flank','pdbid','pairid','isHbond','isVdw','param','value','mismatch']]


# generate plot for each bp params
for idx, param in enumerate(bp_params):

    print("--- Working on %s [%d/%d] ---"%(param,idx,len(bp_params)))

    ymean, ymin, ymax, ylower1, yupper1, ylower2, yupper2, ylower3, yupper3 = envelopes[idx]

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
            param_value = subdf[param].tolist()

        else:

            seqrange = np.array(range(len(bseq1)))
            subdf = df.loc[df['pdbid'] == pdbid].reset_index(drop=True)
            param_value = subdf[param].tolist()
        
        targets = mmdf[pdbid] 
        seqs = tf_seqs[pdbid]

        for target, seq in zip(targets,seqs):
            resi = target[0]
            sequence_star = seq[1]
            mm_pos = sequence_star.index("*") - 1
            sequence = sequence_star.replace("*","")
            sequence_comp = complement(sequence)
            wt_bp = sequence[mm_pos]+sequence_comp[mm_pos]

            pairid_o = subdf.ix[resi]['pair_id']
            resi_list = [resi-1,resi,resi+1]
            pairids = []
            values = []
            flanks = []
            for i, r in enumerate(resi_list):
                pairid = subdf.ix[r]['pair_id']
                value = param_value[r]
                if np.isnan(value) and np.isnan(subdf.loc[subdf.pair_id == pairid][param].values[0]):
                    continue
                pairids.append(pairid)
                values.append(value)
                flanks.append(i-1)

            #print(pairids)
            #print(values)
            #print(flanks)

            mmlist = target[1]

            pairid_list = pairid_o.split('_')
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
            
            for pairid, value, flank, mminfo in zip(pairids,values,flanks,mminfos):

                if value != subdf.loc[subdf.pair_id == pairid][param].values[0]:
                    print("ERROR(): Mis-alignment of sequence!")
                    print(value)
                    print(subdf.loc[subdf.pair_id == pairid][param].values[0])
                    sys.exit()
                

                for mm in mmlist:
                    isSign = mminfo.loc[mminfo.mismatch == mm][param[:-4]].values
                    mm_correct = mm.replace("*","")[0:2]
                    wt_bp = wt_bp

                    if mm_correct[1] == wt_bp[1]:
                        strand1 = sequence
                        strand1 = sequence[:mm_pos] + mm_correct[0] + sequence[mm_pos+1:]
                        strand2 = sequence_comp
                        
                    elif mm_correct[0] == wt_bp[1]:
                        mm_correct = mm_correct[::-1]
                        strand1 = sequence
                        strand1 = sequence[:mm_pos] + mm_correct[0] + sequence[mm_pos+1:]
                        strand2 = sequence_comp
                    else:
                        strand1 = sequence[:mm_pos] + mm_correct[0] + sequence[mm_pos+1:]
                        strand2 = sequence_comp[:mm_pos] + mm_correct[1] + sequence_comp[mm_pos+1:]

                    if isSign == -1.0:
                        if value <= ylower1:
                            output.append([name,sequence_star,strand1,strand2,resi,flank,pdbid,pairid_o,isHbond,isVdw,param,value,mm])
                    elif isSign == 1.0:
                        if value >= yupper1:
                            output.append([name,sequence_star,strand1,strand2,resi,flank,pdbid,pairid_o,isHbond,isVdw,param,value,mm])
                    elif isSign == 2.0:
                        if value <= ylower1 or value >= yupper1:
                            output.append([name,sequence_star,strand1,strand2,resi,flank,pdbid,pairid_o,isHbond,isVdw,param,value,mm])

output_df = pd.DataFrame(output[1:],columns=output[0])
print(output_df)
output_df.to_csv("mimicry_abg_flank_sigma%.1f.csv"%sigma,index=False,float_format='%.3f')
