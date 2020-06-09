#!/usr/bin/python

import os
import sys
import pylab as pl
import numpy as np
import pandas as pd
from commontool import read
from pdblib.abg import *

inpdir = './Helice_update'
bfnhlx = './iBformDNA_final.pdb'
bhlx = Mol(bfnhlx)

pdblist = [pdbf[:-4] for pdbf in os.listdir(inpdir)]
pdblist.sort()

bdic = {'DA':'A','DT':'T','DG':'G','DC':'C',
         '5CM':'c','5HC':'c','5FC':'c','1CC':'c',
         'IGU':'a','DI':'g','5IU':'u'}


output_bseq = [['pdbid','bseq1','bseq2','bidx1','bidx2']]


for idx, pdbid in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))

    mol = Mol('./%s/%s.pdb'%(inpdir,pdbid))
   
    bseq1 = []
    bseq2 = []
    bidx1 = []
    bidx2 = []

    for i in range(len(mol.segs[0].reses)):
        bseq1.append(bdic[mol.segs[0].reses[i].name])
        bseq2.append(bdic[mol.segs[1].reses[i].name])
        bidx1.append(str(mol.segs[0].reses[i].resi))
        bidx2.append(str(mol.segs[1].reses[i].resi))

    bseq1 = "".join(bseq1)
    bseq2 = "".join(bseq2[::-1])
    bidx1 = "|".join(bidx1)
    bidx2 = "|".join(bidx2[::-1])

    output_bseq.append([pdbid,bseq1,bseq2,bidx1,bidx2])

df = pd.DataFrame(output_bseq[1:],columns=output_bseq[0])
df.to_csv('./Data/tf_bseq.csv',index=False)

