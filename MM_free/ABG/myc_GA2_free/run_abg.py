#!/usr/bin/python

import os
from os import system,listdir
from pdblib.base import *
from Curves import *

inpdir = '/mnt/d4/hs189/TF/myc_max/myc_GA2_free'
pdblist = []
for file in os.listdir(inpdir): #read input directory
    if file.endswith(".pdb"):
        pdblist.append(file)

pdblist.sort(key=lambda x: int(x.split('-')[-1][:-4]))

for pdb in pdblist:
    pdbfile = str(pdb)
    pdbname = str(pdbfile)[:-4]
    system('cp %s/%s temp.pdb'%(inpdir,pdb))
    mol = Mol('temp.pdb')
    reses = getreses(mol)
    for res in reses:
    	if res.name == 'DA' or res.name == 'RA' or res.name == 'ADE':
    		res.name = 'A'
    	elif res.name == 'DT' or res.name == 'T' or res.name == 'URA' or res.name == 'THY':
    		res.name = 'U'
    	elif res.name == 'DC' or res.name == 'RC' or res.name == 'CYT':
    		res.name = 'C'
    	elif res.name == 'DG' or res.name == 'RG' or res.name == 'GUA':
    		res.name = 'G'
    mol.segs[0].chid = 'A'
    mol.segs[1].chid = 'B'
    mol.write('output.pdb')
    system('mv output.pdb %s'%pdbfile)
    inp = []
    inp.append("A.1:7    B.30:36 '(((.(((' '))).)))'")
    inp.append("A.2:8    B.29:35 '(((.(((' '))).)))'")
    inp.append("A.3:9    B.28:34 '(((.(((' '))).)))'")
    inp.append("A.4:10   B.27:33 '(((.(((' '))).)))'")
    inp.append("A.5:11   B.26:32 '(((.(((' '))).)))'")
    inp.append("A.6:12   B.25:31 '(((.(((' '))).)))'")
    inp.append("A.7:13   B.24:30 '(((.(((' '))).)))'")
    inp.append("A.8:14   B.23:29 '(((.(((' '))).)))'")
    inp.append("A.9:15   B.22:28 '(((.(((' '))).)))'")
    inp.append("A.10:16  B.21:27 '(((.(((' '))).)))'")
    inp.append("A.11:17  B.20:26 '(((.(((' '))).)))'")
    inp.append("A.12:18  B.19:25 '(((.(((' '))).)))'")
    for i in range(len(inp)):
    	system('perl allAtomMeas_DNA.pl %s %s %s'%(pdbfile,inp[i],1))  #incorporate minor groove width to the ABG output
    system('rm temp.pdb %s'%pdbfile)
