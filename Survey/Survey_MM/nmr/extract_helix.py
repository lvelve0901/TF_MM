#!/usr/bin/python

import json # Handle JSON DSSR files
import re
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import learnna_json as lna_json
from pdblib.base import *
from commontool import read

json_dir = './json'
os.system("mkdir -p %s"%json_dir)
pdb_dir = './raw_pdb'
pdb_dir_list = os.listdir(pdb_dir)

filelist = []  #list to store PDB file names
#Walk data_directory and get list of pdb names
for file in pdb_dir_list:
	if file.endswith(".pdb"):
		filelist.append(file)
filelist.sort()

output = [['pdbid','num_hel','hel_idx','bseq1','bseq2','bidx1','bidx2']]

#Loop over list of files in data dir, process with DSSR to get json file
for idx, pdb in enumerate(filelist):
    pdbname = str(pdb)[:-4]
    pdbid = pdb.split('_')[0][0:4]

    pdb_f = os.path.join(pdb_dir,pdb)  #path to pdb file
    json_f = os.path.join(json_dir,pdb.replace(".pdb",".json"))  #path to json file

    print ("--- Working on [%s] (%d of %d) ---"%(pdb,idx+1,len(filelist)))
    if os.path.isfile(pdb_f) is True and os.path.isfile(json_f) is False:  #DSSR convert PDB to json files
    	
        os.system("x3dna-dssr --json --more -i=%s -o=%s"%(pdb_f,json_f))
        os.system("cp dssr-helices.pdb ./helice/%s.pdb"%pdbid)
        os.system("x3dna-dssr --cleanup")


    najson = lna_json.NA_JSON()  #initialize class objects
    with open(json_f) as json_data:  #read each json file
    	data = json.load(json_data)
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file

    #create a nts dictionary with key of nt_code
    nts_idx = {}  #nts index dictionary with nt_id as key
    nts = najson.json_file['nts']
    for i, nt in enumerate(nts):
        nts_idx[nt['nt_id']] = i
    if len(nts_idx) != len(nts):  #error: there are nucleotides labeled same
        print "ERROR(): Potential same naming of nucleotide"
        sys.exit()

    dbn = najson.json_file['dbn']
    helices = najson.json_file['helices']
    
    num_hel = len(helices)
    #count helices
    if num_hel > 1:
        print("Number of helices: %d"%num_hel)
    
    for hel_idx, helice in enumerate(helices):
        #extract the fasta sequence
        bseq1 = helice['strand1']
        bseq2 = helice['strand2']


        bidx1 = [re.split('(\d+)',str(nts[nts_idx[hel['nt1']]]['nt_id'].split('.')[1]))[-2] for hel in helice['pairs']]
        bidx2 = [re.split('(\d+)',str(nts[nts_idx[hel['nt2']]]['nt_id'].split('.')[1]))[-2] for hel in helice['pairs']]
        bidx1 = "|".join(bidx1)
        bidx2 = "|".join(bidx2[::-1])
    
        output.append([pdbid,num_hel,hel_idx,bseq1,bseq2,bidx1,bidx2])

df = pd.DataFrame(output[1:],columns=output[0])
df = df.sort(['pdbid'])
df.to_csv('./pdbinfo_helice.csv',index=False)
