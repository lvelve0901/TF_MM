#!/usr/bin/python

import json # Handle JSON DSSR files
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import learnna_json as lna_json
from pdblib.base import *


pdb_dir = '/mnt/hs189/TF/myc_max/myc_CG_free/'
oupname = 'myc_CG_free'
target1_nts = ['A.DC7','B.DG30']
#target2_nts = ['A.DA8','B.DT29']
target2_nts = ['A.DC8','B.DG29']

json_dir = './Json/Json_%s'%oupname

if os.path.isfile(json_dir) is False: 
    os.system("mkdir -p %s"%json_dir)

pdblist = []  #list to store PDB file names
#Walk data_directory and get list of pdb names
for file in os.listdir(pdb_dir):
	if file.endswith(".pdb"):
		pdblist.append(file)

pdblist.sort(key=lambda x: int(x.split('-')[-1][:-4]))

# params list:
param_list = [['pdbid','pdbname','step','bp_id','bpname','name','LW','nt1_chain','nt1_i','nt1_name','nt2_chain','nt2_i','nt2_name','shift','slide','rise','tilt','roll','twist','x_dis','y_dis','hrise','inclin','tip','htwist','minorgw1','majorgw1']]


#Loop over list of files in data dir, process with DSSR to get json file
for idx, pdb in enumerate(pdblist):
    pdbname = str(pdb)[:-4]
    pdbid = pdb.split('-')[0] 
    pdb_f = os.path.join(pdb_dir,pdb)  #path to pdb file
    json_f = os.path.join(json_dir,pdb.replace(".pdb",".json"))  #path to json file
    print ("--- Working on [%s] (%d of %d) ---"%(pdb,idx+1,len(pdblist)))
    if os.path.isfile(json_f) is False:
        if os.path.isfile(pdb_f):  #DSSR convert PDB to json files
    	    os.system("x3dna-dssr --json --more -i=%s -o=%s"%(pdb_f,json_f))
    	    os.system("x3dna-dssr --cleanup")

    najson = lna_json.NA_JSON()  #initialize class objects
    with open(json_f) as json_data:  #read each json file
    	data = json.load(json_data)
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file

    nts = najson.json_file['nts']
    bps = najson.json_file['pairs']
    hels = najson.json_file['helices']

    if len(hels) > 1:
        for h in hels:
            print h['pairs'][0]['nt1']
        print("WARNING(): Find more than one helix\n")
        continue

    hels = hels[0]['pairs']

    nts_idx = {}  #nts index dictionary with nt_id as key
    for i, nt in enumerate(nts):
        nts_idx[nt['nt_id']] = i
    if len(nts_idx) != len(nts):  #error: there are nucleotides labeled same
        sys.exit()

    for hel in hels:
        nt1_id = hel['nt1']
        nt2_id = hel['nt2']
        step = 0
        if nt1_id in target1_nts and nt2_id in target1_nts:
            step = 1
        elif nt1_id in target2_nts and nt2_id in target2_nts:
            step = 2
        elif nt1_id not in target1_nts and nt1_id not in target2_nts and nt2_id not in target1_nts and nt2_id not in target2_nts:
            continue
        else:
            print("ERROR(): Wrong match of base pairs\n")
            print(nt1_id + "_+_" + nt2_id)
            #sys.exit()

        nt1 = nts[nts_idx[nt1_id]]
        nt2 = nts[nts_idx[nt2_id]]
        
        bpname = hel['bp'].replace("-","")
        name = hel['name']
        LW = hel['LW']
        
        shift, slide, rise, tilt, roll, twist = hel['step_params']
        x_dis, y_dis, hrise, inclin, tip, htwist = hel['heli_params']
        minorgw1, minorgw2, majorgw1, majorgw2 = hel['groove_widths']

    	nt1_i = nt1['index']
    	nt1_name = nt1['nt_name']
        nt1_chain = nt1['chain_name']
        nt1_resc = nt1['nt_code']
        nt1_resi = nt1_id.split('.')[1].replace(nt1_name,"").replace("/","").split('^')[0] 

        if int(nt1_resi) != int(nt1_i):
            print "(ERROR): Inconsistent nucleotide id! %s %s"%(nt1_resi,nt1_i)
            sys.exit()

    	nt2_i = nt2['index']
    	nt2_name = nt2['nt_name']
        nt2_chain = nt2['chain_name']
        nt2_resc = nt2['nt_code']
        nt2_resi = nt2_id.split('.')[1].replace(nt2_name,"").replace("/","").split('^')[0] 
        
        if int(nt2_resi) != int(nt2_i):
            print "(ERROR): Inconsistent nucleotide id! %s %s"%(nt2_resi,nt2_i)
            sys.exit()

        bp_id = pdbid + "_" + bpname + "_" + nt1_chain + "_" + nt1_name + "_" + nt1_resi + "_" + nt2_chain + "_" + nt2_name + "_" + nt1_resi

        param_list.append([pdbid,pdbname,step,bp_id,bpname,name,LW,nt1_chain,nt1_i,nt1_name,nt2_chain,nt2_i,nt2_name,shift,slide,rise,tilt,roll,twist,x_dis,y_dis,hrise,inclin,tip,htwist,minorgw1,majorgw1])

        
df = pd.DataFrame(param_list[1:],columns=param_list[0])
df.to_csv("./Data/%s_step.csv"%oupname,index=False)
