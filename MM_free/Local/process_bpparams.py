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
target_nts = ['A.DA8','B.DT29']
target_nts = ['A.DC8','B.DG29']

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
param_list = [['pdbid','pdbname','bp_id','bpname','name','LW','nt1_chain','nt1_i','nt1_name','nt2_chain','nt2_i','nt2_name','nt1_chi','nt2_chi','nt1_phase','nt2_phase','nt1_ampli','nt2_ampli','C1C1_dist','hbonds_num','shear','stretch','stagger','buckle','propeller','opening']]


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

    nts_idx = {}  #nts index dictionary with nt_id as key
    for i, nt in enumerate(nts):
        nts_idx[nt['nt_id']] = i
    if len(nts_idx) != len(nts):  #error: there are nucleotides labeled same
        sys.exit()

    for bp in bps:
        nt1_id = bp['nt1']
        nt2_id = bp['nt2']
        if nt1_id not in target_nts or nt2_id not in target_nts:
            continue
        nt1 = nts[nts_idx[nt1_id]]
        nt2 = nts[nts_idx[nt2_id]]
        
        bpname = bp['bp'].replace("-","")
        name = bp['name']
        LW = bp['LW']
        hbonds_num = bp['hbonds_num']
        C1C1_dist = bp['C1C1_dist']
        
        shear, stretch, stagger, buckle, propeller, opening = bp['bp_params']

    	nt1_i = nt1['index']
    	nt1_name = nt1['nt_name']
        nt1_chain = nt1['chain_name']
        nt1_resc = nt1['nt_code']
        nt1_resi = nt1_id.split('.')[1].replace(nt1_name,"").replace("/","").split('^')[0] 

        if int(nt1_resi) != int(nt1_i):
            print "(ERROR): Inconsistent nucleotide id! %s %s"%(nt1_resi,nt1_i)
            sys.exit()

    	nt1_alpha = nt1['alpha']
    	nt1_beta = nt1['beta']
    	nt1_gamma = nt1['gamma']
    	nt1_delta = nt1['delta']
    	nt1_epsilon = nt1['epsilon']
    	nt1_zeta = nt1['zeta']
        nt1_e_z = nt1['epsilon_zeta']
    	nt1_chi = nt1['chi']
        nt1_v0 = nt1['v0']
    	nt1_v1 = nt1['v1']
    	nt1_v2 = nt1['v2']
    	nt1_v3 = nt1['v3']
    	nt1_v4 = nt1['v4']
        nt1_phase = nt1['phase_angle']
        nt1_ampli = nt1['amplitude']
        nt1_pucker = nt1['puckering']

    	nt2_i = nt2['index']
    	nt2_name = nt2['nt_name']
        nt2_chain = nt2['chain_name']
        nt2_resc = nt2['nt_code']
        nt2_resi = nt2_id.split('.')[1].replace(nt2_name,"").replace("/","").split('^')[0] 
        
        if int(nt2_resi) != int(nt2_i):
            print "(ERROR): Inconsistent nucleotide id! %s %s"%(nt2_resi,nt2_i)
            sys.exit()

    	nt2_alpha = nt2['alpha']
    	nt2_beta = nt2['beta']
    	nt2_gamma = nt2['gamma']
    	nt2_delta = nt2['delta']
    	nt2_epsilon = nt2['epsilon']
    	nt2_zeta = nt2['zeta']
        nt2_e_z = nt2['epsilon_zeta']
    	nt2_chi = nt2['chi']
        nt2_v0 = nt2['v0']
    	nt2_v1 = nt2['v1']
    	nt2_v2 = nt2['v2']
    	nt2_v3 = nt2['v3']
    	nt2_v4 = nt2['v4']
        nt2_phase = nt2['phase_angle']
        nt2_ampli = nt2['amplitude']
        nt2_pucker = nt2['puckering']

        bp_id = pdbid + "_" + bpname + "_" + nt1_chain + "_" + nt1_name + "_" + nt1_resi + "_" + nt2_chain + "_" + nt2_name + "_" + nt2_resi

        param_list.append([pdbid,pdbname,bp_id,bpname,name,LW,nt1_chain,nt1_i,nt1_name,nt2_chain,nt2_i,nt2_name,nt1_chi,nt2_chi,nt1_phase,nt2_phase,nt1_ampli,nt2_ampli,C1C1_dist,hbonds_num,shear,stretch,stagger,buckle,propeller,opening])

        
df = pd.DataFrame(param_list[1:],columns=param_list[0])
df.to_csv("./Data/%s_bp.csv"%oupname,index=False)

