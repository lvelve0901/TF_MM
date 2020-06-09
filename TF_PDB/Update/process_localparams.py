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

json_dir = './Json_dssr'
os.system("mkdir -p %s"%json_dir)
pdb_dir = './RawPdb'
pdb_dir_list = os.listdir(pdb_dir)
pdbinfo = pd.read_csv("./Data/refine_pdbinfo_final.csv")
oupname = 'tf_local'

filelist = []  #list to store PDB file names
#Walk data_directory and get list of pdb names
for file in pdb_dir_list:
	if file.endswith(".pdb"):
		filelist.append(file)
filelist.sort()


#Output list
output_helis = [['pdbid','name','pair_id','bp_name','bp_type','bp_LW','nt1_chain','nt2_chain','nt1_name','nt2_name','nt1_resi','nt2_resi','nt1_phase','nt2_phase','nt1_chi','nt2_chi','shear','stretch','stagger','buckle','propeller','opening','shift','slide','rise','tilt','roll','twist','x_dis','y_dis','heli_rise','inclin','tip','heli_twist','minorgw1','c1c1_dist','n1n9_dist','bp_hbnum','hbond_desc1','hbond_len1','hbond_desc2','hbond_len2','hbond_desc3','hbond_len3','hbond_desc4','hbond_len4','hbond_desc5','hbond_len5']]

nt_can = ['DA', 'DG', 'DC', 'DT', '5IU', '5CM', '5HC', '5FC', '1CC', 'IGU', 'DI']
aa_can = ['ASP', 'ASN', 'ARG', 'GLU', 'ALA', 'GLN', 'SER', 'SEC', 'CYS', 'PRO', 'TYR', 'MET', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'PHE', 'THR', 'TRP', 'VAL']

bb_atoms = ["P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"]


#Loop over list of files in data dir, process with DSSR to get json file
for idx, pdb in enumerate(filelist):

    pdbname = str(pdb)[:-4]
    pdbid = pdb.split('_')[0][0:4]
    name = pdbinfo['name'][pdbinfo['pdbid'] == pdbname].values[0]
    pdb_f = os.path.join(pdb_dir,pdb)  #path to pdb file
    json_f = os.path.join(json_dir,pdb.replace(".pdb",".json"))  #path to json file
    print ("--- Working on [%s] (%d of %d) ---"%(pdb,idx+1,len(filelist)))
    najson = lna_json.NA_JSON()  #initialize class objects
    with open(json_f) as json_data:  #read each json file
    	data = json.load(json_data)
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file

    #extract the fasta sequence
    bseq = najson.json_file['helices'][0]['strand1']
    helice = najson.json_file['helices'][0]
    if len(bseq) != len(helice['pairs']):
        print "ERROR(): Inconsistent base pair number in the helix"
        sys.exit()

    #create a nts dictionary with key of nt_code
    nts_idx = {}  #nts index dictionary with nt_id as key
    nts = najson.json_file['nts']
    for i, nt in enumerate(nts):
        nts_idx[nt['nt_id']] = i
    if len(nts_idx) != len(nts):  #error: there are nucleotides labeled same
        print "ERROR(): Potential same naming of nucleotide"
        sys.exit()

    #create a bps dictionary with key of bp_code
    bps_idx = {}  #bps index dictionary with nt_id as key
    bps = najson.json_file['pairs']
    for i, bp in enumerate(bps):
        bps_idx[bp['nt1'] + "_" + bp['nt2']] = i
    if len(bps_idx) != len(bps):  #error: there are nucleotides labeled same
        print "ERROR(): Potential same naming of nucleotide"
        sys.exit()

    #create a hbonds dictionary with key of nt_code
    hbonds_bb_sc = []
    hbonds_sc = []
    hbonds = najson.json_file['hbonds']

    #for i, hbond in enumerate(hbonds):
    #    if hbond['residue_pair'] == 'nt:aa':
    #
    #        atom1_id = hbond['atom1_id']  #atom1 id
    #        atom2_id = hbond['atom2_id']  #atom2 id

    #        # Splitting the atom_id
    #        # atom_id separator
    #        atom1_name = atom1_id.split('@')[0]
    #        atom2_name = atom2_id.split('@')[0]

    #        res1_id = atom1_id.split('@')[1]  #resi1 id for example "A.G17"
    #        res2_id = atom2_id.split('@')[1]  #resi2 id

    #        chain1 = res1_id.split('.')[0]  #chain id 1
    #        chain2 = res2_id.split('.')[0]  #chain id 2
    #        resn_resi1 = res1_id.split('.')[1]  #resn_resi1
    #        resn_resi2 = res2_id.split('.')[1]  #resn_resi2
    #
    #        resi1 = re.split('(\-?\d+)',resn_resi1)[-2]
    #        resn1 = resn_resi1[:-len(resi1)]
    #        resi2 = re.split('(\-?\d+)',resn_resi2)[-2]
    #        resn2 = resn_resi2[:-len(resi2)]


    #        if resn1 in nt_can:
    #            hbonds_bb_sc.append(res1_id)
    #        elif resn2 in nt_can:
    #            hbonds_bb_sc.append(res2_id)
    #        
    #        if resn1 in nt_can and atom1_name not in bb_atoms:
    #            hbonds_sc.append(res1_id)
    #        elif resn2 in nt_can and atom2_name not in bb_atoms:
    #            hbonds_sc.append(res2_id)
    
    #print list(set(hbonds_bb_sc))
    #print list(set(hbonds_sc))

    #create a hbonds dictionary with key of bp_code
    

    #loop the helices
    for heli in helice['pairs'][1:]:  #unpack helical parameters
        nt1_id = heli['nt1']
        nt2_id = heli['nt2']


    #    if nt1_id in hbonds_bb_sc:
    #        contact_bb_sc_1 = True
    #    else:
    #        contact_bb_sc_1 = False

    #    if nt2_id in hbonds_bb_sc:
    #        contact_bb_sc_2 = True
    #    else:
    #        contact_bb_sc_2 = False

    #    if nt1_id in hbonds_sc:
    #        contact_sc_1 = True
    #    else:
    #        contact_sc_1 = False

    #    if nt2_id in hbonds_sc:
    #        contact_sc_2 = True
    #    else:
    #        contact_sc_2 = False



        bp_id = heli['nt1'] + "_" + heli['nt2']
        nt1 = nts[nts_idx[nt1_id]]
        nt2 = nts[nts_idx[nt2_id]]
        try:
            bp = bps[bps_idx[bp_id]]
        except KeyError, e:
            print "WARNING(): Cannot find bp_id:%s"%bp_id
            print "Try to swap nt1 and nt2"
            bp_id = heli['nt2'] + "_" + heli['nt1']
            bp = bps[bps_idx[bp_id]]

        nt1_chain = nt1['chain_name']
        nt2_chain = nt2['chain_name']
        nt1_name = nt1['nt_name']
        nt2_name = nt2['nt_name']
        nt1_resi = nt1_id.split('.')[1].replace(nt1_name,"")
        nt2_resi = nt2_id.split('.')[1].replace(nt2_name,"")
        nt1_phase = nt1['phase_angle']
        nt2_phase = nt2['phase_angle']
        nt1_chi = nt1['chi']
        nt2_chi = nt2['chi']
        bp_name = heli['bp'].replace("-","").replace("+","")
        bp_type = heli['name']
        bp_LW = heli['LW']
        c1c1_dist = bp['C1C1_dist']
        n1n9_dist = bp['N1N9_dist']
        bp_hbnum = bp['hbonds_num']
        hbond_desc = ['','','','','']
        hbond_len = [0,0,0,0,0]
        hbonds = bp['hbonds_desc'].split(',')
        for i, hbond in enumerate(hbonds):
            hbond_desc[i] = hbond.split('[')[0]
            hbond_len[i] = hbond.split('[')[1].replace(']','')
        pair_id = pdbid + "_" + bp_name + "_" + nt1_chain + "_" + nt1_name + "_" + nt1_resi + "_" + nt2_chain + "_" + nt2_name + "_" + nt2_resi 
        if 'groove_widths' in heli.keys():
            shear, stretch, stagger, buckle, propeller, opening = heli['bp1_params']
            shift, slide, rise, tilt, roll, twist = heli['step_params']
            x_dis, y_dis, heli_rise, inclin, tip, heli_twist = heli['heli_params']
            minorgw1, minorgw2, majorgw1, majorgw2 = heli['groove_widths']
            
            if bp_name in ['TA','CG','tA','cG','Ta','Cg','ta','cg','TG','CA','Tg','Ca','tG','cA','tg','ca']:
                shear = -shear
                buckle = -buckle
                shift = -shift
                tilt = -tilt
                y_dis = -y_dis
                tip = -tip

            output_helis.append([pdbid,name,pair_id,bp_name,bp_type,bp_LW,nt1_chain,nt2_chain,nt1_name,nt2_name,nt1_resi,nt2_resi,nt1_phase,nt2_phase,nt1_chi,nt2_chi,shear,stretch,stagger,buckle,propeller,opening,shift,slide,rise,tilt,roll,twist,x_dis,y_dis,heli_rise,inclin,tip,heli_twist,minorgw1,c1c1_dist,n1n9_dist,bp_hbnum,hbond_desc[0],hbond_len[0],hbond_desc[1],hbond_len[1],hbond_desc[2],hbond_len[2],hbond_desc[3],hbond_len[3],hbond_desc[4],hbond_len[4]])

df = pd.DataFrame(output_helis[1:],columns=output_helis[0])
df.to_csv('./Data/%s.csv'%oupname,index=False)


