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
oupname = 'tf_contact'

filelist = []  #list to store PDB file names
#Walk data_directory and get list of pdb names
for file in pdb_dir_list:
	if file.endswith(".pdb"):
		filelist.append(file)
filelist.sort()


#Output list
output_contact = [['pdbid','name','pair_id','bp_name','nt1_chain','nt2_chain','nt1_name','nt2_name','nt1_resi','nt2_resi','hbonds_bb_s1','hb_ats_bb_s1','hbonds_bb_s2','hb_ats_bb_s2','hbonds_sc_s1','hb_ats_sc_s1','hbonds_sc_s2','hb_ats_sc_s2','num_bb_1','num_bb_2','num_sc_1','num_sc_2']]

nt_can = ['DA', 'DG', 'DC', 'DT', '5IU', '5CM', '5HC', '5FC', '1CC', 'IGU', 'DI']
aa_can = ['ASP', 'ASN', 'ARG', 'GLU', 'ALA', 'GLN', 'SER', 'CYS', 'PRO', 'TYR', 'MET', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'PHE', 'THR', 'TRP', 'VAL']
aa_code = {'ASP':'D','ASN':'N','ARG':'R','GLU':'E','ALA':'A','GLN':'Q','SER':'S','CYS':'C','PRO':'P','TYR':'Y','MET':'M','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','PHE':'F','THR':'T','TRP':'W','VAL':'V'}
bb_atoms = ("P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'")


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
    hbonds_bb = {}
    hbonds_sc = {}
    hb_ats_bb = {}
    hb_ats_sc = {}

    hbonds = najson.json_file['hbonds']

    for i, hbond in enumerate(hbonds):

        atom1_id = hbond['atom1_id']  #atom1 id
        atom2_id = hbond['atom2_id']  #atom2 id
        distance = hbond['distance']

        # Splitting the atom_id
        # atom_id separator
        atom1 = atom1_id.split('@')[0]
        atom2 = atom2_id.split('@')[0]

        res1_id = atom1_id.split('@')[1]  #resi1 id for example "A.G17"
        res2_id = atom2_id.split('@')[1]  #resi2 id

        resn_resi1 = res1_id.split('.')[1]  #resn_resi1
        resn_resi2 = res2_id.split('.')[1]  #resn_resi2
    
        resi1 = re.split('(\-?\d+)',resn_resi1)[-2]
        resn1 = resn_resi1[:-len(resi1)]
        resi2 = re.split('(\-?\d+)',resn_resi2)[-2]
        resn2 = resn_resi2[:-len(resi2)]

        # identify whether it is protein DNA contact
        if (resn1 in nt_can and resn2 in aa_can) or (resn1 in aa_can and resn2 in nt_can):
            if resn1 in nt_can:
                nt_id = res1_id
                atom_nt = atom1
                aa_id = res2_id
                aa_name = resn2
                aa_resi = resi2
                code_aa = aa_code[resn2]
                atom_aa = atom2
            else:
                nt_id = res2_id
                atom_nt = atom2
                aa_id = res1_id
                aa_name = resn1
                aa_resi = resi1
                code_aa = aa_code[resn1]
                atom_aa = atom1
            
            if atom_nt in bb_atoms: 
                if nt_id not in hbonds_bb.keys():
                    hbonds_bb[nt_id] = []
                    hb_ats_bb[nt_id] = []
                    hbonds_bb[nt_id].append(code_aa+str(aa_resi))
                    hb_ats_bb[nt_id].append(atom_nt+"-"+atom_aa)
                else:
                    hbonds_bb[nt_id].append(code_aa+str(aa_resi))
                    hb_ats_bb[nt_id].append(atom_nt+"-"+atom_aa)
            else:
                if nt_id not in hbonds_sc.keys():
                    hbonds_sc[nt_id] = []
                    hb_ats_sc[nt_id] = []
                    hbonds_sc[nt_id].append(code_aa+str(aa_resi))
                    hb_ats_sc[nt_id].append(atom_nt+"-"+atom_aa)
                else:
                    hbonds_sc[nt_id].append(code_aa+str(aa_resi))
                    hb_ats_sc[nt_id].append(atom_nt+"-"+atom_aa)
    

    #loop the helices
    for heli in helice['pairs'][1:]:  #unpack helical parameters
        nt1_id = heli['nt1']
        nt2_id = heli['nt2']

        try:
            hbonds_bb_1 = hbonds_bb[nt1_id]
        except KeyError, e:
            hbonds_bb_1 = []

        try:
            hb_ats_bb_1 = hb_ats_bb[nt1_id]
        except KeyError, e:
            hb_ats_bb_1 = []

        try:
            hbonds_bb_2 = hbonds_bb[nt2_id]
        except KeyError, e:
            hbonds_bb_2 = []

        try:
            hb_ats_bb_2 = hb_ats_bb[nt2_id]
        except KeyError, e:
            hb_ats_bb_2 = []

        try:
            hbonds_sc_1 = hbonds_sc[nt1_id]
        except KeyError, e:
            hbonds_sc_1 = []

        try:
            hb_ats_sc_1 = hb_ats_sc[nt1_id]
        except KeyError, e:
            hb_ats_sc_1 = []

        try:
            hbonds_sc_2 = hbonds_sc[nt2_id]
        except KeyError, e:
            hbonds_sc_2 = []

        try:
            hb_ats_sc_2 = hb_ats_sc[nt2_id]
        except KeyError, e:
            hb_ats_sc_2 = []


        num_bb_1 = len(hbonds_bb_1)
        num_bb_2 = len(hbonds_bb_2)
        num_sc_1 = len(hbonds_sc_1)
        num_sc_2 = len(hbonds_sc_2)

        hbonds_bb_s1 = "&".join(hbonds_bb_1)
        hbonds_bb_s2 = "&".join(hbonds_bb_2)
        hbonds_sc_s1 = "&".join(hbonds_sc_1)
        hbonds_sc_s2 = "&".join(hbonds_sc_2)

        hb_ats_bb_s1 = "&".join(hb_ats_bb_1)
        hb_ats_bb_s2 = "&".join(hb_ats_bb_2)
        hb_ats_sc_s1 = "&".join(hb_ats_sc_1)
        hb_ats_sc_s2 = "&".join(hb_ats_sc_2)

        nt1 = nts[nts_idx[nt1_id]]
        nt2 = nts[nts_idx[nt2_id]]

        nt1_chain = nt1['chain_name']
        nt2_chain = nt2['chain_name']
        nt1_name = nt1['nt_name']
        nt2_name = nt2['nt_name']
        nt1_resi = nt1_id.split('.')[1].replace(nt1_name,"")
        nt2_resi = nt2_id.split('.')[1].replace(nt2_name,"")
        bp_name = heli['bp'].replace("-","").replace("+","")
        pair_id = pdbid + "_" + bp_name + "_" + nt1_chain + "_" + nt1_name + "_" + nt1_resi + "_" + nt2_chain + "_" + nt2_name + "_" + nt2_resi 

        if 'groove_widths' in heli.keys():
            output_contact.append([pdbid,name,pair_id,bp_name,nt1_chain,nt2_chain,nt1_name,nt2_name,nt1_resi,nt2_resi,hbonds_bb_s1,hb_ats_bb_s1,hbonds_bb_s2,hb_ats_bb_s2,hbonds_sc_s1,hb_ats_sc_s1,hbonds_sc_s2,hb_ats_sc_s2,num_bb_1,num_bb_2,num_sc_1,num_sc_2])


df = pd.DataFrame(output_contact[1:],columns=output_contact[0])
df.to_csv('./Data/%s.csv'%oupname,index=False)


