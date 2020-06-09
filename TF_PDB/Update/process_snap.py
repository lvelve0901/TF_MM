#!/usr/bin/python

import json # Handle JSON DSSR files
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import learnna_json as lna_json
from pprint import pprint
from pdblib.base import *

json_dir = './Json_snap'
os.system("mkdir -p %s"%json_dir)
pdb_dir = './RawPdb'
pdb_dir_list = os.listdir(pdb_dir)
pdbinfo = pd.read_csv("./Data/refine_pdbinfo_final.csv")
oupname = 'tf_snap'

filelist = []  #list to store PDB file names
#Walk data_directory and get list of pdb names
for file in pdb_dir_list:
	if file.endswith(".pdb"):
		filelist.append(file)
filelist.sort()


#Output list
output_snap_hbonds = [['pdbid','name','hbond_type','distance','nt_atom_id','aa_atom_id','nt_chain','nt_resn','nt_resi','nt_atom','aa_chain','aa_resn','aa_resi','aa_atom']]
output_snap_stacks = [['pdbid','name','plane_angle','vertical_dist','nt_id','nt_chain','nt_resn','nt_resi','aa_id','aa_chain','aa_resn','aa_resi']]

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


    #create a nts dictionary with key of nt_code
    nts_idx = {}  #nts index dictionary with nt_id as key
    nts = najson.json_file['nts']
    for i, nt in enumerate(nts):
        nts_idx[nt['nt_id']] = i
    if len(nts_idx) != len(nts):  #error: there are nucleotides labeled same
        print "ERROR(): Potential same naming of nucleotide"
        sys.exit()

    try:
        phosAA_hbonds = najson.json_file['phosAA_hbonds']
    except KeyError:
        phosAA_hbonds = []

    try:
        sugarAA_hbonds = najson.json_file['sugarAA_hbonds']
    except KeyError:
        sugarAA_hbonds = []

    try:
        baseAA_hbonds = najson.json_file['baseAA_hbonds']
    except KeyError:
        baseAA_hbonds = []

    try:
        ntAA_stacks = najson.json_file['ntAA_stacks']
    except KeyError:
        ntAA_stacks = []

    ntAA_hbonds = phosAA_hbonds + sugarAA_hbonds + baseAA_hbonds

    for i, hbond in enumerate(ntAA_hbonds):

        nt_atom_id = hbond['nt_atom']  #atom1 id
        aa_atom_id = hbond['aa_atom']  #atom2 id
        distance = hbond['dist']
        hbond_type = hbond['type']

        # Splitting the atom_id
        # atom_id separator
        nt_atom = nt_atom_id.split('@')[0]
        aa_atom = aa_atom_id.split('@')[0]

        nt_id = nt_atom_id.split('@')[1]  #resi1 id for example "A.G17"
        aa_id = aa_atom_id.split('@')[1]  #resi2 id

        nt_resn_resi = nt_id.split('.')[1]  #resn_resi1
        aa_resn_resi = aa_id.split('.')[1]  #resn_resi2

        if ":" in nt_id.split('.')[0]:
            nt_chain = nt_id.split('.')[0].split(':')[1]
            aa_chain = aa_id.split('.')[0].split(':')[1]
        else:
            nt_chain = nt_id.split('.')[0]
            aa_chain = aa_id.split('.')[0]
    
        nt_resi = re.split('(\-?\d+)',nt_resn_resi)[-2]
        nt_resn = nt_resn_resi[:-len(nt_resi)]
        aa_resi = re.split('(\-?\d+)',aa_resn_resi)[-2]
        aa_resn = aa_resn_resi[:-len(aa_resi)]

        nt = nts[nts_idx[nt_id]]

        output_snap_hbonds.append([pdbid,name,hbond_type,distance,nt_atom_id,aa_atom_id,nt_chain,nt_resn,nt_resi,nt_atom,aa_chain,aa_resn,aa_resi,aa_atom])

    for i, stack in enumerate(ntAA_stacks):

        nt_id = stack['nt']  #atom1 id
        aa_id = stack['aa']  #atom2 id
        plane_angle = stack['plane_angle']
        vertical_dist = stack['vertical_dist']

        nt_resn_resi = nt_id.split('.')[1]  #resn_resi1
        aa_resn_resi = aa_id.split('.')[1]  #resn_resi2

        if ":" in nt_id.split('.')[0]:
            nt_chain = nt_id.split('.')[0].split(':')[1]
            aa_chain = aa_id.split('.')[0].split(':')[1]
        else:
            nt_chain = nt_id.split('.')[0]
            aa_chain = aa_id.split('.')[0]
    
        nt_resi = re.split('(\-?\d+)',nt_resn_resi)[-2]
        nt_resn = nt_resn_resi[:-len(nt_resi)]
        aa_resi = re.split('(\-?\d+)',aa_resn_resi)[-2]
        aa_resn = aa_resn_resi[:-len(aa_resi)]

        nt = nts[nts_idx[nt_id]]

        output_snap_stacks.append([pdbid,name,plane_angle,vertical_dist,nt_id,nt_chain,nt_resn,nt_resi,aa_id,aa_chain,aa_resn,aa_resi])

df_snap_hbonds = pd.DataFrame(output_snap_hbonds[1:],columns=output_snap_hbonds[0])
df_snap_stacks = pd.DataFrame(output_snap_stacks[1:],columns=output_snap_stacks[0])
df_snap_hbonds.to_csv('./Data/%s_hbonds.csv'%oupname)
df_snap_stacks.to_csv('./Data/%s_stacks.csv'%oupname)

