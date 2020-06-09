
import json # Handle JSON DSSR files
import re
import os
import sys
import pandas as pd
import numpy as np
import learnna_json as lna_json

inpdir = './Json_DNAproDB'
filenames = os.listdir(inpdir)
filenames.sort()

nt_moiety_map = {'pp':'backbone','sr':'sugar','sg':'minorgw','wg':'majorgw','None':'base','bs':'base'}

output = [['pdbid','int_type','distance','res_pair','int_pair','int_moiety','nt_code','nt_moiety','nt_chain','nt_resn','nt_resi','aa_moiety','aa_chain','aa_resi']]

for idx, filename in enumerate(filenames):
    
    pdbid = filename.split('.')[0]
    print("Process %s [%d/%d]"%(pdbid,idx,len(filenames)))
    json_f = os.path.join(inpdir,filename)

    najson = lna_json.NA_JSON()  #initialize class objects
    with open(json_f) as json_data:  #read each json file
        data = json.load(json_data)

    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file

    chain_maps = {}
    dna = najson.json_file['dna']
    protein = najson.json_file['protein']

    for chain in dna['chains']:
        chain_id_new = chain['id']
        chain_id_old = chain['au_chain_id']
        if chain_id_new in chain_maps.keys():
            print("ERROR(): exsiting chain id: %s"%chain_id_new)
            sys.exit()
        else:
            chain_maps[chain_id_new] = chain_id_old

    for chain in protein['chains']:
        chain_id_new = chain['id']
        chain_id_old = chain['au_chain_id']
        if chain_id_new in chain_maps.keys():
            print("ERROR(): exsiting chain id: %s"%chain_id_new)
            sys.exit()
        else:
            chain_maps[chain_id_new] = chain_id_old

    interfaces = najson.json_file['interfaces']
    model = interfaces['models'][0]
    ntaa_ints = []
    for m in model:
        ntaa_ints = ntaa_ints + m['nucleotide-residue_interactions']

    for ntaa_int in ntaa_ints:

        nt_chain = chain_maps[ntaa_int['nuc_chain']]
        #nt_chain = ntaa_int['nuc_chain']
        nt_resn = ntaa_int['nuc_name']
        nt_resi = ntaa_int['nuc_number']
        nt_code = nt_chain + '.' + nt_resn + str(nt_resi)
    
        aa_chain = chain_maps[ntaa_int['res_chain']]
        #aa_chain = ntaa_int['res_chain']
        aa_resn = ntaa_int['res_name']
        aa_resi = ntaa_int['res_number']

        hbonds = ntaa_int['hbonds']
        vdws = ntaa_int['vdw_interactions']


        for hbond in hbonds:
            int_type = "hbond"
            nt_atom = hbond['nuc_atom']
            nt_moiety = nt_moiety_map[str(hbond['nuc_moiety'])]
            aa_atom = hbond['res_atom']
            aa_moiety = hbond['res_moiety']
            distance = hbond['distance']
            res_pair = nt_chain + "." + nt_resn + str(nt_resi) + ":" + aa_chain + "." + aa_resn + str(aa_resi) 
            int_pair = nt_atom + "@" + nt_chain + "." + nt_resn + str(nt_resi) + ":" + aa_atom + "@" + aa_chain + "." + aa_resn + str(aa_resi) 
            int_moiety = str(nt_moiety) + ":" + str(aa_moiety)
            output.append([pdbid,int_type,distance,res_pair,int_pair,int_moiety,nt_code,nt_moiety,nt_chain,nt_resn,nt_resi,aa_moiety,aa_chain,aa_resi])

        for vdw in vdws:
            int_type = "vdw"
            nt_atom = vdw['nuc_atom']
            nt_moiety = nt_moiety_map[str(vdw['nuc_moiety'])]
            aa_atom = vdw['res_atom']
            aa_moiety = vdw['res_moiety']
            distance = vdw['distance']
            res_pair = nt_chain + "." + nt_resn + str(nt_resi) + ":" + aa_chain + "." + aa_resn + str(aa_resi) 
            int_pair = nt_atom + "@" + nt_chain + "." + nt_resn + str(nt_resi) + ":" + aa_atom + "@" + aa_chain + "." + aa_resn + str(aa_resi) 
            int_moiety = str(nt_moiety) + ":" + str(aa_moiety)
            output.append([pdbid,int_type,distance,res_pair,int_pair,int_moiety,nt_code,nt_moiety,nt_chain,nt_resn,nt_resi,aa_moiety,aa_chain,aa_resi])

output_df = pd.DataFrame(output[1:],columns=output[0])
output_df.to_csv("./Data/tf_DNAproDB.csv",index=False)
            
