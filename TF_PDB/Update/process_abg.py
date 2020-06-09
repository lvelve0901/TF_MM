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

output_abg = [['pdbid','pair_id','bp_name','nt1_chain','nt2_chain','nt1_name','nt2_name','nt1_resi','nt2_resi','c1c1_dist','rmsd1_h','rmsd2_h','alpha_h_b','beta_h_b','gamma_h_b','zeta_h_b','minorgw','majorgw']]

ntcode_df = {'DA':'A','DT':'T','DG':'G','DC':'C','5CM':'c','5HC':'c','5FC':'c','1CC':'c','IGU':'a','DI':'g','5IU':'u'}

for idx, pdbid in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))

    mol = Mol('./%s/%s.pdb'%(inpdir,pdbid))
   
    for i in range(len(mol.segs[0].reses)):
        res1 = mol.segs[0].reses[i]
        res2 = mol.segs[1].reses[-i-1]
        nt1_chain = res1.chid
        nt2_chain = res2.chid
        nt1_name = res1.name
        nt2_name = res2.name
        nt1_resi = str(res1.resi)
        nt2_resi = str(res2.resi)
        bp_name = ntcode_df[nt1_name] + ntcode_df[nt2_name]
        pair_id = pdbid + "_" + bp_name + "_" + nt1_chain + "_" + nt1_name + "_" + nt1_resi + "_" + nt2_chain + "_" + nt2_name + "_" + nt2_resi 
        
        c1c1_dist = np.nan
        alpha_h = np.nan
        beta_h = np.nan
        gamma_h = np.nan
        zeta_h = np.nan
        rmsd1_h = np.nan
        rmsd2_h = np.nan
        minorgw = np.nan
        majorgw = np.nan

        # calculation of C1'-C1' distance
        isC1C1dist = True
        atom_k_1 = None
        atom_k_2 = None
       
        try:
            res_i_1 = mol.segs[0].reses[i]
            atom_i_1 = res_i_1.getat("C1'")
        except IndexError:
            isC1C1dist = False
        try:
            res_i_2 = mol.segs[1].reses[len(mol.segs[1].reses)-i-1]
            atom_i_2 = res_i_2.getat("C1'")
        except IndexError:
            isC1C1dist = False

        if isC1C1dist == True:
            c1c1_dist = round(dist(atom_i_1,atom_i_2),3)

        # calculation of ABG angle with 2bp helix
        #hlen = 2
        #if i > hlen-1 and i < len(mol.segs[0].reses) - hlen:
        #    mresl1 = mol.segs[0].reses[i-2:i-2+2] + mol.segs[1].reses[-i:][0:2]
        #    mresl2 = mol.segs[0].reses[i-2+2+1:i-2+2+3] + mol.segs[1].reses[-i-3:-i-3+2]
        #    rmsd1_h, rmsd2_h, alpha_h, beta_h, gamma_h = getabgB2(bhlx,mresl1,mresl2)  
        #    
        #    if beta_h < 0 and beta_h != np.nan:
        #        beta_h = -beta_h
        #        if alpha_h >= -180 and alpha_h < 0:
        #            alpha_h = alpha_h + 180
        #        else:
        #            alpha_h = alpha_h - 180
        #        if gamma_h >= -180 and gamma_h < 0:
        #            gamma_h = gamma_h + 180
        #        else:
        #            gamma_h = gamma_h - 180
        #    if beta_h != np.nan:
        #        zeta_h = alpha_h + gamma_h
        #        if zeta_h < -180:
        #            zeta_h = zeta_h + 360

        # calculation of ABG angle with 3bp helix
        hlen = 3
        if i > hlen-1 and i < len(mol.segs[0].reses) - hlen:
            mresl1 = mol.segs[0].reses[i-hlen:i-hlen+hlen] + mol.segs[1].reses[-i:][0:hlen]
            mresl2 = mol.segs[0].reses[i-hlen+hlen+1:i-hlen+hlen+hlen+1] + mol.segs[1].reses[-i-hlen-1:-i-hlen-1+hlen]
            rmsd1_h, rmsd2_h, alpha_h, beta_h, gamma_h = getabgB2(bhlx,mresl1,mresl2)  
            
            if beta_h < 0 and beta_h != np.nan:
                beta_h = -beta_h
                if alpha_h >= -180 and alpha_h < 0:
                    alpha_h = alpha_h + 180
                else:
                    alpha_h = alpha_h - 180
                if gamma_h >= -180 and gamma_h < 0:
                    gamma_h = gamma_h + 180
                else:
                    gamma_h = gamma_h - 180
            if beta_h != np.nan:
                zeta_h = alpha_h + gamma_h
                if zeta_h < -180:
                    zeta_h = zeta_h + 360

       
        # calculation of minor groove width
        isMinorgw = True
        atom_i_1 = None
        atom_i_2 = None
        atom_j_1 = None
        atom_j_2 = None

        isMajorgw = True
        atom_k_1 = None
        atom_k_2 = None
       
        try:
            res_i_1 = mol.segs[0].reses[i+2]
            atom_i_1 = res_i_1.getat("P")
        except IndexError:
            isMinorgw = False
        try:
            res_i_2 = mol.segs[0].reses[i+3]
            atom_i_2 = res_i_2.getat("P")
        except IndexError:
            isMinorgw = False
        try:
            res_j_1 = mol.segs[1].reses[len(mol.segs[1].reses)-i+1]
            atom_j_1 = res_j_1.getat("P")
        except IndexError:
            isMinorgw = False
        try:
            res_j_2 = mol.segs[1].reses[len(mol.segs[1].reses)-i]
            atom_j_2 = res_j_2.getat("P")
        except IndexError:
            isMinorgw = False

        if i > 0 and i < len(mol.segs[0].reses)-2:
            try:
                res_k_1 = mol.segs[0].reses[i-1]
                atom_k_1 = res_k_1.getat("P")
            except IndexError:
                isMajorgw = False
            try:
                res_k_2 = mol.segs[1].reses[len(mol.segs[1].reses)-i-3]
                atom_k_2 = res_k_2.getat("P")
            except IndexError:
                isMajorgw = False
        
        if (atom_i_1 is None) or (atom_i_2 is None) or (atom_j_1 is None) or (atom_j_2 is None):
            isMinorgw = False

        if (atom_k_1 is None) or (atom_k_2 is None):
            isMajorgw = False

        if isMinorgw == True:
            minorgw = round((dist(atom_i_1,atom_j_1) + dist(atom_i_2,atom_j_2)) / 2.,3)
        elif isMinorgw == False:
            minorgw = np.nan

        if isMajorgw == True:
            majorgw = round(dist(atom_k_1,atom_k_2),3)
        elif isMajorgw == False:
            majorgw = np.nan

        output_abg.append([pdbid,pair_id,bp_name,nt1_chain,nt2_chain,nt1_name,nt2_name,nt1_resi,nt2_resi,c1c1_dist,rmsd1_h,rmsd2_h,alpha_h,beta_h,gamma_h,zeta_h,minorgw,majorgw])

df = pd.DataFrame(output_abg[1:],columns=output_abg[0])
df.to_csv('./Data/tf_abg.csv',index=False)
        
