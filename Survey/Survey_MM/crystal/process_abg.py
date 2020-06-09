#!/usr/bin/python

import os
import sys
import pylab as pl
import numpy as np
import pandas as pd
from commontool import read
from pdblib.abg import *

inpdir = './helice_refine'
bfnhlx = './iBformDNA_final.pdb'
bhlx = Mol(bfnhlx)

pdblist = [pdbf[:-4] for pdbf in os.listdir(inpdir)]
pdblist.sort()

output_abg = [['pdbf','pair_id','bp_name','nt1_chain','nt2_chain','nt1_name','nt2_name','nt1_resi','nt2_resi','rmsd1_h','rmsd2_h','alpha_h_b','beta_h_b','gamma_h_b','zeta_h_b','minorgw']]



for idx, pdbf in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbf,idx+1,len(pdblist)))
    pdbid = pdbf.split('_')[0]

    mol = Mol('./%s/%s.pdb'%(inpdir,pdbf))
   
    for i in range(len(mol.segs[0].reses)):
        res1 = mol.segs[0].reses[i]
        res2 = mol.segs[1].reses[-i-1]
        nt1_chain = res1.chid
        nt2_chain = res2.chid
        nt1_name = res1.name
        nt2_name = res2.name
        nt1_resi = str(res1.resi)
        nt2_resi = str(res2.resi)
        bp_name = nt1_name[-1] + nt2_name[-1]
        pair_id = pdbid + "_" + bp_name + "_" + nt1_chain + "_" + nt1_name + "_" + nt1_resi + "_" + "^" + "_" + nt2_chain + "_" + nt2_name + "_" + nt2_resi + "_" + "^"
        
        # calculation of ABG angle
        alpha_h = np.nan
        beta_h = np.nan
        gamma_h = np.nan
        zeta_h = np.nan
        rmsd1_h = np.nan
        rmsd2_h = np.nan
        if i > 1 and i < len(mol.segs[0].reses) - 2:
            mresl1 = mol.segs[0].reses[i-2:i-2+2] + mol.segs[1].reses[-i:][0:2]
            mresl2 = mol.segs[0].reses[i-2+2+1:i-2+2+3] + mol.segs[1].reses[-i-3:-i-3+2]
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
        
        
        if (atom_i_1 is None) or (atom_i_2 is None) or (atom_j_1 is None) or (atom_j_2 is None):
            isMinorgw = False

        if isMinorgw == True:
            minorgw = round((dist(atom_i_1,atom_j_1) + dist(atom_i_2,atom_j_2)) / 2.,3)
        elif isMinorgw == False:
            minorgw = np.nan

        output_abg.append([pdbid,pair_id,bp_name,nt1_chain,nt2_chain,nt1_name,nt2_name,nt1_resi,nt2_resi,rmsd1_h,rmsd2_h,alpha_h,beta_h,gamma_h,zeta_h,minorgw])

df = pd.DataFrame(output_abg[1:],columns=output_abg[0])
df.to_csv('./pdbinfo_abg.csv',index=False)
        
