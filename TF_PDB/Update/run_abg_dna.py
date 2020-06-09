#!/usr/bin/python

from pdblib.num import *
from supimp import *
import numpy as np
from copy import deepcopy

# define output folder
oupdir = './Ets1'
# define input pdb
inpdir = './Helice_update'
fnmol = '1dux.pdb'
#fnmol = '1k79.pdb'
# define ideal helix pdb file
fnhlx = './iBformDNA_final.pdb'

# define which atoms to align
atns = ("O3'","P","O5'","C5'","C4'","C3'","C2'","C1'","O4'")

# define region of helix
mol = Mol(os.path.join(inpdir,fnmol))
seg1 = mol.segs[0]
seg2 = mol.segs[1]

#helix 1 segment list
mseg1 = [seg1,seg1,seg1,seg2,seg2,seg2]
#mseg1 = [seg1,seg1,seg2,seg2]
#helix 1 resi list
mresi1 = [4,5,6,8,9,10]
#mresi1 = [6,7,8,9,10,11]
#mresi1 = [6,7,7,8]
#helix 1 icode list
mic1 = [' ',' ',' ',' ',' ',' ']
#mic1 = [' ',' ',' ',' ']
#helix 2 segment list
mseg2 = [seg1,seg1,seg1,seg2,seg2,seg2]
#mseg2 = [seg1,seg1,seg2,seg2]
#helix 2 resi list
mresi2 = [8,9,10,4,5,6]
#mresi2 = [10,11,12,5,6,7]
#mresi2 = [9,10,4,5]
#helix 2 icode list
mic2 = [' ',' ',' ',' ',' ',' ']
#mic2 = [' ',' ',' ',' ']

# build ref
hlx = Mol(fnhlx)
hresi1 = [8,9,10,31,32,33]
hresi2 = [11,12,13,28,29,30]
#hresi1 = [9,10,31,32]
#hresi2 = [11,12,29,30]

# helix 1
mats1 = []
hats1 = []
m1 = Mol()
m1.segs = [Segment()]
for i in range(len(mresi1)):
    mres1 = getres(mseg1[i], mresi1[i], mic1[i])
    hres1 = getres(hlx, hresi1[i])
    m1.segs[0].reses.append(mres1)
    for atn in atns:
        mat1 = mres1.getat(atn)
        hat1 = hres1.getat(atn)
        if mat1:
            mats1.append(mat1)
            hats1.append(hat1)
M1 = getmat(mats1).astype('float32')
ref1 = getmat(hats1).astype('float32')
ref1 -= mean(ref1, axis=0)

# helix 2
mats2 = []
hats2 = []
m2 = Mol()
m2.segs = [Segment()]
for i in range(len(mresi2)):
    mres2 = getres(mseg2[i], mresi2[i], mic2[i])
    hres2 = getres(hlx, hresi2[i])
    m2.segs[0].reses.append(mres2)
    for atn in atns:
        mat2 = mres2.getat(atn)
        hat2 = hres2.getat(atn)
        if mat2:
            mats2.append(mat2)
            hats2.append(hat2)
M2 = getmat(mats2).astype('float32')
ref2 = getmat(hats2).astype('float32')
ref2 -= mean(ref2, axis=0)

# align region 1 to ideal helix 1
mc_M1 = mean(M1, axis=0)
M1 -= mc_M1
rot,rmsd1 = matvec.lsqfit(M1, ref1)
# apply resulted rot matrix to region 2
M2 -= mc_M1
M2 = dot(M2, rot)
# align region 2 to ideal helix 2
M2 -= mean(M2, axis=0)
rot,rmsd2 = matvec.lsqfit(M2, ref2)
# calculate Euler angles
a = np.arctan2(rot[2,1], rot[2,0])
g = np.arctan2(rot[1,2], -rot[0,2])
b = np.arctan2(rot[1,2]/np.sin(g), rot[2,2])
# solve the degeneracy issue
A = array([a, a-pi, a+pi, a+pi, a-pi])
B = array([b, -b, -b, -b, -b])
G = array([g, g+pi, g-pi, g+pi, g-pi])
M = A**2 + B**2 + G**2
idx = np.argmin(M)
a = A[idx]
b = B[idx]
g = G[idx]

print('alpha: %8.3f'%(a*180./pi))
print('beta : %8.3f'%(b*180./pi))
print('gamma: %8.3f'%(g*180./pi))
print('rmsd1: %8.3f'%(rmsd1))
print('rmsd2: %8.3f'%(rmsd2))

#superimpose

m1.renumber(1,1)
m2.renumber(1,1)

m1.write('m1.pdb')
m2.write('m2.pdb')

m1 = Mol('m1.pdb')
m2 = Mol('m2.pdb')

hlx2 = Mol('iBformDNA_final.pdb')
#hlx2.segs[0].reses = hlx2.segs[0].reses[10:20] + hlx2.segs[0].reses[20:30]

m1.renumber(1,1)
m2.renumber(1,1)

supimpose(hlx,m1,'8-10,31-33','1-6')
supimpose(hlx2,m2,'11-13,28-30','1-6')
#supimpose(hlx,m1,'9-10,31-32','1-4')
#supimpose(hlx2,m2,'11-12,29-30','1-4')

ctbp = Mol()
ctbp.segs = [Segment(),Segment()]
ctbp.segs[0].reses = [deepcopy(hlx.segs[0].reses[9])]
ctbp.segs[1].reses = [deepcopy(hlx.segs[0].reses[30])]

mol2 = Mol(os.path.join(inpdir,fnmol))

pdb = Pdb()
hlx.mdid = 1
hlx2.mdid = 2
ctbp.mdid = 3
mol2.mdid = 4
hlx.segs[0].chid = 'A'
hlx2.segs[0].chid = 'B'
ctbp.segs[0].chid = 'C'
ctbp.segs[1].chid = 'C'
pdb.mds.append(hlx)
pdb.mds.append(hlx2)
pdb.mds.append(ctbp)
pdb.mds.append(mol2)

os.system("rm m1.pdb m2.pdb")

pdb.write(os.path.join(oupdir,fnmol[:-4]+"_align.pdb"))
