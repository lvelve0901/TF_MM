
import os
from pdblib.base import *

inpdir = './Pdb'
oupdir = './Pdb_refine'

pdblist = os.listdir(inpdir)
pdblist.sort()
print(pdblist)

seg_dic = {'1k79':[0,1,4],'1dux':[0,1,4]}

for pdb in pdblist:
    
    pdbid = pdb[:-4]
    x,y,z = seg_dic[pdbid]

    ipdbf = os.path.join(inpdir,pdb)
    opdbf = os.path.join(oupdir,pdb)
    mol = Mol(ipdbf)
    newmol = Mol()
    newmol.segs = [Segment(),Segment(),Segment()]
    newmol.segs[0] = mol.segs[x]
    newmol.segs[1] = mol.segs[y]
    newmol.segs[2] = mol.segs[z]
    newmol.write(opdbf)


