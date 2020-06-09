
import os
from pdblib.base import *

inpdir = './Pdb'
oupdir = './Pdb_refine'

pdblist = os.listdir(inpdir)
pdblist.sort()
print(pdblist)

seg_dic = {'1qne':[0,1,2],'6njq':[3,4,5],'xac1':[0,1,2],'xcc2':[0,1,2],'xcca':[3,4,5],'xccb':[0,1,2]}
aalist = [26,27,29,31,56,57,58,61,63,68,70,72,74,76,78,80,82,83,85,116,117,119,147,148,149,152,154,161,163,165,167,169,171,173,176]

for pdb in pdblist:
    
    pdbid = pdb[:-4]
    x,y,z = seg_dic[pdbid]

    ipdbf = os.path.join(inpdir,pdb)
    opdbf = os.path.join(oupdir,pdb)
    oaaf = os.path.join(oupdir,pdbid+"_aa.pdb")
    mol = Mol(ipdbf)
    newmol = Mol()
    newmol.segs = [Segment(),Segment(),Segment()]
    newmol.segs[0] = mol.segs[x]
    newmol.segs[1] = mol.segs[y]
    newmol.segs[2] = mol.segs[z]
    newmol.write(opdbf)

    aa = Mol()
    aa.segs = [Segment()]
    aa.segs[0].reses = [newmol.segs[0].reses[i-12] for i in aalist]
    aa.write(oaaf)



