
from pymol import cmd

inpdir = 'Pdb_refine'
pdb1 = '1k79'
chain1 = ['A','B','C']
helix_aa = '386,387,388,389,390,391,392,393,394,395,396,397'
sele_aa = '391,381,337,336,399,394,388,404,395,384,375,410,379,397,386,408,396'
bp_resi = ['8,9,10,11','5,6,7,8,9']

cmd.load('./%s/%s.pdb'%(inpdir,pdb1),object=pdb1)
cmd.create('%s_protein'%pdb1,'/%s//%s//'%(pdb1,chain1[0]))
cmd.create('%s_helix'%pdb1,'/%s//%s/%s/'%(pdb1,chain1[0],helix_aa))
cmd.create('%s_aa'%pdb1,'/%s//%s/%s/'%(pdb1,chain1[0],sele_aa))
cmd.create('%s_dna'%pdb1,'/%s//%s,%s//'%(pdb1,chain1[1],chain1[2]))
#cmd.create('%s_bp'%pdb1,'((/%s//%s/%s/) or (/%s//%s/%s/))'%(pdb1,chain1[1],bp_resi[0],pdb1,chain1[2],bp_resi[1]))
cmd.create('%s_bp1'%pdb1,'/%s//%s/%s/'%(pdb1,chain1[1],bp_resi[0]))
cmd.create('%s_bp2'%pdb1,'/%s//%s/%s/'%(pdb1,chain1[2],bp_resi[1]))
cmd.hide('line','(all)')
cmd.show('cartoon','%s_helix'%pdb1)
cmd.show('cartoon','%s_protein'%pdb1)
cmd.show('cartoon','%s_dna'%pdb1)
cmd.show('sticks','%s_aa'%pdb1)
cmd.show('sticks','%s_bp1'%pdb1)
cmd.show('sticks','%s_bp2'%pdb1)
cmd.hide('cartoon','(/%s_protein//%s/%s/)'%(pdb1,chain1[0],helix_aa))
cmd.color('magenta','%s_helix'%pdb1)
cmd.color('gray50','%s_protein'%pdb1)
cmd.color('gray50','%s_aa'%pdb1)
cmd.color('green','%s_dna'%pdb1)
cmd.color('green','%s_bp1'%pdb1)
cmd.color('green','%s_bp2'%pdb1)

cmd.center('(chain %s and resi %s)'%(chain1[1],bp_resi[0]))

cmd.bg_color(color='white')

cmd.save(filename='./%s/%s_overlay.pse'%(inpdir,pdb1),selection='(all)')

