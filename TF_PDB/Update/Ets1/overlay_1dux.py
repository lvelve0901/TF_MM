
from pymol import cmd

inpdir = 'Pdb_refine'
pdb1 = '1dux'
chain1 = ['C','A','B']
helix_aa = '56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71'
sele_aa = '6,7,46,50,52,55,57,58,59,62,65,66,67,68,75,79,80,81'
bp_resi = ['6,7,8,9','4,5,6,7,8']

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

