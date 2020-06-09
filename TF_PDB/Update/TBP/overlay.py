
from pymol import cmd

inpdir = 'Pdb_refine'
pdb1 = '1qne'
pdb2 = 'xac1'
pdb3 = 'xcc2'
chain1 = ['A','C','D']
chain2 = ['A','B','C']
chain3 = ['B','C','D']
sele_aa = '26,27,29,31,56,57,58,61,63,68,70,72,74,76,78,80,82,83,85,116,117,119,121,147,148,149,152,154,161,163,165,167,169,171,173,174,176'
bp_resi = ['209,210','219,220']

cmd.load('./%s/%s.pdb'%(inpdir,pdb1),object=pdb1)
cmd.load('./%s/%s.pdb'%(inpdir,pdb2),object=pdb2)
cmd.load('./%s/%s.pdb'%(inpdir,pdb3),object=pdb3)
cmd.align(pdb2,pdb1,cutoff=2.0,cycles=5)
cmd.align(pdb3,pdb1,cutoff=2.0,cycles=5)
cmd.create('%s_aa'%pdb1,'/%s//%s/%s/'%(pdb1,chain1[0],sele_aa))
cmd.create('%s_aa'%pdb2,'/%s//%s/%s/'%(pdb2,chain2[0],sele_aa))
cmd.create('%s_aa'%pdb3,'/%s//%s/%s/'%(pdb3,chain3[0],sele_aa))
cmd.create('%s_dna'%pdb1,'/%s//%s,%s//'%(pdb1,chain1[1],chain1[2]))
cmd.create('%s_dna'%pdb2,'/%s//%s,%s//'%(pdb2,chain2[1],chain2[2]))
cmd.create('%s_dna'%pdb3,'/%s//%s,%s//'%(pdb3,chain3[1],chain3[2]))
cmd.create('%s_bp'%pdb1,'((/%s//%s/%s/) or (/%s//%s/%s/))'%(pdb1,chain1[1],bp_resi[0],pdb1,chain1[2],bp_resi[1]))
cmd.create('%s_bp'%pdb2,'((/%s//%s/%s/) or (/%s//%s/%s/))'%(pdb2,chain2[1],bp_resi[0],pdb2,chain2[2],bp_resi[1]))
cmd.create('%s_bp'%pdb3,'((/%s//%s/%s/) or (/%s//%s/%s/))'%(pdb3,chain3[1],bp_resi[0],pdb3,chain3[2],bp_resi[1]))
cmd.hide('line','(all)')
cmd.show('cartoon',pdb1)
cmd.show('cartoon',pdb2)
cmd.show('cartoon',pdb3)
cmd.show('sticks','%s_aa'%pdb1)
cmd.show('sticks','%s_aa'%pdb2)
cmd.show('sticks','%s_aa'%pdb3)
cmd.show('sticks','%s_bp'%pdb1)
cmd.show('sticks','%s_bp'%pdb2)
cmd.show('sticks','%s_bp'%pdb3)
cmd.color('green',pdb1)
cmd.color('orange',pdb2)
cmd.color('cyan',pdb3)
cmd.color('green',"%s_dna"%pdb1)
cmd.color('orange',"%s_dna"%pdb2)
cmd.color('cyan',"%s_dna"%pdb3)

cmd.bg_color(color='white')

cmd.save(filename='./%s/%s_overlay.pse'%(inpdir,pdb1),selection='(all)')

