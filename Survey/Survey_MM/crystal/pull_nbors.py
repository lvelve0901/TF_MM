import pandas as pd
import os
import sys

filename = 'DNA_mismatch_GT_automatic.csv'
bf=pd.read_csv(filename)
pdbdir = '/mnt/hs189/NAfinder_2.0_crystal/Crystal/'
#command = "((chain " + str(bf.chid_1[ele]) + " and resid " + str(bf.resi_1[ele] - 2) + " " + str(bf.resi_1[ele] - 1) + " " + str(bf.resi_1[ele]    ) + " " + str(bf.resi_1[ele] + 1) + " " + str(bf.resi_1[ele] + 2) + ") or (chain  " + str(bf.chid_2[ele]) + " and resid " + str(bf.resi_2[ele] - 2) + " " + str(bf.resi_2[ele] - 1) + " " + str(bf.resi_2[ele]    ) + " " + str(bf.resi_2[ele] + 1) + " " + str(bf.resi_2[ele] + 2) + ")) "

for ele in range(0, len(bf.pdbid)):
    if bf.status_bp[ele] == "Duplex":
        script_name = str(bf.pdbid[ele]) + ".tcl"
        pdb_name =  pdbdir + str(bf.pdbid[ele]) + ".pdb"
        f = open(script_name, "w")
    
        # Write the VMD commands
        # Load the molecule
        f.write("mol new " + pdb_name + "\n")
    
        # General display commands
        f.write("display projection orthographic\n")
        f.write("display reset view\n")
        f.write("display depthcue off\n")
        f.write("display nearclip set 0.01\n")
        f.write("display height 4.0\n")
        f.write("axes location off\n")

        dna_command = "((chain " + str(bf.chid_1[ele]) + " and resid " + str(bf.resi_1[ele]) + ") or (chain  " + str(bf.chid_2[ele]) + " and resid " + str(bf.resi_2[ele]) + "))"
        dna_nbor_command = "all and same residue as within 5 of ((chain " + str(bf.chid_1[ele]) + " and resid " + str(bf.resi_1[ele]) + ") or (chain  " + str(bf.chid_2[ele]) + " and resid " + str(bf.resi_2[ele]) + "))"
    
        # Representation commands
        f.write("mol delrep 0 top\n")
        f.write("mol representation Lines 3.0\n")
        f.write("mol selection {" + str(dna_command) + "}\n")
        f.write("mol addrep top\n")
        f.write("mol representation Lines 1.0\n")
        f.write("mol selection {" + str(dna_nbor_command) + "}\n")
        f.write("mol addrep top\n")
    
        f.close()
    
        # Now execute vmd
        os.system("vmd -e " + str(script_name))

        sys.exit(0)
