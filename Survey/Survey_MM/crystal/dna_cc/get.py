import os
import sys
import pandas as pd

filename = str(sys.argv[1])
bf = pd.read_csv(filename)

# Crystal pdb path
crystal_pdbpath = "/mnt/hs189/NAfinder_2.0_crystal/Crystal/"
# get duplex
bf_duplex = bf.loc[bf['status_bp']=='Duplex']

# Get pdbs
pdblist = []
for dummy in range(len(bf_duplex['pdbid'])):
    pdbname = bf_duplex['pdbid'].iloc[dummy]
    os.system("cp " + str(crystal_pdbpath) + str(pdbname) + ".pdb .")
