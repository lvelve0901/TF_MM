#!/usr/bin/python

import json # Handle JSON DSSR files
import re
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import learnna_json as lna_json
from pdblib.base import *
from commontool import read

mismatches = ['aa','ac','ag','cc','ct','tt','gg','gt']
pdblist = []

os.system("mkdir -p raw_pdb")
os.system("mkdir -p helice")

for mismatch in mismatches:

    df = pd.read_csv("./dna_%s/DNA_mismatch_%s_nmr_automatic.csv"%(mismatch,mismatch.upper())) 
    sub_pdblist = df.pdbid.unique().tolist()
    pdblist = pdblist + sub_pdblist

pdblist = list(set(pdblist))
pdblist.sort()
print(pdblist)

for pdbid in pdblist:
    os.system("cp /mnt/hs189/NAfinder_2.0_nmr/Nmr/%s.pdb ./raw_pdb/%s.pdb"%(pdbid,pdbid))


