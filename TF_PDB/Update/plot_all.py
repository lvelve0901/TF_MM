
import os, sys
import pandas as pd
import numpy as np
import pylab as pl


pdbinfo = pd.read_csv("./Data/refine_pdbinfo_final.csv")
namelist = pdbinfo.name.unique().tolist()

#namelist = ['TBP']
#namelist = ['Ets1']

for idx, name in enumerate(namelist):
    print("<<< Working on %s [%d/%d] >>>"%(name,idx,len(namelist)))
    #os.system("python plot_bp.py %s"%name)
    os.system("python plot_abg.py %s"%name)
    #os.system("python plot_minorgw.py %s"%name)
    #os.system("python plot_step.py %s"%name)
    #os.system("python plot_contact.py %s"%name)

