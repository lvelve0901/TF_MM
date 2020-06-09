
#!/usr/bin/python

import os
import sys
import time
import pandas as pd
from pdblib.base import *
from commontool import read, readchar


def addTER(filename1, filename2):  #Add TER in the dssr-stem file 
    os.system('cp %s %s'%(filename1,filename2))
    strings = read(filename1)  #inport stem file into string list
    chars = readchar(filename1)  #import stem file into char list
    ters = []  #line index that should add ter
    models = []  #number of MODEL line
    bpnums = [] #number of bps in each MODEL
    count = 0  #initialize count of residue id
    bpnum = -1  #initialize num of bps
    for i in range(len(strings)):
        if strings[i][0] == 'REMARK':
            bpnum = int(strings[i][2][4:])  #read num of bps
            bpnums.append(bpnum)
            models.append(strings[i-1][1])  #read model index
            count = 0  #reinitialize count of residue id
        if strings[i][0] in ['ATOM','HETATM'] and strings[i+1][0] in ['ATOM','HETATM']:
            if chars[i][17:27] != chars[i+1][17:27]:  #read column of resi
                count += 1  #count + 1 if the resi change
                if count == bpnum:  #if count of resi = num of bp
                    ters.append(i)
    for i in range(len(ters)):
        os.system("sed -i '%d i\TER' %s"%(ters[i]+i+2,filename2))
    return ters, bpnums, models






inpdir = './helice'
oupdir = './helice_refine'

os.system('mkdir -p %s'%oupdir)
curdir = os.getcwd()
stemdic = {}  #dictionary to store pdb -- stem information
pdblist = [pdbf.split('.')[0] for pdbf in os.listdir(inpdir)]  #list to store pdb file names
pdblist.sort()  #sort the pdblist

tic = time.time()  #timer: start

#if len(os.listdir(oupdir)) != 0:
#    print ("Files already in output folder! Please check.")
#    sys.exit()

old_pdblist = [pdbf.split('_')[0] for pdbf in os.listdir(oupdir)]
old_pdblist = list(set(old_pdblist))
old_pdblist.sort()


for idx, pdbid in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))
    
    if pdbid in old_pdblist:
        continue
    else:
        print("Now process the helical structure")

    pdb = pdbid + ".pdb"
    stemdic[pdbid.split('-')[0]] = ''

    #Step1: Create output dssr-stems file
    ipdbf = os.path.join(inpdir,pdb)
    tpdbf = os.path.join(curdir,pdb)
    if pdb not in os.listdir(inpdir):  #check if we already process the pdbfile
        continue

    #Step2: Add TER between two strands for each stem in stem file
    ter, bpnums, model = addTER(ipdbf,tpdbf)
    if len(ter) != len(model) or len(ter) != len(bpnums):  #check point: ter bpnums and model has same number
        print(">>>Abort! Inconsistent model, ter and bps number: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
        sys.exit()

    #Step3: Split stem file into each stem using pdblib
    stems = Pdb(tpdbf)
    for i, md in enumerate(stems.mds):
        if len(md.segs) > 2:
            print(">>>Abort! There are more than 2 chains: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
            sys.exit()
        if len(md.segs[0].reses) != len(md.segs[1].reses) or len(md.segs[0].reses) != int(bpnums[i]):
            print(">>>In model %d, seg0 length: %d seg1 length: %d bpnums: %d"%(i,len(md.segs[0].reses),len(md.segs[1].reses),int(bpnums[i])))
            print(">>>Abort! Wrong position of TER: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
            sys.exit()
        opdbf = os.path.join(oupdir,pdbid+"_"+str(i)+"_"+str(len(md.segs[0].reses))+".pdb")  #pdbid-stemid-stemlen.pdb
        stemdic[pdbid.split('-')[0]] += " "+pdbid+"_"+str(i)+"_"+str(len(md.segs[0].reses))
        md.write("%s"%opdbf)
    print(stemdic[pdbid.split('-')[0]])
    os.system("rm %s"%tpdbf)

toc = time.time()  #timer: end
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))


