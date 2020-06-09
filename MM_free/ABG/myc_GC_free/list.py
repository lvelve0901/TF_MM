#!/usr/bin/python

import numpy as np
from ABG import *
from commontool import read

log = read('log.txt')
finalist = []

for line in log:
    pdbname = line[0][:-4]
    pdbid = pdbname.split('-')[-1]  
    num = line[1].split(',')
    resi = str(int(num[0].split(':')[1])+1)
    RMSD1 = float(line[6])
    RMSD2 = float(line[8])
    alpha = float(line[10])
    beta = float(line[11])
    gamma = float(line[12])
    zeta = alpha + gamma
    if zeta > 180:
        zeta = zeta - 360
    elif zeta < -180:
        zeta = zeta + 360
    if beta >= 0:
    	if gamma > -90 and gamma < 90:
            finalist.append([pdbname, pdbid, resi, RMSD1, RMSD2, alpha, beta, gamma, zeta, 'Major'])
    	else:
            finalist.append([pdbname, pdbid, resi, RMSD1, RMSD2, alpha, beta, gamma, zeta, 'Minor'])
    if beta < 0:
    	beta = -beta
    	alpha = alpha + 180
    	gamma = gamma + 180
    	if alpha > 180:
    		alpha = alpha - 360
    	if gamma > 180:
    		gamma = gamma - 360
    	zeta = alpha + gamma
        if zeta > 180:
            zeta = zeta - 360
        elif zeta < -180:
            zeta = zeta + 360
        if gamma >= -90 and gamma <= 90:
    		finalist.append([pdbname, pdbid, resi, RMSD1, RMSD2, alpha, beta, gamma, zeta, 'Major'])
    	else:
    		finalist.append([pdbname, pdbid, resi, RMSD1, RMSD2, alpha, beta, gamma, zeta, 'Minor'])

file = open('abglist.txt','wt')
file.write('          pdbname pdbid resi RMSD1 RMSD2   alpha    beta   gamma    zeta  direction  \n')
for final in finalist:
	file.write('%17s %5s %4s %5.2f %5.2f  %6.1f  %6.1f  %6.1f  %6.1f  %9s\n'%(final[0], final[1], final[2], final[3], final[4], final[5], final[6], final[7], final[8],final[9]))

file.close()

