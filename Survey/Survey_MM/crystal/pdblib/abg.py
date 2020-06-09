from pdblib.num import *

#===============================================================================
class Mol(Mol):
    pass

### Public Functions ###
#===============================================================================
def getabgA1(hlx,mol,mresi1,mic1,mresi2,mic2):

    # build ref
    hresi1 = Ah1dic[str(len(mresi1)/2)]
    hresi2 = Ah2dic[str(len(mresi2)/2)]
    
    # helix 1
    mats1 = []
    hats1 = []
    for i in range(len(mresi1)):
        mres1 = getres(mol, mresi1[i], mic1[i])
        hres1 = getres(hlx, hresi1[i])
        for atn in atns:
            mat1 = mres1.getat(atn)
            hat1 = hres1.getat(atn)
            if mat1:
                mats1.append(mat1)
                hats1.append(hat1)
    M1 = getmat(mats1).astype('float32')
    ref1 = getmat(hats1).astype('float32')
    ref1 -= mean(ref1, axis=0)
    
    # helix 2
    mats2 = []
    hats2 = []
    for i in range(len(mresi2)):
        mres2 = getres(mol, mresi2[i], mic2[i])
        hres2 = getres(hlx, hresi2[i])
        for atn in atns:
            mat2 = mres2.getat(atn)
            hat2 = hres2.getat(atn)
            if mat2:
                mats2.append(mat2)
                hats2.append(hat2)
    M2 = getmat(mats2).astype('float32')
    ref2 = getmat(hats2).astype('float32')
    ref2 -= mean(ref2, axis=0)
    
    # align region 1 to ideal helix 1
    mc_M1 = mean(M1, axis=0)
    M1 -= mc_M1
    rot,rmsd1 = matvec.lsqfit(M1, ref1)
    # apply resulted rot matrix to region 2
    M2 -= mc_M1
    M2 = dot(M2, rot)
    # align region 2 to ideal helix 2
    M2 -= mean(M2, axis=0)
    rot,rmsd2 = matvec.lsqfit(M2, ref2)
    # calculate Euler angles
    a = np.arctan2(rot[2,1], rot[2,0])
    g = np.arctan2(rot[1,2], -rot[0,2])
    b = np.arctan2(rot[1,2]/np.sin(g), rot[2,2])
    # solve the degeneracy issue
    A = array([a, a-pi, a+pi, a+pi, a-pi])
    B = array([b, -b, -b, -b, -b])
    G = array([g, g+pi, g-pi, g+pi, g-pi])
    M = A**2 + B**2 + G**2
    idx = np.argmin(M)
    a = A[idx]
    b = B[idx]
    g = G[idx]
    
    return round(rmsd1,2), round(rmsd2,2), round(-a*180./pi,1), round(b*180./pi,1), round(-g*180./pi,1)

def getabgA2(hlx,mresl1,mresl2):

    # build ref
    hresi1 = Ah1dic[str(len(mresl1)/2)]
    hresi2 = Ah2dic[str(len(mresl2)/2)]
    
    # helix 1
    mats1 = []
    hats1 = []
    for i in range(len(mresl1)):
        mres1 = mresl1[i]
        hres1 = getres(hlx, hresi1[i])
        for atn in atns:
            mat1 = mres1.getat(atn)
            hat1 = hres1.getat(atn)
            if mat1:
                mats1.append(mat1)
                hats1.append(hat1)
    M1 = getmat(mats1).astype('float32')
    ref1 = getmat(hats1).astype('float32')
    ref1 -= mean(ref1, axis=0)
    
    # helix 2
    mats2 = []
    hats2 = []
    for i in range(len(mresl2)):
        mres2 = mresl2[i]
        hres2 = getres(hlx, hresi2[i])
        for atn in atns:
            mat2 = mres2.getat(atn)
            hat2 = hres2.getat(atn)
            if mat2:
                mats2.append(mat2)
                hats2.append(hat2)
    M2 = getmat(mats2).astype('float32')
    ref2 = getmat(hats2).astype('float32')
    ref2 -= mean(ref2, axis=0)
    
    # align region 1 to ideal helix 1
    mc_M1 = mean(M1, axis=0)
    M1 -= mc_M1
    rot,rmsd1 = matvec.lsqfit(M1, ref1)
    # apply resulted rot matrix to region 2
    M2 -= mc_M1
    M2 = dot(M2, rot)
    # align region 2 to ideal helix 2
    M2 -= mean(M2, axis=0)
    rot,rmsd2 = matvec.lsqfit(M2, ref2)
    # calculate Euler angles
    a = np.arctan2(rot[2,1], rot[2,0])
    g = np.arctan2(rot[1,2], -rot[0,2])
    b = np.arctan2(rot[1,2]/np.sin(g), rot[2,2])
    # solve the degeneracy issue
    A = array([a, a-pi, a+pi, a+pi, a-pi])
    B = array([b, -b, -b, -b, -b])
    G = array([g, g+pi, g-pi, g+pi, g-pi])
    M = A**2 + B**2 + G**2
    idx = np.argmin(M)
    a = A[idx]
    b = B[idx]
    g = G[idx]
    
    return round(rmsd1,2), round(rmsd2,2), round(-a*180./pi,1), round(b*180./pi,1), round(-g*180./pi,1)

def getabgB1(hlx,mol,mresi1,mic1,mresi2,mic2):

    # build ref
    hresi1 = Bh1dic[str(len(mresi1)/2)]
    hresi2 = Bh2dic[str(len(mresi2)/2)]
    
    # helix 1
    mats1 = []
    hats1 = []
    for i in range(len(mresi1)):
        mres1 = getres(mol, mresi1[i], mic1[i])
        hres1 = getres(hlx, hresi1[i])
        for atn in atns:
            mat1 = mres1.getat(atn)
            hat1 = hres1.getat(atn)
            if mat1:
                mats1.append(mat1)
                hats1.append(hat1)
    M1 = getmat(mats1).astype('float32')
    ref1 = getmat(hats1).astype('float32')
    ref1 -= mean(ref1, axis=0)
    
    # helix 2
    mats2 = []
    hats2 = []
    for i in range(len(mresi2)):
        mres2 = getres(mol, mresi2[i], mic2[i])
        hres2 = getres(hlx, hresi2[i])
        for atn in atns:
            mat2 = mres2.getat(atn)
            hat2 = hres2.getat(atn)
            if mat2:
                mats2.append(mat2)
                hats2.append(hat2)
    M2 = getmat(mats2).astype('float32')
    ref2 = getmat(hats2).astype('float32')
    ref2 -= mean(ref2, axis=0)
    
    # align region 1 to ideal helix 1
    mc_M1 = mean(M1, axis=0)
    M1 -= mc_M1
    rot,rmsd1 = matvec.lsqfit(M1, ref1)
    # apply resulted rot matrix to region 2
    M2 -= mc_M1
    M2 = dot(M2, rot)
    # align region 2 to ideal helix 2
    M2 -= mean(M2, axis=0)
    rot,rmsd2 = matvec.lsqfit(M2, ref2)
    # calculate Euler angles
    a = np.arctan2(rot[2,1], rot[2,0])
    g = np.arctan2(rot[1,2], -rot[0,2])
    b = np.arctan2(rot[1,2]/np.sin(g), rot[2,2])
    # solve the degeneracy issue
    A = array([a, a-pi, a+pi, a+pi, a-pi])
    B = array([b, -b, -b, -b, -b])
    G = array([g, g+pi, g-pi, g+pi, g-pi])
    M = A**2 + B**2 + G**2
    idx = np.argmin(M)
    a = A[idx]
    b = B[idx]
    g = G[idx]
    
    return round(rmsd1,2), round(rmsd2,2), round(a*180./pi,1), round(b*180./pi,1), round(g*180./pi,1)

def getabgB2(hlx,mresl1,mresl2):

    # build ref
    hresi1 = Bh1dic[str(len(mresl1)/2)]
    hresi2 = Bh2dic[str(len(mresl2)/2)]
    
    # helix 1
    mats1 = []
    hats1 = []
    for i in range(len(mresl1)):
        mres1 = mresl1[i]
        hres1 = getres(hlx, hresi1[i])
        for atn in atns:
            mat1 = mres1.getat(atn)
            hat1 = hres1.getat(atn)
            if mat1 and hat1:
                mats1.append(mat1)
                hats1.append(hat1)
    M1 = getmat(mats1).astype('float32')
    ref1 = getmat(hats1).astype('float32')
    ref1 -= mean(ref1, axis=0)
    
    # helix 2
    mats2 = []
    hats2 = []
    for i in range(len(mresl2)):
        mres2 = mresl2[i]
        hres2 = getres(hlx, hresi2[i])
        for atn in atns:
            mat2 = mres2.getat(atn)
            hat2 = hres2.getat(atn)
            if mat2 and hat2:
                mats2.append(mat2)
                hats2.append(hat2)
    M2 = getmat(mats2).astype('float32')
    ref2 = getmat(hats2).astype('float32')
    ref2 -= mean(ref2, axis=0)
    
    # align region 1 to ideal helix 1
    mc_M1 = mean(M1, axis=0)
    M1 -= mc_M1
    rot,rmsd1 = matvec.lsqfit(M1, ref1)
    # apply resulted rot matrix to region 2
    M2 -= mc_M1
    M2 = dot(M2, rot)
    # align region 2 to ideal helix 2
    M2 -= mean(M2, axis=0)
    rot,rmsd2 = matvec.lsqfit(M2, ref2)
    # calculate Euler angles
    a = np.arctan2(rot[2,1], rot[2,0])
    g = np.arctan2(rot[1,2], -rot[0,2])
    b = np.arctan2(rot[1,2]/np.sin(g), rot[2,2])
    # solve the degeneracy issue
    A = array([a, a-pi, a+pi, a+pi, a-pi])
    B = array([b, -b, -b, -b, -b])
    G = array([g, g+pi, g-pi, g+pi, g-pi])
    M = A**2 + B**2 + G**2
    idx = np.argmin(M)
    a = A[idx]
    b = B[idx]
    g = G[idx]
    
    return round(rmsd1,2), round(rmsd2,2), round(a*180./pi,1), round(b*180./pi,1), round(g*180./pi,1)


### Global Variable ###
#===============================================================================

# define which atoms to align
atns = ("O3'","P","O5'","C5'","C4'","C3'","C2'","C1'","O4'","O2'")

# define region of reference helix
Ah1dic = {'2':[10,11,34,35],'3':[9,10,11,34,35,36],'4':[8,9,10,11,34,35,36,37]}
Ah2dic = {'2':[12,13,32,33],'3':[12,13,14,31,32,33],'4':[12,13,14,15,30,31,32,33]}
Bh1dic = {'2':[9,10,31,32],'3':[8,9,10,31,32,33],'4':[7,8,9,10,31,32,33,34]}
Bh2dic = {'2':[11,12,29,30],'3':[11,12,13,28,29,30],'4':[11,12,13,14,27,28,29,30]}

