
import os
import sys
import numpy as np
import pandas as pd

#df = pd.read_csv("./Data/final_final_MD_mimic_minorgw.csv")
df = pd.read_csv("./Data/final_final_MD_mimic_majorgw.csv")
tfs = df.tf.unique().tolist()

for idx, tf in enumerate(tfs):
    subdf = df.loc[df.tf == tf].reset_index(drop=True)

    subdf['pos1_dif'] = subdf['pos1_avg'] - subdf['pos1_avg'][0]
    subdf['pos2_dif'] = subdf['pos2_avg'] - subdf['pos2_avg'][0]
    subdf['pos3_dif'] = subdf['pos3_avg'] - subdf['pos3_avg'][0]
    subdf['pos4_dif'] = subdf['pos4_avg'] - subdf['pos4_avg'][0]

    subdf['isMimic'] = "0000"

    for i in range(2,len(subdf)):
        if (np.abs(subdf['pos1_dif'][i]) < np.abs(subdf['pos1_dif'][1])) and (np.abs(subdf['pos1_avg'][i]-subdf['pos1_avg'][1]) > subdf['pos1_std'][1]):
            subdf['isMimic'][i] = 'T' + subdf['isMimic'][i][1:4]
        if (np.abs(subdf['pos2_dif'][i]) < np.abs(subdf['pos2_dif'][1])) and (np.abs(subdf['pos2_avg'][i]-subdf['pos2_avg'][1]) > subdf['pos2_std'][1]):
            subdf['isMimic'][i] = subdf['isMimic'][i][0:1] + 'T' + subdf['isMimic'][i][2:4]
        if (np.abs(subdf['pos3_dif'][i]) < np.abs(subdf['pos3_dif'][1])) and (np.abs(subdf['pos3_avg'][i]-subdf['pos3_avg'][1]) > subdf['pos3_std'][1]):
            subdf['isMimic'][i] = subdf['isMimic'][i][0:2] + 'T' + subdf['isMimic'][i][3:4]
        if (np.abs(subdf['pos4_dif'][i]) < np.abs(subdf['pos4_dif'][1])) and (np.abs(subdf['pos4_avg'][i]-subdf['pos4_avg'][1]) > subdf['pos4_std'][1]):
            subdf['isMimic'][i] = subdf['isMimic'][i][0:3] + 'T'
        if (np.abs(subdf['pos1_dif'][i]) > np.abs(subdf['pos1_dif'][1])) and (np.abs(subdf['pos1_avg'][i]-subdf['pos1_avg'][1]) > subdf['pos1_std'][1]):
            subdf['isMimic'][i] = 'F' + subdf['isMimic'][i][1:4]
        if (np.abs(subdf['pos2_dif'][i]) > np.abs(subdf['pos2_dif'][1])) and (np.abs(subdf['pos2_avg'][i]-subdf['pos2_avg'][1]) > subdf['pos2_std'][1]):
            subdf['isMimic'][i] = subdf['isMimic'][i][0:1] + 'F' + subdf['isMimic'][i][2:4]
        if (np.abs(subdf['pos3_dif'][i]) > np.abs(subdf['pos3_dif'][1])) and (np.abs(subdf['pos3_avg'][i]-subdf['pos3_avg'][1]) > subdf['pos3_std'][1]):
            subdf['isMimic'][i] = subdf['isMimic'][i][0:2] + 'F' + subdf['isMimic'][i][3:4]
        if (np.abs(subdf['pos4_dif'][i]) > np.abs(subdf['pos4_dif'][1])) and (np.abs(subdf['pos4_avg'][i]-subdf['pos4_avg'][1]) > subdf['pos4_std'][1]):
            subdf['isMimic'][i] = subdf['isMimic'][i][0:3] + 'F'

    print(subdf[['tf_name','tf','mm_id','pos1_dif','pos2_dif','pos3_dif','pos4_dif','isMimic']])





#df = pd.read_csv("./Data/final_final_MD_mimic_beta.csv")
#df = pd.read_csv("./Data/final_final_MD_mimic_gamma.csv")
#df = pd.read_csv("./Data/final_final_MD_mimic_zeta.csv")

#tfs = df.tf.unique().tolist()
#
#for idx, tf in enumerate(tfs):
#    subdf = df.loc[df.tf == tf].reset_index(drop=True)
#
#    subdf['pos1_dif'] = subdf['pos1_avg'] - subdf['pos1_avg'][0]
#    subdf['pos2_dif'] = subdf['pos2_avg'] - subdf['pos2_avg'][0]
#    subdf['pos3_dif'] = subdf['pos3_avg'] - subdf['pos3_avg'][0]
#
#    subdf['isMimic'] = "000"
#
#    for i in range(2,len(subdf)):
#        if (np.abs(subdf['pos1_dif'][i]) < np.abs(subdf['pos1_dif'][1])) and (np.abs(subdf['pos1_avg'][i]-subdf['pos1_avg'][1]) > subdf['pos1_std'][1]):
#            subdf['isMimic'][i] = 'T' + subdf['isMimic'][i][1:4]
#        if (np.abs(subdf['pos2_dif'][i]) < np.abs(subdf['pos2_dif'][1])) and (np.abs(subdf['pos2_avg'][i]-subdf['pos2_avg'][1]) > subdf['pos2_std'][1]):
#            subdf['isMimic'][i] = subdf['isMimic'][i][0:1] + 'T' + subdf['isMimic'][i][2:4]
#        if (np.abs(subdf['pos3_dif'][i]) < np.abs(subdf['pos3_dif'][1])) and (np.abs(subdf['pos3_avg'][i]-subdf['pos3_avg'][1]) > subdf['pos3_std'][1]):
#            subdf['isMimic'][i] = subdf['isMimic'][i][0:2] + 'T' + subdf['isMimic'][i][3:4]
#        if (np.abs(subdf['pos1_dif'][i]) > np.abs(subdf['pos1_dif'][1])) and (np.abs(subdf['pos1_avg'][i]-subdf['pos1_avg'][1]) > subdf['pos1_std'][1]):
#            subdf['isMimic'][i] = 'F' + subdf['isMimic'][i][1:4]
#        if (np.abs(subdf['pos2_dif'][i]) > np.abs(subdf['pos2_dif'][1])) and (np.abs(subdf['pos2_avg'][i]-subdf['pos2_avg'][1]) > subdf['pos2_std'][1]):
#            subdf['isMimic'][i] = subdf['isMimic'][i][0:1] + 'F' + subdf['isMimic'][i][2:4]
#        if (np.abs(subdf['pos3_dif'][i]) > np.abs(subdf['pos3_dif'][1])) and (np.abs(subdf['pos3_avg'][i]-subdf['pos3_avg'][1]) > subdf['pos3_std'][1]):
#            subdf['isMimic'][i] = subdf['isMimic'][i][0:2] + 'F' + subdf['isMimic'][i][3:4]
#
#    print(subdf[['tf_name','tf','mm_id','pos1_dif','pos2_dif','pos3_dif','isMimic']])

