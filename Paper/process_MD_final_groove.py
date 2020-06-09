
from Bio.Seq import complement
import re
import sys
import pandas as pd
import pylab as pl
import numpy as np
from commontool import read
from pprint import pprint

pd.options.display.float_format = '{:,.2f}'.format



tf_dic = {
    
    'myc_max'    : ['myc_max','specific','CTCTGCCACGTGGGTCGT',8,'A-T','bound/free/free/free','A-T/A-T/G-T/C-T','//+/-','AT_bound/AT_free/GT_free/CT_free'],
    'max_max'    : ['max_max','specific','GTAGATCACGTGAATGAA',8,'A-T','bound/free/free/free','A-T/A-T/G-T/C-T','//+/-','AT_bound/AT_free/GT_free/CT_free'], 
    'mad_max'    : ['mad_max','specific','GTAGATCACGTGAATGAA',8,'A-T','bound/free/free/free','A-T/A-T/G-T/C-T','//+/-','AT_bound/AT_free/GT_free/CT_free'],
    'atf1'       : ['atf1','specific','GTGCCAACGTCACCTTCT',7,'A-T','bound/free/free/free','A-T/A-T/G-T/C-T','//+/-','AT_bound/AT_free/GT_free/CT_free'],
    'ctcf-1'     : ['ctcf','specific','CAGCGCCCTCTACTGGC',6,'C-G','bound/free/free/free/free','C-G/C-G/T-G/G*-G/G-G*','//+/-/-','GC_bound/GC_free/GT_free/GG1_free/GG2_free'],
    'runx1'      : ['runx1','specific','TCTTATGTGGTTTTTG',5,'A-T','bound/free/free/free/free','A-T/A-T/T-T/C-T/G-T','//+/+/-','AT_bound/AT_free/TT_free/CT_free/GT_free'],
    'egr1'       : ['egr1','specific','GTGGCGTGGGCGATA',11,'C-G','bound/free/free/free/free','C-G/C-G/G-G*/G*-G/G-T','//+/+/-','GC_bound/GC_free/GG1_free/GG2_free/GT_free'],
    'rela'       : ['rela','specific','CCCGGAAATTCCCC',4,'G-C','bound/free/free/free','G-C/G-C/G-G*/G*-G','//+/+','GC_bound/GC_free/GG1_free/GG2_free'],
    'stat3-1'    : ['stat3','specific','AAAGTTCCTGGAATTTC',9,'T-A','bound/free/free','T-A/T-A/C-A','//+','AT_bound/AT_free/AC_free'],
    'gr'         : ['gr','specific','GTCAGAACATGATGTTCTCAAA',6,'A-T','bound/free/free/free/free','A-T/A-T/G-T/T-T/C-T','//+/+/-','AT_bound/AT_free/GT_free/TT_free/CT_free'], 
    'p53-1'      : ['p53','specific','AGACATGCCCGGGCATGCCT',5,'A-T','bound/free/free/free/free','A*-T/A-T/C-T/T-T/G-T','//+/+/-','HG_bound/WC_free/CT_free/TT_free/GT_free'],
    'elk1'       : ['elk1','specific','ATACCGGAACAGGG',7,'G-C','bound/free/free/free/free/free','G-C/G-C/G-A/G-A*/G*-G/G-G*','//+/+/-/-','GC_bound/GC_free/GA2_free/GA1_free/GG1_free/GG2_free'],
    'ets1-5'     : ['ets1','specific','ATACCGGAACAGG',7,'G-C','bound/free/free/free/free/free/free','G-C/G-C/G-A/G-A*/G-T/G*-G/G-G*','//+/+/+/-/-','GC_bound/GC_free/GA2_free/GA1_free/GT_free/GG1_free/GG2_free'],
    'ets1-2'     : ['ets1','non-specific','GAGGGGTAAGCGG',7,'T-A','bound/free/free/free','T-A/T-A/G-A/G-A*','//+/+','AT_bound/AT_free/GA2_free/GA1_free'],
    'tbp'        : ['tbp','specific','GCTATAAAAGGGCA',10,'G-C','bound/free/free','G-C/G-C/C-C','//+','AdMLP_bound/AdMLP_free/CC2_free']

    }

tf_list = ['p53-1','ets1-5','ets1-2','elk1','myc_max','max_max','mad_max','atf1','runx1','stat3-1','rela','egr1','ctcf-1','gr','tbp']
params = ['minorgw','majorgw']


for i, param in enumerate(params):

    print("--- Working on parameters: %s [%d/%d] ---"%(param,i,len(params)))
    output = [['tf_name','tf','seq_type','core','strand1','strand2','position','sim_type','mismatch','mm_id','binding','pos1_avg','pos1_std','pos2_avg','pos2_std','pos3_avg','pos3_std','pos4_avg','pos4_std']]

    for j, tf in enumerate(tf_list):
        
        print(">>> Processing protein: %s [%d/%d] >>>"%(tf,j,len(tf_list)))

        tf_name = tf_dic[tf][0]
        seq_type = tf_dic[tf][1]
        core = tf_dic[tf][2]
        position = tf_dic[tf][3]
        wt_bp = tf_dic[tf][4]
        sim_types = tf_dic[tf][5].split('/')
        mismatches = tf_dic[tf][6].split('/')
        bindings = tf_dic[tf][7].split('/')
        folders = tf_dic[tf][8].split('/')

        core = core[0:position] + "*" + core[position:]
        nt1_list = range(position-2,position+2)

        df_temp = pd.DataFrame([],columns=['nt1_i','value'])
        df_temp['nt1_i'] = pd.Series(nt1_list)

        for i in range(len(mismatches)):

            sim_type = sim_types[i]
            mismatch = mismatches[i]
            binding = bindings[i]
            folder = folders[i]
            mm_id = folder.split('_')[0]
            base1 = mismatch.split('-')[0].replace("*","")
            base2 = mismatch.split('-')[1].replace("*","")

            strand1 = core.replace("*","")
            strand2 = complement(strand1)
            strand1 = "5'-" + strand1[0:position-1] + base1 + strand1[position:] + "-3'"
            strand2 = "3'-" + strand2[0:position-1] + base2 + strand2[position:] + "-5'"

            df = pd.read_csv("/home/hs189/MMproj/TF_MD/Local/%s/%s_%s_steps.csv"%(tf,tf_name,folder))
            subdf = df.loc[(df.nt1_i.isin(nt1_list))]

            df_mean = subdf.groupby(['nt1_i'])[param].mean().reset_index().rename(columns={param:"mean"}) 
            df_mean = df_mean - 5.8
            df_std = subdf.groupby(['nt1_i'])[param].std().reset_index().rename(columns={param:"std"})
            df_param = pd.concat([df_mean,df_std['std']],axis=1)
            t = df_temp.merge(df_param,left_on='nt1_i',right_on='nt1_i',how='right')
            mean = t['mean'].values
            std = t['std'].values
            pos1_avg, pos2_avg, pos3_avg, pos4_avg = mean
            pos1_std, pos2_std, pos3_std, pos4_std = std

            output.append([tf_name,tf,seq_type,core,strand1,strand2,position,sim_type,mismatch,mm_id,binding,pos1_avg,pos1_std,pos2_avg,pos2_std,pos3_avg,pos3_std,pos4_avg,pos4_std])
            
    output_df = pd.DataFrame(output[1:],columns=output[0])
    print(output_df)
    output_df.to_csv("./Data/final_final_MD_mimic_%s.csv"%param,index=False,float_format='%.2f')
        
        
