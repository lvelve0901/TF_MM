
import sys
import pandas as pd
import numpy as np


for sigma in [0.5,1.0]:

    df1 = pd.read_csv("./Data/final_cutoff_mm_value_bp_step_sigma%s.csv"%sigma)
    df2 = pd.read_csv("./Data/cutoff_mm_value_abg_sigma%s_flank_-1.csv"%sigma,usecols=range(1,4)).rename(columns={'beta':'beta(-1.0)','gamma':'gamma(-1.0)','zeta':'zeta(-1.0)'})
    df3 = pd.read_csv("./Data/cutoff_mm_value_abg_sigma%s_flank_0.csv"%sigma,usecols=range(1,4)).rename(columns={'beta':'beta(0.0)','gamma':'gamma(0.0)','zeta':'zeta(0.0)'})
    df4 = pd.read_csv("./Data/cutoff_mm_value_abg_sigma%s_flank_1.csv"%sigma,usecols=range(1,4)).rename(columns={'beta':'beta(1.0)','gamma':'gamma(1.0)','zeta':'zeta(1.0)'})
    df5 = pd.read_csv("./Data/cutoff_mm_value_groove_sigma%s_flank_-2.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(-2.0)','majorgw':'majorgw(-2.0)'})
    df6 = pd.read_csv("./Data/cutoff_mm_value_groove_sigma%s_flank_-1.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(-1.0)','majorgw':'majorgw(-1.0)'})
    df7 = pd.read_csv("./Data/cutoff_mm_value_groove_sigma%s_flank_0.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(0.0)','majorgw':'majorgw(0.0)'})
    df8 = pd.read_csv("./Data/cutoff_mm_value_groove_sigma%s_flank_1.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(1.0)','majorgw':'majorgw(1.0)'})
    
    df = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8],axis=1)
    print(df)
    df.to_csv("./Data/final_cutoff_mm_value_all_sigma%.1f.csv"%sigma,index=False,float_format='%.0f')

    df1 = pd.read_csv("./Data/final_cutoff_mm_string_bp_step_sigma%s.csv"%sigma)
    df2 = pd.read_csv("./Data/cutoff_mm_string_abg_sigma%s_flank_-1.csv"%sigma,usecols=range(1,4)).rename(columns={'beta':'beta(-1.0)','gamma':'gamma(-1.0)','zeta':'zeta(-1.0)'})
    df3 = pd.read_csv("./Data/cutoff_mm_string_abg_sigma%s_flank_0.csv"%sigma,usecols=range(1,4)).rename(columns={'beta':'beta(0.0)','gamma':'gamma(0.0)','zeta':'zeta(0.0)'})
    df4 = pd.read_csv("./Data/cutoff_mm_string_abg_sigma%s_flank_1.csv"%sigma,usecols=range(1,4)).rename(columns={'beta':'beta(1.0)','gamma':'gamma(1.0)','zeta':'zeta(1.0)'})
    df5 = pd.read_csv("./Data/cutoff_mm_string_groove_sigma%s_flank_-2.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(-2.0)','majorgw':'majorgw(-2.0)'})
    df6 = pd.read_csv("./Data/cutoff_mm_string_groove_sigma%s_flank_-1.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(-1.0)','majorgw':'majorgw(-1.0)'})
    df7 = pd.read_csv("./Data/cutoff_mm_string_groove_sigma%s_flank_0.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(0.0)','majorgw':'majorgw(0.0)'})
    df8 = pd.read_csv("./Data/cutoff_mm_string_groove_sigma%s_flank_1.csv"%sigma,usecols=range(1,3)).rename(columns={'minorgw':'minorgw(1.0)','majorgw':'majorgw(1.0)'})
    
    df = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8],axis=1)
    print(df)
    df.to_csv("./Data/final_cutoff_mm_string_all_sigma%.1f.csv"%sigma,index=False)


