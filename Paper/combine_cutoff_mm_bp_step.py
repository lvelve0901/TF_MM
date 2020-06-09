
import pandas as pd
import numpy as np


for sigma in [0.5,1.0]:

    df1 = pd.read_csv("./Data/cutoff_mm_value_bp_sigma%.1f.csv"%sigma)
    df2 = pd.read_csv("./Data/cutoff_mm_value_step_sigma%.1f.csv"%sigma,usecols=range(1,7))
    
    df = pd.concat([df1,df2],axis=1)
    columns = ['mismatch','shear','stretch','stagger','buckle','propeller','opening','shift','slide','rise','tilt','roll','twist','C1C1_dist']
    df = df[columns]
    
    for param in columns[1:1+12]:
        for i in range(7,10):
            df.set_value(i,param,None)
    
    df.set_value(3,'shear',2.0)
    df.set_value(5,'shear',2.0)
    
    df = df.rename(columns={"C1C1_dist":"c1c1_dist"})
    print(df)
    df.to_csv("./Data/final_cutoff_mm_value_bp_step_sigma%.1f.csv"%sigma,index=False,float_format='%.0f')
    
    df = df.replace(1.0,'Pos')
    df = df.replace(-1.0,'Neg')
    df = df.replace(2.0,'Pos/Neg')
    
    print(df)
    df.to_csv("./Data/final_cutoff_mm_string_bp_step_sigma%.1f.csv"%sigma,index=False)


