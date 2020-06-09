
import sys
import pandas as pd

df_param_free = pd.read_csv("../Survey/Survey_goldWC/stem/free/Canonical_all_bp_params.csv")
df_bend_free = pd.read_csv("../Survey/Survey_bend/free/Bend_Canonical_all_abg_params.csv")

df_free = pd.concat([df_param_free,df_bend_free],axis=0).reset_index(drop=True)
df_free = df_free[:-2]

df_free.at[7,'param_mean'] = df_free.at[7,'param_mean'] - 5.8
df_free.at[7,'param_max'] = df_free.at[7,'param_max'] - 5.8
df_free.at[7,'param_min'] = df_free.at[7,'param_min'] - 5.8
df_free.at[8,'param_mean'] = df_free.at[8,'param_mean'] - 5.8
df_free.at[8,'param_max'] = df_free.at[8,'param_max'] - 5.8
df_free.at[8,'param_min'] = df_free.at[8,'param_min'] - 5.8

df_param_bound = pd.read_csv("../Survey/Survey_goldWC/stem/bound/Canonical_all_bp_params.csv")
df_bend_bound = pd.read_csv("../Survey/Survey_bend/bound/Bend_Canonical_all_abg_params.csv")

df_bound = pd.concat([df_param_bound,df_bend_bound],axis=0).reset_index(drop=True)
df_bound = df_bound[:-2]

df_bound.at[7,'param_mean'] = df_bound.at[7,'param_mean'] - 5.8
df_bound.at[7,'param_max'] = df_bound.at[7,'param_max'] - 5.8
df_bound.at[7,'param_min'] = df_bound.at[7,'param_min'] - 5.8
df_bound.at[8,'param_mean'] = df_bound.at[8,'param_mean'] - 5.8
df_bound.at[8,'param_max'] = df_bound.at[8,'param_max'] - 5.8
df_bound.at[8,'param_min'] = df_bound.at[8,'param_min'] - 5.8


df_free.set_value(6,'param','c1c1_dist')
df_free.set_value(15,'param','beta')
df_free.set_value(16,'param','gamma')
df_free.set_value(17,'param','zeta')

df_bound.set_value(6,'param','c1c1_dist')
df_bound.set_value(15,'param','beta')
df_bound.set_value(16,'param','gamma')
df_bound.set_value(17,'param','zeta')

df_free.to_csv("./Data/final_envelope_free.csv",index=False,float_format='%.3f')
df_bound.to_csv("./Data/final_envelope_bound.csv",index=False,float_format='%.3f')


