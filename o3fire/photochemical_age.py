import numpy as np
import pandas as pd

R = 8.314472E-03

def bimolecular_rate_law(A, n, B, T=298):
    k = A * ( T / 298 )**n * np.exp( -B / (R*T) )
    return k


# X (compound that reacts with OH more quickly)
# Y (compound that reacts with OH less quickly)
def photochemical_age(df, measurements_df, X, Y, T, OH):
    """
    returns age in [s]
    """
    
    AX = eval(measurements_df.loc[measurements_df['Core_Name'] == X]['A'].values[0])
    nX = measurements_df.loc[measurements_df['Core_Name'] == X]['n'].values[0]
    BX = measurements_df.loc[measurements_df['Core_Name'] == X]['B'].values[0]
    EFX = measurements_df.loc[measurements_df['Core_Name'] == X]['EF'].values[0]
    
    AY = eval(measurements_df.loc[measurements_df['Core_Name'] == Y]['A'].values[0])
    nY = measurements_df.loc[measurements_df['Core_Name'] == Y]['n'].values[0]
    BY = measurements_df.loc[measurements_df['Core_Name'] == Y]['B'].values[0]
    EFY = measurements_df.loc[measurements_df['Core_Name'] == Y]['EF'].values[0]
    
    kX = bimolecular_rate_law(AX, nX, BX, T=T)
    kY = bimolecular_rate_law(AY, nY, BY, T=T)
    
    ER = EFX / EFY
    
    age = ( 1 / ( OH * (kX - kY ) ) * ( ( np.log(ER) ) -  ( np.log( df[X] / df[Y] ) ) ) )
    
    return age