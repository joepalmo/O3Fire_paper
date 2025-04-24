import numpy as np
import pandas as pd
from .utils import ppx_to_numden

R = 8.314472E-03

def Pa_to_atm(data):
    atm = data / 101325
    return atm

def hPa_to_Pa(data):
    Pa = data * 100
    return Pa

def C_to_K(data):
    return (data+273.15)

def bimolecular_rate_law(A, n, B, T=298):
    k = A * ( T / 298 )**n * np.exp( -B / (R*T) )
    return k

def compute_ohr(data_df, measurements_df):
    ohr_df = pd.DataFrame(index=data_df.index, columns=measurements_df['Core_Name'].to_numpy())
    
    # initialize OHR column
    ohr_df['OHR'] = 0

    isomers = {}

    tmp = data_df.copy(deep=True)
    for i,m in measurements_df.iterrows():

        if m['Molecular_Formula'] in isomers.keys():
            continue
        
        # save the name of the column for this measurement
        col = m['Core_Name']
        
        #check for isomer
        if (m['A'] == 'isomer'):
            isomers[m['Molecular_Formula']] = m['Core_Name']
            continue    

        # tmp[col] = data_df[m['Core_Name']]

        # create temporary dataframe to house VOC and rxn rate info for each OH reaction
        tmp_df = pd.DataFrame(index=data_df.index)

        # create VOC column to be filled
        tmp_df['VOC'] = 0.0

        # set negative values to 0
        tmp.loc[(tmp[col] < 0), col] = 0
        
        # convert to number density
        tmp_col = ppx_to_numden(tmp[col], 'ppt', P=data_df['pressure'], T=data_df['temperature'])

        # fill nans with 0 so that we can do addition
        tmp_col = tmp_col.fillna(0)

        # add VOC contribution
        tmp_df['VOC'] = tmp_df['VOC']+tmp_col
        
        # get the reaction rate depending on temperature
        try:
            tmp_df['k'] = bimolecular_rate_law(eval(m['A']), m['n'], m['B'], T=tmp['temperature'])
        except:
            tmp_df['k'] = bimolecular_rate_law(m['A'], m['n'], m['B'], T=tmp['temperature'])

        # compute the OHR contribution from this reaction
        tmp_df['OHR'] = tmp_df['VOC']*tmp_df['k']

        # save the OHR contribution for this reaction under the core name identifier
        ohr_df[m['Core_Name']] = tmp_df['OHR']

        # add the OHR contribution to the total
        ohr_df['OHR'] = ohr_df['OHR'] + tmp_df['OHR']
    
    return ohr_df