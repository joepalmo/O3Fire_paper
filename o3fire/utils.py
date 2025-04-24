import numpy as np
from scipy.stats import pearsonr
import statsmodels.api as sm
import pandas as pd

def ppx_to_numden( data, ppx, P=1, T=298 ):
    """
    ppx with x = {m,b,t}
    P in atm
    T in K
    """
    ppx = ppx.lower().replace('v', '')
    ppx = ppx.replace('arb', '')
    ppx = ppx.replace('dry mole fraction mixing ratio', '')
    ppx = ppx.replace('(', '')
    ppx = ppx.replace(')', '')
    ppx = ppx.replace(' ', '')
    
    if ppx == 'missing':
        ppx = 'ppb'
    
    ppx_dict = {'ppm': 2.46e13,
                'ppb': 2.46e10,
                'ppt': 2.46e7,}
    
    numden = data * ppx_dict[ppx] * (P) * (298/T)
    
    return numden

def unit_conversion_ppx(unit, target_unit):
    unit = unit.lower().replace('v', '')
    target_unit = target_unit.lower().replace('v', '')
    # if unit ppb and target_unit ppt
    if ((unit == 'ppb') and (target_unit == 'ppt')) or ((unit == 'ppm') and (target_unit == 'ppb')):
        ratio_val = 1000
    elif ((unit == 'ppt') and (target_unit == 'ppb')) or ((unit == 'ppb') and (target_unit == 'ppm')):
        ratio_val = 0.001
    elif (unit == 'ppt') and (target_unit == 'ppm'):
        ratio_val = 1e-6
    elif (unit == 'ppm') and (target_unit == 'ppt'):
        ratio_val = 1e6
    else:
        ratio_val = 1
        
    return ratio_val

def get_season(month, lat):
    s_and_s_months = [4,5,6,7,8,9]
    f_and_w_months = [1,2,3,10,11,12]
        
    if month in s_and_s_months:
        if lat < 0:
            return 'F&W'
        else:
            return 'S&S'
    else:
        if lat < 0:
            return 'S&S'
        else:
            return 'F&W'
        
def percentile(n):
    def percentile_(x):
        return x.quantile(n)
    return percentile_

def reg_coef(x,y):
    x_clean = x[~np.isnan(x) & ~np.isnan(y)]
    y_clean = y[~np.isnan(x) & ~np.isnan(y)]

    try:
        r,p = pearsonr(x[~np.isnan(x) & ~np.isnan(y)], y[~np.isnan(x) & ~np.isnan(y)],)
        N = len(x[~np.isnan(x) & ~np.isnan(y)])
    except:
        r = np.nan
        N = np.nan
    return r, N

def NMB(M, O):
    if np.nanmean(M) >= np.nanmean(O):
        return (np.sum(M-O) / np.sum(O))
    else:
        return (np.sum(O-M) / np.sum(M))
    
def rma_regression_slope(x, y):
    try:
        slope_a = sm.OLS(y, x, missing='drop').fit().params[0]
        slope_b = sm.OLS(x, y, missing='drop').fit().params[0]
    except ValueError:  #raised if `y` is empty.
        return np.nan
    
    slope_b = 1 / slope_b
    # Check if correlated in same direction
    if np.sign(slope_a) != np.sign(slope_b):
        raise RuntimeError('Type I regressions of opposite sign.')
    
    # Compute Reduced Major Axis Slope
    slope = np.sign(slope_a) * np.sqrt(slope_a * slope_b)
    
    return slope

#### separate airborne campaign data into different flights
def separate_flights(df: pd.DataFrame, column_name = 'flight_no'):    
    flight_boundaries = df.loc[df.index.to_series().diff().dt.seconds > 180].index
    
    df['flight_no'] = 0
    for i,boundary in enumerate(flight_boundaries):
        if i==0:
            df.loc[:flight_boundaries[i], 'flight_no'] = 1
        elif i==len(flight_boundaries)-1:
            df.loc[flight_boundaries[i-1]:flight_boundaries[i], 'flight_no'] = i+1
            df.loc[flight_boundaries[i]:, 'flight_no'] = i+2
        else:
            df.loc[flight_boundaries[i-1]:flight_boundaries[i], 'flight_no'] = i+1
            
    return df