import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import statsmodels.api as sm

def percentile(n):
    def percentile_(x):
        return x.quantile(n)
    return percentile_

def split_clean(df, method='percentile', percentile=0.4):
    """
    split data into clean and polluted
    """
    #check that method exists
    possible_methods = ['percentile', 'tracers',]
    if method not in possible_methods:
        raise ValueError("method must be one of {} ; not '{}'".format(possible_methods, method))
    
    if method == 'percentile':
        clean_df = df.loc[(df['CO'] < df['CO'].quantile(percentile))]
    elif method == 'tracers':
        clean_df = df.loc[(df['CH3CN'] < df['CH3CN'].quantile(0.5))&(df['C2Cl4'] < df['C2Cl4'].quantile(0.5))]
        
    return clean_df

def subtract_baseline(df, clean_df, method='clean', baseline_percentile=None, remote=False, arctic=False):
    """
    make baseline subtraction and return enhancements_df
    """
    #check that method exists
    possible_methods = ['clean', 'altitude_dependent']
    if method not in possible_methods:
        raise ValueError("method must be one of {} ; not '{}'".format(possible_methods, method))
    
    if method == 'clean':
        # constant baseline
        trop_baseline = clean_df.agg(np.nanmedian)
        enhancements_df = df - trop_baseline
        
    elif method == 'altitude_dependent':
        # altitude dependent baseline
        troposphere_alt_bins = np.array([0, 2, 4, 6, 8, 10, 12])*1000
        troposphere_alt_ranges = troposphere_alt_bins[0:-1]

        ##################****************
        if remote:
            ft_tmp = clean_df[clean_df['campaign']=='ATom'].copy(deep=True)
        elif arctic:
            # clean_df = split_clean(df, method=clean_method, percentile=percentile)
            ft_tmp = clean_df[((clean_df['lat'])>49)].copy(deep=True)
        else:
            ft_tmp = clean_df[clean_df['campaign']!='ATom'].copy(deep=True)
            
        ft_tmp['alt_range'] =  pd.cut(ft_tmp['altitude'], troposphere_alt_bins)
        if baseline_percentile:
            trop_baseline = ft_tmp.groupby('alt_range').quantile(baseline_percentile)
        else:
            trop_baseline = ft_tmp.groupby('alt_range').agg(np.nanmedian)
        trop_baseline = trop_baseline.replace(np.nan,0)
        
        cols = list(trop_baseline.keys())
        cols.remove('altitude')
        
        enhancements_df = df.copy(deep=True)
        enhancements_df = enhancements_df[~enhancements_df['altitude'].isna()]
        enhancements_df['alt_range'] =  pd.cut(enhancements_df['altitude'], troposphere_alt_bins)
        
        # Subtract baseline
        for col in cols:
            for interval in ft_tmp['alt_range'].unique():
                if interval is np.nan:
                    pass
                else:
                    save_interval = interval
                    enhancements_df.loc[enhancements_df["alt_range"] == interval, col] = enhancements_df.loc[enhancements_df["alt_range"] == interval, col] - trop_baseline.loc[interval, col]
            enhancements_df.loc[enhancements_df["altitude"] > 12000, col] = enhancements_df.loc[enhancements_df["altitude"] > 12000, col] - trop_baseline.loc[save_interval, col]
        enhancements_df.drop('alt_range', axis=1, inplace=True)
        
    return trop_baseline, enhancements_df