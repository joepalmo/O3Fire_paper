import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import statsmodels.api as sm
from .utils import NMB, reg_coef, rma_regression_slope

# get ICARTT file information
icartt_df = pd.read_csv('../../data/icartt.csv')

def get_measurement_list(df, rename_columns={}):
    measurements = {'Column_Name': [],
                'Core_Name': [],
                'Definition': [],
                'Molecular_Formula': [],
                'CAS_Number': [],
                'Specificity': [],
                'Category': []}
    
    tmp = df.copy(deep=True)
    
    for col in tmp.columns:
        if col in rename_columns.keys():
            column_name = rename_columns[col]
        else:
            column_name = col
        match, core_name = match_to_ict(column_name)
        if match:
            
            definition = icartt_df.loc[(icartt_df['CoreName'] == core_name), "Definition"].values[0]
            mf = icartt_df.loc[(icartt_df['CoreName'] == core_name), "Chemical Formula"].values[0]
            specificity = icartt_df.loc[(icartt_df['CoreName'] == core_name), "Specificity"].values[0]
            cas_no = icartt_df.loc[(icartt_df['CoreName'] == core_name), "CAS Number"].values[0]
            category = icartt_df.loc[(icartt_df['CoreName'] == core_name), "Category"].values[0]

            measurements['Column_Name'].append(col)
            measurements['Core_Name'].append(core_name)
            measurements['Definition'].append(definition)
            measurements['Molecular_Formula'].append(mf)
            measurements['CAS_Number'].append(cas_no)
            measurements['Specificity'].append(specificity)    
            measurements['Category'].append(category)
        
    return measurements

def match_to_ict(column_name):
    # i- --> iso
    mod_column_name = column_name.replace('i-', 'iso')
    for i,r in icartt_df.iterrows():
        core = r["CoreName"]
        definition = r["Definition"]
        if core == column_name:
            return True, core
        elif (core+"_" in column_name) and (column_name.startswith(core)):
            return True, core
        elif (definition.lower()+"_" in mod_column_name.lower()) and (mod_column_name.lower().startswith(definition.lower())):
            return True, core
        else:
            pass
    return False, None

def check_for_duplicates(full_measurements_df, core, col):
    duplicates = full_measurements_df.loc[(full_measurements_df['Core_Name']==core), 'Column_Name'].to_list()
    if col not in duplicates:
        duplicates.append(col)
    if len(duplicates)<1:
        return False, duplicates
    else:
        return True, duplicates
    

def unit_conversion(x, y):
    x = x.copy()
    y = y.copy()
    
    #cleaning
    x.loc[x == -888888] = np.nan
    y.loc[y == -888888] = np.nan
    
    x_clean = x[~np.isnan(x) & ~np.isnan(y)]
    y_clean = y[~np.isnan(x) & ~np.isnan(y)]
    
    ratio_val = 1
    
    if len(x_clean)>0:
        # find where 1:1 line should be
        ratio = np.nanmedian(y_clean/x_clean)
        # if y ppb and x ppt
        if ratio > 1e-6 and ratio < 0.01:
            ratio_val = 0.001
        # if x ppt and y ppb
        elif ratio > 20 and ratio < 2000:
            ratio_val = 1000
        # if x and y same
        elif ratio > 0.2 and ratio < 20:
            ratio_val = 1
    
    return ratio_val

def rank_columns(data_df, full_measurements_df, measurements_df, core, r_threshold=0.75, nmb_threshold=0.3, N_threshold=50, override=False):
    duplicate_dict = {'col': [],
                          'R2': [],
                          'N': [],
                          'NMB': []}
    
    
    col = measurements_df.loc[(measurements_df['Core_Name']==core), 'Column_Name'].values[0]
    # if column has duplicates
    duplicates_bool, duplicates = check_for_duplicates(full_measurements_df, core, col)
    if duplicates_bool:
        # loop through duplicates
        x = data_df[col].copy()
        rs = []
        ns = []
        nmbs = []
        rma_slopes = []
        duplicates.remove(col)
        for d in duplicates[:]:
            y = data_df[d].copy()
            # get unit conversion
            units = unit_conversion(x, y)
            r, n = reg_coef(x,(y/units))
            nmb = NMB(x, (y/units))
            rma_slope = rma_regression_slope(x, (y/units))
            # if suitable measurement agreement (R^2, NMB, N>?)
            # if (r > r_threshold) and (np.abs(nmb) < nmb_threshold) and n>N_threshold:
            if (r > r_threshold) and n>N_threshold and (rma_slope > 0.66) and (rma_slope < 1.5):
                rs.append(r)
                ns.append(n)
                nmbs.append(nmb)
                rma_slopes.append(rma_slope)
            else:
                duplicates.remove(d)
            duplicate_dict = {'col': duplicates,
                          'R2': rs,
                          'N': ns,
                          'NMB': nmbs,
                          'RMA_slope': rma_slopes}
            
        #rank the dict by R2
        duplicate_df = pd.DataFrame(duplicate_dict).sort_values('R2', axis=0, ascending=False,)
        return duplicate_df.to_dict('list')
    else:
        return duplicate_dict
    

def completeness(data_df, full_measurements_df, measurements_df, core, r_threshold=0.75, nmb_threshold=0.3, override=False):
    tmp = data_df.copy(deep=True)
    #add exceptions for PTR aggregate columns (monoterpenes, dimefurans)
    
    # get default col
    col = measurements_df.loc[(measurements_df['Core_Name']==core), 'Column_Name'].values[0]
    # get column rankings
    duplicate_dict = rank_columns(tmp, full_measurements_df, measurements_df, core, r_threshold=r_threshold, nmb_threshold=nmb_threshold, override=override)
    # replace NaNs in each column
    for c in duplicate_dict['col']:
        # get unit conversion
        units = unit_conversion(tmp[col], tmp[c])
        tmp[col].fillna((tmp[c] / units), inplace=True)
        
    return tmp[col]


def aggregate_campaign(data_df, full_measurements_df, measurements_df, extra_columns):
    
    incomplete_df = data_df.copy(deep=True)
    incomplete_df = incomplete_df[extra_columns.keys()]
    incomplete_df.rename(columns=extra_columns, inplace=True)
    complete_df = incomplete_df.copy(deep=True)
    
    for i,r in measurements_df.iterrows():
        try:
            incomplete_df[r['Core_Name']] = data_df[r['Column_Name']]
            complete_df[r['Core_Name']] = completeness(data_df, full_measurements_df, measurements_df, r['Core_Name'])
        except:
            print('Error with column: ', r['Core_Name'])
            pass
        
    return (incomplete_df, complete_df)