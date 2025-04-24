import o3fire
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

############################################
############ Load and clean data ###########
############################################

#load incomplete cleaned data
firex_df_incomplete = pd.read_csv('../../data/raw/firex_incomplete.csv', parse_dates=['timestamp'], index_col='timestamp')
wecan_df_incomplete = pd.read_csv('../../data/raw/wecan_incomplete.csv', parse_dates=['timestamp'], index_col='timestamp')
arctas_df_incomplete = pd.read_csv('../../data/raw/arctas_incomplete.csv', parse_dates=['timestamp'], index_col='timestamp')
atom_df_incomplete = pd.read_csv('../../data/raw/atom_incomplete.csv', parse_dates=['timestamp'], index_col='timestamp')
dc3_60s_df_incomplete = pd.read_csv('../../data/raw/dc3_incomplete.csv', parse_dates=['timestamp'], index_col='timestamp')

#load complete cleaned data (missing value treatment)
firex_df_clean = pd.read_csv('../../data/gap_filled/firex_complete.csv', parse_dates=['timestamp'], index_col='timestamp')
wecan_df_clean = pd.read_csv('../../data/gap_filled/wecan_complete.csv', parse_dates=['timestamp'], index_col='timestamp')
arctas_df_clean = pd.read_csv('../../data/gap_filled/arctas_complete.csv', parse_dates=['timestamp'], index_col='timestamp')
atom_df_clean = pd.read_csv('../../data/gap_filled/atom_complete.csv', parse_dates=['timestamp'], index_col='timestamp')
dc3_60s_df_clean = pd.read_csv('../../data/gap_filled/dc3_complete.csv', parse_dates=['timestamp'], index_col='timestamp')

# add WECAN NOy
wecan_df_incomplete['NOy'] = wecan_df_incomplete[['NO', 'NO2', 'PAN', 'PPN', 'HNO3', ]].sum(axis=1)
wecan_df_clean['NOy'] = wecan_df_clean[['NO', 'NO2', 'PAN', 'PPN', 'HNO3', ]].sum(axis=1)

# altitude correction - to meters
wecan_df_clean['altitude'] = wecan_df_clean['altitude']*1000 
arctas_df_clean['altitude'] = arctas_df_clean['altitude']*1000
dc3_60s_df_clean['altitude'] = dc3_60s_df_clean['altitude']*1000

wecan_df_incomplete['altitude'] = wecan_df_incomplete['altitude']*1000
arctas_df_incomplete['altitude'] = arctas_df_incomplete['altitude']*1000
dc3_60s_df_incomplete['altitude'] = dc3_60s_df_incomplete['altitude']*1000

#load units data
firex_units = pd.read_csv("../../data/units/FIREX_TOGAmerge_units.csv", index_col=False)
wecan_units = pd.read_csv("../../data/units/WECAN_TOGAmerge_units.csv", index_col=False)
arctas_units = pd.read_csv("../../data/units/ARCTAS_TOGAmerge_units.csv", index_col=False)
atom_units = pd.read_csv("../../data/units/ATOM_TOGAmerge_units.csv", index_col=False)
dc3_units = pd.read_csv("../../data/units/DC3_TOGAmerge_units.csv", index_col=False)
dc3_60s_units = pd.read_csv("../../data/units/DC3_units.csv", index_col=False)

#load final measurement data
firex_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_firex_measurement_list.csv')
wecan_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_wecan_measurement_list.csv')
dc3_60s_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_dc3_60s_measurement_list.csv')
arctas_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_arctas_measurement_list.csv')
atom_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_atom_measurement_list.csv')

#remove stratospheric influence
firex_df_clean = firex_df_clean[((firex_df_clean['O3'] / firex_df_clean['H2O'])<1) & ((firex_df_clean['O3'] / firex_df_clean['H2O'])>0.003)]
atom_df_clean = atom_df_clean[((atom_df_clean['O3'] / atom_df_clean['H2O'])<1) & ((atom_df_clean['O3'] / atom_df_clean['H2O'])>0.003)]
wecan_df_clean = wecan_df_clean[((wecan_df_clean['O3'] / wecan_df_clean['H2O'])<1) & ((wecan_df_clean['O3'] / wecan_df_clean['H2O'])>0.003)]
arctas_df_clean = arctas_df_clean[((arctas_df_clean['O3'] / arctas_df_clean['H2O'])<1) & ((arctas_df_clean['O3'] / arctas_df_clean['H2O'])>0.003)]
dc3_60s_df_clean = dc3_60s_df_clean[((dc3_60s_df_clean['O3'] / dc3_60s_df_clean['H2O'])<1) & ((dc3_60s_df_clean['O3'] / dc3_60s_df_clean['H2O'])>0.003)]

#firex temperature correction
firex_df_clean['temperature'] = firex_df_clean['temp']+273.15

# ATom Seasonality -- keep only spring and summer months
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
        
atom_df_clean['season'] = atom_df_clean.reset_index().apply((lambda x: get_season(x['timestamp'].month, x['lat'])), axis=1).to_numpy()
# drop F&W
atom_df_clean = atom_df_clean[atom_df_clean['season'] == 'S&S']

### Combine campaigns

# Combine dataframes
firex_df_clean['campaign'] = 'FIREX'
wecan_df_clean['campaign'] = 'WECAN'
atom_df_clean['campaign'] = 'ATom'
arctas_df_clean['campaign'] = 'ARCTAS'
dc3_60s_df_clean['campaign'] = 'DC3'

measurements_df_dict = {'FIREX': firex_measurements_df,
                 'WECAN': wecan_measurements_df,
                 'ATom': atom_measurements_df,
                 'ARCTAS': arctas_measurements_df,
                 'DC3': dc3_60s_measurements_df,}
units_df_dict = {'FIREX': firex_units,
                 'WECAN': wecan_units,
                 'ATom': atom_units,
                 'ARCTAS': arctas_units,
                 'DC3': dc3_60s_units,}

full_df = pd.concat([firex_df_clean, wecan_df_clean, atom_df_clean, arctas_df_clean, dc3_60s_df_clean])

# correct units -- convert everything to ppt
columns = full_df.columns.to_list()
columns.remove('campaign')
columns.remove('lat')
columns.remove('long')
grouped_full_df = full_df.groupby('campaign')

print("Columns that couldn't be unit corrected:")
print()

for col in columns:
    #set everything to ppt
    target_unit = 'ppt'
    for i, (campaign, group) in enumerate(grouped_full_df):
        column_used = measurements_df_dict[campaign].loc[(measurements_df_dict[campaign]['Core_Name'] == col), 'Column_Name'].values
        if len(column_used)>0:
            try:
                unit = units_df_dict[campaign][column_used[0]].values[0] 
                if unit == 'missing':
                    unit = 'ppb'
                full_df.loc[(full_df['campaign']==campaign), col] = full_df.loc[(full_df['campaign']==campaign), col]*o3fire.utils.unit_conversion_ppx(unit, target_unit)
            except:
                print(campaign, col, column_used)
        else:
            pass


# derive extra columns
full_df['NOx'] = full_df['NO']+full_df['NO2']
full_df['Ox'] = full_df['O3']+full_df['NO2']

full_df = full_df.drop_duplicates()
full_df.sort_index(inplace=True)

#####################################
####### Baseline Subtraction ########
#####################################

#split up clean data -- 40th percentile of CO
clean_df = o3fire.baseline.split_clean(full_df)

# altitude dependent baseline
troposphere_alt_bins = np.array([0, 2, 4, 6, 8, 10, 12])*1000
troposphere_alt_ranges = troposphere_alt_bins[0:-1]

# subtract baseline
atom_trop_baseline, atom_enhancements_df = o3fire.baseline.subtract_baseline(full_df[full_df['campaign']=='ATom'], clean_df, method='altitude_dependent', baseline_percentile=0.25, remote=True)
arctas_trop_baseline, arctas_enhancements_df = o3fire.baseline.subtract_baseline((full_df[(full_df['campaign'] == 'ARCTAS')&(full_df['lat'] > 49)]), clean_df, method='altitude_dependent', baseline_percentile=0.25, arctic=True)
# arctas_trop_baseline, arctas_enhancements_df = o3fire.baseline.subtract_baseline((full_df[(full_df['lat'] > 49)]), clean_df, method='altitude_dependent', clean_method='percentile', baseline_percentile=0.25, arctic=True)
other_trop_baseline, other_enhancements_df = o3fire.baseline.subtract_baseline(full_df[(full_df['campaign']!='ATom')&(full_df['lat'] < 49)], clean_df, method='altitude_dependent', baseline_percentile=0.25,)

# concat to enhancements_df
campaign_enhancements_df = pd.concat([atom_enhancements_df, arctas_enhancements_df, other_enhancements_df]).reindex(full_df.index)

##########################################
############ Define Regimes ##############
##########################################
tmp_full = full_df.copy(deep=True)
tmp_enhancements = campaign_enhancements_df.copy(deep=True)

tmp_full['regime'] = o3fire.regimes.regime_definition(tmp_full, tmp_enhancements, method='relative_tracer', clean_method='percentile')            

#derive extra ratio columns

#raw data
tmp_full['HNO3/H2O2'] = tmp_full['HNO3']/tmp_full['H2O2']
tmp_full['O3/CO_raw'] = tmp_full['O3']/tmp_full['CO']
tmp_full['HNO3+PAN/H2O2'] = (tmp_full['HNO3']+tmp_full['PAN'])/tmp_full['H2O2']
tmp_full['HNO3/NOy'] = (tmp_full['HNO3'] / tmp_full['NOy'])
tmp_full['Ox/CO_raw'] = tmp_full['Ox']/tmp_full['CO']
tmp_full['NOx/NOy'] = tmp_full['NOx']/tmp_full['NOy']
tmp_full['PAN/NOy'] = tmp_full['PAN'] / tmp_full['NOy']

#enhancements data -- delO3/delCO
tmp_enhancements['O3/CO'] = tmp_enhancements['O3']/tmp_enhancements['CO']

df = tmp_full.merge(tmp_enhancements[['O3', 'Ox', 'CO', 'O3/CO', 'NOx', 'CH3CN', 'HCN', 'C2Cl4', 'CH2Cl2']], suffixes=['', '_delta'], left_index=True, right_index=True)

##########################################
######### Photochemical Age ##############
##########################################
df['age_phenol_benzene'] = (o3fire.photochemical_age.photochemical_age(df, firex_measurements_df, 'Phenol', 'Benzene', (df['temperature']), 1e6)/3600)
df['age_furan_benzene'] = (o3fire.photochemical_age.photochemical_age(df, firex_measurements_df, 'Furan', 'Benzene', (df['temperature']), 1e6)/3600)
df['age_toluene_benzene'] = (o3fire.photochemical_age.photochemical_age(df, firex_measurements_df, 'Toluene', 'Benzene', (df['temperature']), 1e6)/3600)
df['age'] = np.nan
df.loc[(df['regime'] == 'fire'),'age'] = df.loc[(df['regime'] == 'fire'), 'age_toluene_benzene']
df.loc[(df['regime'] == 'fire')&(df['campaign'] == 'WECAN'), 'age'] = df.loc[(df['regime'] == 'fire')&(df['campaign'] == 'WECAN'), 'age_furan_benzene']
df.loc[(df['regime'] == 'fire')&(df['campaign'] == 'FIREX'), 'age'] = df.loc[(df['regime'] == 'fire')&(df['campaign'] == 'FIREX'), 'age_phenol_benzene']

##########################################
########### OHR Calculation ##############
##########################################
df['pressure_atm'] = o3fire.ohr.hPa_to_Pa(o3fire.ohr.Pa_to_atm(df['pressure']))
df_ohr = o3fire.ohr.compute_ohr(df, firex_measurements_df,)

# compute theta prime
NOx_cols = ['NO', 'NO2']
ROC_cols = ['CO']
df_ohr['theta'] = (df_ohr[NOx_cols].sum(axis=1) / df_ohr[ROC_cols].sum(axis=1))
df = df.merge(df_ohr['theta'], left_index=True, right_index=True)

# save data
df.to_csv('../../data/processed/complete_data.csv')
df_ohr.to_csv('../../data/processed/ohr_data.csv')
