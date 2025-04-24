import o3fire
import numpy as np
import pandas as pd

# Load raw data -- TOGA merge files available online for each campaign
# get urls

# replace with path to raw .ict data
path_to_firex = '/home/jpalmo/fs09/Projects/O3Fire/data/raw/FIREX_instruments_raw/TOGA_merge/'
path_to_wecan = '/home/jpalmo/fs09/Projects/O3Fire/data/raw/WECAN_TOGA_merge/'
path_to_arctas = '/home/jpalmo/fs09/Projects/O3Fire/data/raw/ARCTAS_TOGA_merge/'
path_to_dc3 = '/home/jpalmo/fs09/Projects/O3Fire/data/raw/DC3_TOGA_merge/'
path_to_atom = '/home/jpalmo/fs09/Projects/O3Fire/data/raw/ATom_TOGA_merge/'
path_to_dc3_60s = '/home/jpalmo/fs09/Projects/O3Fire/data/raw/DC3/'

# read in data
firex_df = o3fire.io.read_data_ict(path_to_firex)
arctas_df = o3fire.io.read_data_ict(path_to_arctas, tscol='UTC', file_pattern="*R14.ict")
wecan_df = o3fire.io.read_data_ict(path_to_wecan, tscol=' UTC_mid', file_pattern="*R4.ict")
dc3_df = o3fire.io.read_data_ict(path_to_dc3, tscol=' UTC', file_pattern="*R6.ict")
dc3_60s_df = o3fire.io.read_data_ict(path_to_dc3_60s, tscol=' UTC')
atom_df = o3fire.io.read_data_ict(path_to_atom, tscol='UTC_Start')

firex_units_df = o3fire.io.get_units_ict(path_to_firex+'firexaq-mrgTOGA-dc8_merge_20190722_R1.ict', firex_df.columns)
firex_units_df.to_csv("../../data/units/FIREX_TOGAmerge_units.csv")
wecan_units_df = o3fire.io.get_units_ict(path_to_wecan+'wecan-mrgTOGA-c130_merge_20180724_R4.ict', wecan_df.columns)
wecan_units_df.to_csv("../../data/units/WECAN_TOGAmerge_units.csv")
arctas_units_df = o3fire.io.get_units_ict(path_to_arctas+'ARCTAS-mrgTOGA-dc8_merge_20080401_R14.ict', arctas_df.columns)
arctas_units_df.to_csv("../../data/units/ARCTAS_TOGAmerge_units.csv")
atom_units_df = o3fire.io.get_units_ict(path_to_atom+'MER-TOGA_DC8_20160729_R22.ict', atom_df.columns)
atom_units_df.to_csv("../../data/units/ATOM_TOGAmerge_units.csv")
dc3_60s_units_df = o3fire.io.get_units_ict(path_to_dc3_60s+'dc3-mrg60-dc8_merge_20120518_R10.ict', dc3_60s_df.columns)
dc3_60s_units_df.to_csv("../../data/units/DC3_units.csv")
dc3_units_df = o3fire.io.get_units_ict(path_to_dc3+'dc3-mrgTOGA-gV_merge_20120518_R6.ict', dc3_df.columns)
dc3_units_df.to_csv("../../data/units/DC3_TOGAmerge_units.csv")

#cleaning
firex_df[firex_df == -888888] = np.nan
wecan_df[wecan_df == -8888888] = np.nan
wecan_df[wecan_df == -9.999999e6] = np.nan
arctas_df[arctas_df == -999999999.0] = np.nan
arctas_df[arctas_df == -888888888] = np.nan
dc3_df[dc3_df == -888888] = np.nan
dc3_60s_df[dc3_60s_df == -888888] = np.nan
atom_df[atom_df == -8888] = np.nan
atom_df[atom_df == -7777] = np.nan
atom_df[atom_df == -99999] = np.nan

# use ICARTT file information to rename columns from each campaign to common names -- the ICARTT CoreName

# First rename columns that don't adhere to ICARTT standards
# Rename problematic OHR columns
dc3_rename = {'OH_Reactivity': 'OHR', 'OH_Reactivity_TubeTemp':'OHRTubeTemp', 'OH_Reactivity_TubePress':'OHRTubePress'}
# Rename Wennberg columns we need that are named irregularly
firex_rename = {'H2O_DLH_DISKIN':'H2OMRV_DISKIN',
                        'HNO3-1Hz_CIT_WENNBERG': 'HNO3_1Hz_CIT_WENNBERG',
                        'HCN-1Hz_CIT_WENNBERG': 'HCN_1Hz_CIT_WENNBERG',
                         'BUTENE-HN-1Hz_CIT_WENNBERG': 'C4H9NO4_1Hz_CIT_WENNBERG',
                         'ETHENE-HN-1Hz_CIT_WENNBERG': 'x2HydEthONO2_1Hz_CIT_WENNBERG',
                         'PHENOL-1Hz_CIT_WENNBERG': 'Phenol_1Hz_CIT_WENNBERG',
                         'H2O2-1Hz_CIT_WENNBERG': 'H2O2_1Hz_CIT_WENNBERG',
                         'ISOPN-1Hz_CIT_WENNBERG': 'ISOPN_1Hz_CIT_WENNBERG',
                         'PROPENE-HN-1Hz_CIT_WENNBERG': 'C3H7NO4_HN-1Hz_CIT_WENNBERG',}
arctas_rename = {'GT_PAN':'PAN_GT', 'GT_PPN':'PPN_GT', 'GT_APAN':'APAN_GT', 'GT_PIBN':'PIBN_GT', 'Carbon_Monoxide_mixing_ratio':'CO'}
atom_rename = {'OHReactivity_OHR': 'OHR'}

# Get measurement list for each campaign
firex_measurements_list = o3fire.completeness.get_measurement_list(firex_df, rename_columns=firex_rename)
firex_measurements_df = pd.DataFrame(firex_measurements_list)
firex_measurements_df[firex_measurements_df['Core_Name']=='iButene']
wecan_measurements_list = o3fire.completeness.get_measurement_list(wecan_df, rename_columns={})
wecan_measurements_df = pd.DataFrame(wecan_measurements_list)
arctas_measurements_list = o3fire.completeness.get_measurement_list(arctas_df, rename_columns=arctas_rename)
arctas_measurements_df = pd.DataFrame(arctas_measurements_list)
dc3_measurements_list = o3fire.completeness.get_measurement_list(dc3_df, rename_columns=dc3_rename)
dc3_measurements_df = pd.DataFrame(dc3_measurements_list)
dc3_60s_measurements_list = o3fire.completeness.get_measurement_list(dc3_60s_df, rename_columns=dc3_rename)
dc3_60s_measurements_df = pd.DataFrame(dc3_60s_measurements_list)
atom_measurements_list = o3fire.completeness.get_measurement_list(atom_df, rename_columns=atom_rename)
atom_measurements_df = pd.DataFrame(atom_measurements_list)
# remove GMI model output data
atom_measurements_df = atom_measurements_df[~atom_measurements_df['Column_Name'].str.contains('GMI')]

# Add water vapor columns
measurement_cols = wecan_measurements_df.columns.to_list()
wecan_water_cols = [['H2O_PICARRO', 'H2O_QCL'], ['H2O', 'H2O'], ['water vapor', 'water vapor'], ['H2O', 'H2O'], [np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]]
full_wecan_measurements_df = pd.concat([wecan_measurements_df, pd.DataFrame(dict(zip(measurement_cols, wecan_water_cols)))], ignore_index=True)
arctas_water_cols = [['H2O(v)'], ['H2O'], ['water vapor'], ['H2O'], [np.nan],[np.nan],[np.nan]]
full_arctas_measurements_df = pd.concat([arctas_measurements_df, pd.DataFrame(dict(zip(measurement_cols, arctas_water_cols)))], ignore_index=True)
dc3_60s_water_cols = [['H2O_vapor_DLH'], ['H2O'], ['water vapor'], ['H2O'], [np.nan],[np.nan],[np.nan]]
full_dc3_60s_measurements_df = pd.concat([dc3_60s_measurements_df, pd.DataFrame(dict(zip(measurement_cols, dc3_60s_water_cols)))], ignore_index=True)
dc3_water_cols = [['H2O_MixingRatio_VCSEL'], ['H2O'], ['water vapor'], ['H2O'], [np.nan],[np.nan],[np.nan]]
full_dc3_measurements_df = pd.concat([dc3_measurements_df, pd.DataFrame(dict(zip(measurement_cols, dc3_water_cols)))], ignore_index=True)
atom_water_cols = [['H2O_DLH'], ['H2O',], ['water vapor',], ['H2O',], [np.nan,], [np.nan,], [np.nan,]]
full_atom_measurements_df = pd.concat([atom_measurements_df, pd.DataFrame(dict(zip(measurement_cols, atom_water_cols)))], ignore_index=True)
full_firex_measurements_df = firex_measurements_df

# final measurement lists -- these contain the first choice columns for each campaign, chosen based on coverage and measurement quality
# pared down measurement lists
firex_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_firex_measurement_list.csv')
wecan_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_wecan_measurement_list.csv')
dc3_60s_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_dc3_60s_measurement_list.csv')
arctas_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_arctas_measurement_list.csv')
atom_measurements_df = pd.read_csv('../../data/final_measurement_lists/final_atom_measurement_list.csv')

# print(full_firex_measurements_df)

#### FIREX-AQ ####
firex_extra_columns = {'Latitude_YANG': 'lat', 'Longitude_YANG': 'long', 'MSL_GPS_Altitude_YANG':'altitude',
                 'jNO2_NO_O3P_CAFS_HALL': 'jNO2', 'jO3_O2_O1D_CAFS_HALL': 'jO3', 'Static_Pressure_YANG': 'pressure', 'Static_Air_Temp_YANG':'temperature',
                'Wind_Speed_YANG': 'wind_speed', 'Wind_Direction_YANG': 'wind_direction', 'OA_PM1_AMS_JIMENEZ': 'OA_AMS'}

firex_incomplete_df, firex_complete_df = o3fire.completeness.aggregate_campaign(firex_df, full_firex_measurements_df, firex_measurements_df, firex_extra_columns)

#### WE-CAN ####
wecan_extra_columns = {'LATITUDE': 'lat', 'LONGITUDE': 'long', 'GPS_ALT':'altitude',
                 'j[NO2->NO+O(3P)]': 'jNO2', 'j[O3->O2+O(1D)]': 'jO3', 'PRESSURE': 'pressure', 'TEMPERATURE':'temperature',
                'WNS': 'wind_speed', 'WND': 'wind_direction', 'OA_AMS': 'OA_AMS'}

wecan_incomplete_df, wecan_complete_df = o3fire.completeness.aggregate_campaign(wecan_df, full_wecan_measurements_df, wecan_measurements_df, wecan_extra_columns)

#### ARCTAS ####
arctas_extra_columns = {'LATITUDE': 'lat', 'LONGITUDE': 'long', 'GPS_Altitude':'altitude',
                 'J[NO2->NO+O(3P)]': 'jNO2', 'J[O3->O2+O(1D)]': 'jO3', 'PRESSURE': 'pressure', 'TEMPERATURE':'temperature',
                'WNS': 'wind_speed', 'WND': 'wind_direction', 'Organics<213mz': 'OA_AMS'}

arctas_incomplete_df, arctas_complete_df = o3fire.completeness.aggregate_campaign(arctas_df, full_arctas_measurements_df, arctas_measurements_df, arctas_extra_columns)

#### ATOM ####
atom_extra_columns = {'G_LAT': 'lat', 'G_LONG': 'long', 'G_ALT':'altitude',
                 'jNO2_NO_O3P_CAFS': 'jNO2', 'jO3_O2_O1D_CAFS': 'jO3', 'P': 'pressure', 'T':'temperature',
                'Wind_Speed': 'wind_speed', 'Wind_Direction': 'wind_direction', 'OA_PM1_AMSSD': 'OA_AMS'}

atom_incomplete_df, atom_complete_df = o3fire.completeness.aggregate_campaign(atom_df, full_atom_measurements_df, atom_measurements_df, atom_extra_columns)

#### DC3 ####
dc3_60s_extra_columns = {'LATITUDE': 'lat', 'LONGITUDE': 'long', 'GPS_ALT':'altitude',
                 'J[NO2->NO+O(3P)]': 'jNO2', 'J[O3->O2+O(1D)]': 'jO3', 'PRESSURE': 'pressure', 'TEMPERATURE':'temperature',
                'WNS': 'wind_speed', 'WND': 'wind_direction', 'Org_lt_1um_AMS': 'OA_AMS'}

dc3_60s_incomplete_df, dc3_60s_complete_df = o3fire.completeness.aggregate_campaign(dc3_60s_df, full_dc3_60s_measurements_df, dc3_60s_measurements_df, dc3_60s_extra_columns)
# resample to 1min45s to match other campaigns (TOGA merge files)
dc3_60s_complete_df = dc3_60s_complete_df.resample('105s').apply(np.nanmean)
dc3_60s_incomplete_df = dc3_60s_incomplete_df.resample('105s').apply(np.nanmean)
# remove rows which are all NaN
dc3_60s_complete_df.dropna(how='all', inplace=True)
dc3_60s_incomplete_df.dropna(how='all', inplace=True)

## Save complete dataframes
firex_complete_df.to_csv('../../data/gap_filled/firex_complete.csv')
wecan_complete_df.to_csv('../../data/gap_filled/wecan_complete.csv')
arctas_complete_df.to_csv('../../data/gap_filled/arctas_complete.csv')
atom_complete_df.to_csv('../../data/gap_filled/atom_complete.csv')
dc3_60s_complete_df.to_csv('../../data/gap_filled/dc3_complete.csv')

## save incomplete dataframes
firex_incomplete_df.to_csv('../../data/raw/firex_incomplete.csv')
wecan_incomplete_df.to_csv('../../data/raw/wecan_incomplete.csv')
arctas_incomplete_df.to_csv('../../data/raw/arctas_incomplete.csv')
atom_incomplete_df.to_csv('../../data/raw/atom_incomplete.csv')
dc3_60s_incomplete_df.to_csv('../../data/raw/dc3_incomplete.csv')