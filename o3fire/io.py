import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime
import csv
import matplotlib.pyplot as plt


def read_data_ict(campaign, tscol="Time_Start", cols_to_keep=[], extra_variables = [], file_pattern = "*.ict"):
    """
    Read in ict files from a campaign into a pandas DataFrame. Credit: Olivia Norman (MIT)
    """
    path = campaign + "/"
    
    df_list = []
    
    for i,file in enumerate(glob.glob(path+file_pattern)):

        lineStart = pd.read_csv(file, nrows = 1, header=None, encoding='latin1')[0][0] - 1

        single_flight_data = pd.read_csv(file, skiprows = lineStart, encoding='latin1')
            
        # get fillValue
        fillValue = pd.read_csv(file, skiprows = 11, nrows = 1, header = None, encoding='latin1')[0][0]

        # replace fillValue with NA
        single_flight_data[single_flight_data == fillValue] = np.nan
        
        # read in date of flight
        read_date = pd.read_csv(file, skiprows = 6, nrows = 1, header = None, encoding='latin1')
        date_info = [read_date[0][0],read_date[1][0],read_date[2][0]]
        flight_date = '-'.join([str(n) for n in date_info])
        
        # update timestamp to datetime -> timestamp is seconds from beginning of day
        time_col_name = tscol
        # compute timestamp first
        timestamp = pd.to_datetime(single_flight_data[time_col_name], unit='s', origin=pd.Timestamp(flight_date))

        # use concat instead of direct assignment to avoid fragmentation
        single_flight_data = pd.concat(
            [single_flight_data.copy(), pd.DataFrame({"timestamp": timestamp})],
            axis=1
        )

        # single_flight_data["timestamp"] = pd.to_datetime(single_flight_data[time_col_name], unit='s', origin=pd.Timestamp(flight_date))#datetime.utcfromtimestamp(single_flight_data[' UTC']).strftime('%Y-%m-%d %H:%M:%S')
        
        df_list.append(single_flight_data)
        
    # concat all the single flight dataframes into one
    campaign_data = pd.concat(df_list)
    campaign_data.set_index("timestamp", drop=True, inplace=True)
    campaign_data.sort_index(inplace=True)
    
    # remove special character
    campaign_data.columns = campaign_data.columns.str.replace(' ', '')
     
    return campaign_data


def get_units_ict(file, columns):
    unit_dict = {}

    # Open the file.
    with open(file, "r") as f:  # Get number of header rows
        header_row = int(f.readlines()[0].split(",")[0]) - 1

    with open(file, "r") as f:
        reader = f.readlines()
        ln_num = 0  # intitalize line counting var.
        for row in reader:
            line = row  # read line by line.

            for col in columns:
                if col in line:
                    try:
                        unit_dict[(line.split(',')[0].strip())] = [(line.split(',')[1].strip())]
                    except:
                        pass


            if ln_num > header_row - 2:
                break  # top once you reach data.
            ln_num = ln_num + 1
            
    return pd.DataFrame(unit_dict)