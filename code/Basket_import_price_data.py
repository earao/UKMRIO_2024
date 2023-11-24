# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:09:45 2023

@author: geolki
"""

import pandas as pd
from sys import platform
import os

# set working directory
# make different path depending on operating system

if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
data_filepath = wd + 'data/'

years = range(2010, 2022)

# import item prices
price_dir = data_filepath + 'raw/Basket_data/Prices/'
files = os.listdir(price_dir)

for year in years:
    for file in files:
        if str(year) in file:
            date = str(year) + file.split(str(year))[1].replace('_', '')[:2]
            if '.c' in date:
                pass
            else:
                price = pd.read_csv(price_dir + file); 
                price.columns = [x.upper().replace(' ', '') for x in price.columns]  
                
                # add date
                date = '20' + file.split('20')[1].replace('_', '')[:4]
                price['FILE_DATE'] = date
                price['FILE'] = file
               
                # calculate average
                price['PRICE_W'] = price['STRATUM_WEIGHT'] * price['PRICE']
                price = price.groupby(['QUOTE_DATE', 'ITEM_ID', 'FILE', 'FILE_DATE']).sum()[['STRATUM_WEIGHT', 'PRICE_W']]
                price['PRICE'] = price['PRICE_W'] / price['STRATUM_WEIGHT']
                
                # save file
                price.to_csv(data_filepath + 'processed/Basket_data/Prices/average_prices_' + date + '.csv')
                
                print('Done: ' + date)
