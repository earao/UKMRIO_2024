# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:09:45 2023

@author: geolki
"""

# Basket data (not used here, does not contain all items): https://www.ons.gov.uk/economy/inflationandpriceindices/articles/shoppingpricescomparisontool/2023-05-03
# Price and weights: https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceindicescpiandretailpricesindexrpiitemindicesandpricequotes


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

years = list(range(2010, 2022))

# import item prices
prices_dir = data_filepath + 'processed/Basket_data/Prices/'
files = os.listdir(prices_dir)
prices = pd.DataFrame()
for file in files:
    prices = prices.append(pd.read_csv(prices_dir + file).fillna(0));

# import item weights
weight_dir = data_filepath + 'raw/Basket_data/Weights/'
files = os.listdir(weight_dir)
weights = pd.DataFrame()
for file in files:
    for year in years:
        if str(year) in file:
            temp = pd.read_csv(weight_dir + file)
            temp.columns = [x.upper().replace(' ', '') for x in temp.columns]
            weights = weights.append(temp);
weights = weights.drop(['IMPUTATION_FLAG', 'IMPUTATION_FLAGS'], axis=1).dropna(how='all', axis=0).dropna(how='all', axis=1)

# match items and weights
prices['ITEM_ID'] = [str(x).split('.')[0] for x in prices['ITEM_ID']]
prices['INDEX_DATE'] = [str(x).split('.')[0] for x in prices['QUOTE_DATE']]
# fix missing years for weights
missing_years = ['201703', '201704', '201705', '201706', '201707', '201708', '201709', '201710']
prices.loc[prices['INDEX_DATE'].isin(missing_years) == True, 'INDEX_DATE'] = '201711'

weights['ITEM_ID'] = [str(x).split('.')[0] for x in weights['ITEM_ID']]
weights['INDEX_DATE'] = [str(x).split('.')[0] for x in weights['INDEX_DATE']]

# merge prices and weights
keep = ['INDEX_DATE', 'ITEM_ID', 'ITEM_DESC', 'QUOTE_DATE', 'PRICE', 'COICOP_WEIGHT', 'ITEM_WEIGHT', 'CPIH_COICOP_WEIGHT', 'ITEM_INDEX']
basket = prices.merge(weights, on=['ITEM_ID', 'INDEX_DATE'], how='outer').fillna(0)[keep]
basket['spend'] = basket['PRICE'] * basket['CPIH_COICOP_WEIGHT']
new_desc = []
for item in basket['ITEM_DESC']:
    new = str(item)
    while new[-1] == ' ':
        new = new[:-1]
    new_desc.append(new)
basket['ITEM_DESC'] = new_desc

# iport coicop lookup
lookup = pd.read_csv(data_filepath + 'lookups/basket_to_ccp4.csv')
lookup['ITEM_ID'] = lookup['ITEM_ID'].astype(str)

# calculate basket price by coicop
basket = basket.merge(lookup[['ITEM_ID', 'lcfs']], on='ITEM_ID', how='left')\
    .groupby(['INDEX_DATE', 'lcfs']).sum()[['spend', 'CPIH_COICOP_WEIGHT']].reset_index()
basket['spend'] = basket['spend'] / basket['CPIH_COICOP_WEIGHT'] # check method of how CPIH weights are used for CPI

basket['year'] = [int(x[:4]) for x in basket['INDEX_DATE']]
basket = basket.set_index(['lcfs', 'year', 'INDEX_DATE'])  
basket.to_csv(data_filepath + 'processed/Basket_data/Basket_spends.csv')
  