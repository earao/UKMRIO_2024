# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:09:45 2023

@author: geolki
"""

import pandas as pd
from sys import platform
import pickle
import copy as cp
import os

# set working directory
# make different path depending on operating system

if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
data_filepath = wd + 'data/'
outputs_filepath = wd + 'outputs/'

# load data
hhd_ghg = pickle.load(open(outputs_filepath + 'results_2024/GHG_by_hhds.p', 'rb')) # emissions by household in survey
multipliers = pickle.load(open(outputs_filepath + 'results_2024/GHG_multipliers.p', 'rb'))

years = list(hhd_ghg.keys())

# fix multipliers
# fill in missing multipliers
multipliers_all = pd.DataFrame(index = multipliers[years[0]].index)
for year in list(multipliers.keys()):
    multipliers_all = multipliers_all.join(multipliers[year][['multipliers']].rename(columns={'multipliers':year}))

for year in multipliers_all.columns[1:]:
    multipliers_all.loc[multipliers_all[multipliers_all.columns[0]].isna() == True, multipliers_all.columns[0]] = multipliers_all[year]
for year in reversed(multipliers_all.columns[:-1]):
    multipliers_all.loc[multipliers_all[multipliers_all.columns[-1]].isna() == True, multipliers_all.columns[-1]] = multipliers_all[year]
for year in multipliers_all.columns[1:-1]:
    multipliers_all.loc[multipliers_all[year].isna() == True, year] = multipliers_all[[year-1, year+1]].mean(axis=1)


###########################
## Equivalised Household ##
###########################

equ_hhd = pd.DataFrame(columns=hhd_ghg[years[0]].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].columns)
for year in years:
    temp = cp.copy(hhd_ghg[year])
    if temp['OECD_mod'].sum() == 0:
        print(str(year) + ' missing OECD modifier')
    temp['pop'] = temp['OECD_mod'] * temp['weight']
    temp.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'] = temp.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.']\
        .apply(lambda x:x*temp['weight'])
    temp = temp.sum(0)
    pop = temp['pop']
    temp = (temp.loc['1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].astype(float) / pop.astype(float)) * 1.5
    temp['total_ghg'] = temp.sum()
    temp['year'] = year
    equ_hhd = equ_hhd.append(pd.DataFrame(temp).T)
equ_hhd = equ_hhd.set_index('year')   

#############################
## Carbon multiplier index ##
#############################

# CPI data: https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceinflation
# import CPI data
cpi = pd.read_excel(data_filepath + 'raw/CPI/202311_consumerpriceinflationdetailedreferencetables.xlsx', sheet_name='Table 9', index_col=2, header=5)\
    .drop(['Unnamed: 0', 'Unnamed: 1'], axis=1).dropna(how='all')
# map LCFS column names
lookup = pd.read_csv(data_filepath + 'lookups/cpi_to_lcfs.csv').set_index('CPI')
cpi = cpi.join(lookup).set_index('LCFS').mean(axis=0, level=0)
# change base year to 2010
cpi_comp = cpi[2010]
for year in cpi.columns:
    cpi[year] = cpi[year] / cpi_comp * 100

# calculate index
cm_index = pd.DataFrame(index = cpi.index)
for year in range(cpi.columns.min(), 2022):
    temp = cpi[[year]]
    temp['multiplier'] = multipliers_all[year]
    temp[year] = temp[year] * temp['multiplier']
    cm_index = cm_index.join(temp[[year]])

#####################
## Basket of goods ##
#####################

# Basket data (not used here, does not contain all items): https://www.ons.gov.uk/economy/inflationandpriceindices/articles/shoppingpricescomparisontool/2023-05-03
# Price and weights: https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceindicescpiandretailpricesindexrpiitemindicesandpricequotes

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
#basket['spend'] = basket['spend'] / basket['CPIH_COICOP_WEIGHT']

basket['year'] = [int(x[:4]) for x in basket['INDEX_DATE']]
basket = basket.set_index(['lcfs', 'year', 'INDEX_DATE'])    
   
# calculate basket emissions
temp = pd.DataFrame(multipliers_all.unstack())
temp.index.names = ['year', 'lcfs']; temp.columns = ['multiplier']

basket_ghg = basket.join(temp, how='left')
basket_ghg['ghg'] = basket_ghg['spend'] * basket_ghg['multiplier']

# get ghg emissions
basket_ghg = basket_ghg[['ghg']].unstack(level=['lcfs']).droplevel(axis=1, level=0).fillna(0)
order = []
for item in multipliers_all.index.tolist():
    if item in basket_ghg.columns:
        order.append(item)
        
basket_ghg = basket_ghg[order]
basket_ghg['TOTAL'] = basket_ghg.sum(1)
basket_ghg = basket_ghg.T

# calculate basket_change compared to 2015
basket_change = basket_ghg.apply(lambda x: x / basket_ghg[(2015, '201501')] * 100)

# at coicop 1 level
basket_ghg_ccp1 = cp.copy(basket_ghg)
basket_ghg_ccp1.index = [x.split('.')[0] for x in basket_ghg_ccp1.index]
basket_ghg_ccp1 = basket_ghg_ccp1.sum(axis=0, level=0)

# calculate basket_change compared to 2015
basket_change_ccp1 = basket_ghg_ccp1.apply(lambda x: x / basket_ghg_ccp1[(2015, '201501')] * 100).fillna(0)

##############
## Save All ##
##############

equ_hhd.to_csv(outputs_filepath + 'basket_2024/equivalised_household.csv')
cm_index.to_csv(outputs_filepath + 'basket_2024/carbon_multiplier_index.csv')
basket_change.to_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change.csv')
