# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:09:45 2023

@author: geolki
"""

import pandas as pd
from sys import platform
import pickle
import copy as cp
import calendar
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

# make dictionaries
prices_dict = {}; check_p = pd.DataFrame(index=[1], columns=[1])
for date in prices[['QUOTE_DATE']].drop_duplicates()['QUOTE_DATE']:
    prices_dict[date] = prices.loc[prices['QUOTE_DATE'] == date]
    temp = prices_dict[date].set_index('ITEM_ID')
    temp[date] = 1
    check_p = check_p.join(temp[[date]], how='outer')
check_p = check_p.drop(1, axis=0).drop(1, axis=1)
check_p.index = [str(x).split('.')[0] for x in check_p.index]
check_p.columns = [str(x).split('.')[0] for x in check_p.columns]
    
weights_dict = {}; check_w = pd.DataFrame(index=[1], columns=[1])
for date in weights[['INDEX_DATE']].drop_duplicates()['INDEX_DATE']:
    weights_dict[date] = weights.loc[weights['INDEX_DATE'] == date]
    temp = weights_dict[date]
    temp[date] = 1
    temp = temp[[date, 'ITEM_ID']].drop_duplicates().set_index(['ITEM_ID'])
    temp[date] = 1
    check_w = check_w.join(temp[[date]], how='outer')
check_w = check_w.drop(1, axis=0).drop(1, axis=1)
check_w.index = [str(x).split('.')[0] for x in check_w.index]
check_w.columns = [str(x).split('.')[0] for x in check_w.columns]

check_all = check_p.join(check_w, lsuffix='_p', rsuffix = '_w', how='outer').T
new_idx = []
for item in check_all.index:
    if '_' in item:
        new_idx.append(item)
check_all = check_all.loc[new_idx]
check_all['date'] = [x.split('_')[0] for x in check_all.index.tolist()]
check_all['data'] = [x.split('_')[1] for x in check_all.index.tolist()]
check_all = check_all.set_index(['date', 'data']).stack().dropna(how='all').unstack('data')
check_missing = check_all[check_all.isna().any(axis=1)]



# match items and weights
prices['ITEM_ID'] = [str(x).split('.')[0] for x in prices['ITEM_ID']]
prices['INDEX_DATE'] = [str(x).split('.')[0] for x in prices['QUOTE_DATE']]
weights['ITEM_ID'] = [str(x).split('.')[0] for x in weights['ITEM_ID']]
weights['INDEX_DATE'] = [str(x).split('.')[0] for x in weights['INDEX_DATE']]

weights['weight'] = 'weight'
prices['prices'] = 'prices'

basket = prices.merge(weights, on=['ITEM_ID', 'INDEX_DATE'], how='outer')

aaa = basket[['ITEM_ID', 'INDEX_DATE', 'weight', 'prices', 'ITEM_DESC']].drop_duplicates()
aaa = aaa[aaa.isna().any(axis=1)]


check_date = aaa[['INDEX_DATE', 'weight', 'prices']].drop_duplicates()
check_date = check_date[check_date.isna().any(axis=1)].fillna('')
check_date['check'] = check_date['prices'] + check_date['weight']
check_date['count'] = 1
check_date = check_date.set_index(['check', 'INDEX_DATE'])[['count']].unstack('check')
check_date = check_date[check_date.isna().any(axis=1)]


check_item = aaa[['ITEM_ID', 'ITEM_DESC', 'weight', 'prices']].drop_duplicates()
check_item = check_item[check_item.isna().any(axis=1)].fillna('')
check_item['check'] = check_item['prices'] + check_item['weight']
check_item['count'] = 1
check_item = check_item.set_index(['check', 'ITEM_ID', 'ITEM_DESC'])[['count']].unstack('check')
check_item = check_item[check_item.isna().any(axis=1)]

# import lookup to LCFS
lookup = pd.read_csv(data_filepath + 'lookups/basket_id_lookup.csv').fillna('')
lookup['Product'] = lookup['ITEM_DESC'] + ' ' + lookup['WEIGHT\SIZE']
lookup['ITEM_ID'] = [str(x).replace(' ', '') for x in lookup['ITEM_ID']]

# add lookup to weights
weights['ITEM_ID'] = [str(x).replace(' ', '') for x in weights['ITEM_ID']]
weights2 = weights.dropna(axis=1, how='all').merge(lookup, on='ITEM_ID', how='left')

wieghts3 = weights2[['order', 'ITEM_ID', 'ITEM_DESC_x', 'ITEM_DESC_y', 'lcfs']].drop_duplicates()
wieghts3 = wieghts3[wieghts3.isna().any(axis=1)]

# calculate basket emissions
temp = temp.join(lookup[['Product', 'lcfs', 'order']]).set_index(['lcfs', 'Product'])
multipliers_all.index.names = ['lcfs']
basket_price = temp[['order']]
for year in range(2018, 2022):
    temp2 = temp[[str(year)]].rename(columns={str(year):'Basket'}).join(multipliers_all[[year]].rename(columns={year:'multiplier'}))
    temp2[year] = temp2['multiplier'] * temp2['Basket']
    basket_price = basket_price.join(temp2[[year]])

basket_price = basket_price.sort_values('order').drop('order', axis=1).dropna(how='all')

basket_change = pd.DataFrame(index=basket_price.index)
comp_year = cp.copy(basket_price)
for year in years[1:]:
    comp_year.loc[comp_year[2018].isna() == True, 2018] = comp_year[year]
comp_year = comp_year[2018]

for year in years:
    basket_change[year] = basket_price[year] / comp_year * 100

##############
## Save All ##
##############

equ_hhd.to_csv(outputs_filepath + 'basket_2024/equivalised_household.csv')
cm_index.to_csv(outputs_filepath + 'basket_2024/carbon_multiplier_index.csv')
basket_change.to_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change.csv')
