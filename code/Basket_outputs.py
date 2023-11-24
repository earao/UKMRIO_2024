# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:09:45 2023

@author: geolki
"""

import pandas as pd
from sys import platform
import pickle
import copy as cp

# set working directory
# make different path depending on operating system

if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
data_filepath = wd + 'data/'
results_filepath = wd + 'outputs/results_2024/'
lcf_filepath = wd + 'data/raw/LCFS/'
model_inputs = wd + 'data/model_inputs/'

# load data
hhd_ghg = pickle.load(open(results_filepath + 'GHG_by_hhds.p', 'rb')) # emissions by household in survey
multipliers = pickle.load(open(results_filepath + 'GHG_multipliers.p', 'rb'))

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

# Basket data: https://www.ons.gov.uk/economy/inflationandpriceindices/articles/shoppingpricescomparisontool/2023-05-03

# import basket price
temp = pd.read_excel(data_filepath + 'raw/Basket_data/2023_basket.xlsx', sheet_name='averageprice', index_col=0).dropna(axis=0, how='all')
# make yearly average
temp.columns = [str(x)[:4] for x in temp.columns]
temp = temp.mean(axis=1, level=0, skipna=True)

# import lookup to LCFS
lookup = pd.read_csv(data_filepath + 'lookups/basket_id_lookup.csv').set_index('ITEM_ID').fillna('')
lookup['Product'] = lookup['ITEM_DESC'] + ' ' + lookup['WEIGHT\SIZE']

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
comp_year = basket_price[2018]
for year in range(2018, 2022):
    basket_change[year] = basket_price[year] / comp_year * 100
    