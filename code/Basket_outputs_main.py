# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:09:45 2023

@author: geolki

Need to run before:
    1. Basket_emissions_main.py
    2. LCFS_data_main_2024.py
    3. Basket_import_price_data.py
    4. Basket_make_spend.py
"""

import pandas as pd
from sys import platform
import pickle
import copy as cp
import numpy as np
import Basket_outputs_functions as bof

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

years = list(multipliers.keys())

# fix multipliers
# fill in missing multipliers
multipliers_all = pd.DataFrame(index = multipliers[years[0]].index)
for year in years:
    multipliers_all = multipliers_all.join(multipliers[year][['multipliers']].rename(columns={'multipliers':year}))

for year in multipliers_all.columns[1:]:
    multipliers_all.loc[multipliers_all[multipliers_all.columns[0]].isna() == True, multipliers_all.columns[0]] = multipliers_all[year]
for year in reversed(multipliers_all.columns[:-1]):
    multipliers_all.loc[multipliers_all[multipliers_all.columns[-1]].isna() == True, multipliers_all.columns[-1]] = multipliers_all[year]
for year in multipliers_all.columns[1:-1]:
    multipliers_all.loc[multipliers_all[year].isna() == True, year] = multipliers_all[[year-1, year+1]].mean(axis=1)


# make index dictionary from multipliers_all
idx_dict = dict(zip([x.split(' ')[0] for x in multipliers_all.index], multipliers_all.index.tolist()))

## import CPI, which is needed for SDA and CM analysis
# CPI data: https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceinflation
# import CPI data
cpi = bof.import_cpi(data_filepath, idx_dict)

###########################
## Equivalised Household ##
###########################

years = list(hhd_ghg.keys())

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

years = list(hhd_ghg.keys())

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

# import basket data
basket = pd.read_csv(data_filepath + 'processed/Basket_data/Basket_spends.csv')
basket['lcfs'] = [x.split(' ')[0] for x in basket['lcfs']]
basket = basket.set_index(['lcfs', 'year', 'INDEX_DATE'])   

basket_test = basket[['spend']].unstack(level='lcfs').T.mean(axis=1, level=0)

# calculate basket emissions
temp = pd.DataFrame(multipliers_all.unstack())
temp['lcfs'] = [x[1].split(' ')[0] for x in temp.index.tolist()]
temp.index.names = ['year', 'full_label']
temp = temp.set_index('lcfs', append=True); temp.columns = ['multiplier']

basket_ghg = basket.join(temp, how='left')
basket_ghg['ghg'] = basket_ghg['spend'] * basket_ghg['multiplier'] 

# get ghg emissions
basket_ghg = basket_ghg[['ghg']].unstack(level=['lcfs', 'full_label']).droplevel(axis=1, level=[0, 1]).fillna(0)
order = []
for item in multipliers_all.index.tolist():
    if item in basket_ghg.columns:
        order.append(item)
        
basket_ghg = basket_ghg[order]
basket_ghg['TOTAL'] = basket_ghg.sum(1)
basket_ghg = basket_ghg.T

# calculate basket_change compared to 2015
basket_change = basket_ghg.apply(lambda x: x / basket_ghg[(2015, 201501)] * 100)

#########
## SDA ##
#########

# Population data: https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland

# use only household emission
# Structure: pop (1x1) * emission intensities, deflated (1x105) * prop spend from yhh (105x1) * total spend from yhh, deflated (1x1)

years = list(hhd_ghg.keys())

# get total spend to make other data (from yhh)
spend = pickle.load(open(outputs_filepath + 'results_2024/SPEND_yhh.p', 'rb'))
# get toal spend per year by product
total_spend = bof.make_total(spend, years, idx_dict)
# get toal emissions per year per product
temp = {year:hhd_ghg[year].apply(lambda x: pd.to_numeric(x, errors='coerce')*hhd_ghg[year]['weight'])[list(idx_dict.values())] 
        for year in years}
total_co2 = bof.make_total(hhd_ghg, years, idx_dict)

# proportional spend          
prop_spend = total_spend.apply(lambda x: x/x.sum())

# deflated total spend
defl_spend = pd.DataFrame(index = total_spend.index)
for year in years:
    defl_spend[year] = total_spend[year] * (cpi[year] / 100)
defl_tot_spend = defl_spend.sum(axis=0)

# deflated emission 
defl_co2 = pd.DataFrame(total_co2 / (defl_spend + 0.001)).astype(float)
defl_co2 = defl_co2.loc[multipliers_all.index]

# population
population = pd.read_csv(data_filepath + 'raw/Population/series-281123.csv', index_col=0, header=7)
population.columns = ['population']

# make SDA variables
sda_vars = {}
footprint = {}
for year in years:
    temp = {}
    # pop (1x1)
    temp['population'] = np.array(population.loc[year, 'population'])
    # emission intensities, deflated (1x105)
    temp['ghg_deflated'] = np.array(defl_co2.loc[:, year:year].T)
    # prop spend from yhh (105x1)
    temp['proportional_spend'] = np.array(prop_spend.loc[:, year:year])
    # total spend from yhh, deflated (1x1)
    temp['total_spend_deflated'] = np.array(defl_tot_spend[year])
    
    sda_order = list(temp.keys())
    
    # emissions
    foot = cp.copy(temp[sda_order[0]])
    for item in sda_order[1:]:
        foot = np.dot(foot, temp[item])
    footprint[year] = foot[0][0]
    
    # format to match function
    sda_vars[year] = {}
    for i in range(len(sda_order)):
        sda_vars[year][i] = temp[sda_order[i]]

# Run Analysis   
sda = {}
base_year = 2001
sda_years = cp.copy(years); sda_years.remove(base_year)
for year in years:
    sda[year] = bof.sda(sda_vars[year], sda_vars[base_year])
    sda[year].columns = ['total'] + sda_order

##############
## Save All ##
##############

equ_hhd.to_csv(outputs_filepath + 'basket_2024/equivalised_household.csv')
cm_index.to_csv(outputs_filepath + 'basket_2024/carbon_multiplier_index.csv')
basket_change.to_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change.csv')
sda.to_csv(outputs_filepath + 'basket_2024/SDA.csv')