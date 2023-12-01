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
import os
import io_functions as io

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
cpi = {}
cpi['2008-2022'] = pd.read_excel(data_filepath + 'raw/CPI/202311_consumerpriceinflationdetailedreferencetables.xlsx', sheet_name='Table 9', index_col=2, header=5)\
    .drop(['Unnamed: 0', 'Unnamed: 1'], axis=1).dropna(how='all')
    
cpi['2001-2007'] = pd.read_excel(data_filepath + 'raw/CPI/202311_consumerpriceinflationdetailedreferencetables_tcm77-4192423.xls', sheet_name='Table 9', index_col=2, header=5)\
    .drop(['Unnamed: 0', 'Unnamed: 1'], axis=1).dropna(how='all')

lookup = pd.read_csv(data_filepath + 'lookups/cpi_to_lcfs.csv').set_index('CPI')

for item in list(cpi.keys()):
    # macth lcfs columns
    cpi[item] = cpi[item].join(lookup, how='right').set_index('LCFS').mean(axis=0, level=0).fillna(0)
    # change to 2010 as base year
    cpi_comp = cpi[item][2010]
    cpi[item] = cpi[item].apply(lambda x: x / cpi_comp * 100)
cpi = cpi['2001-2007'].loc[:, :2007].join(cpi['2008-2022']).fillna(0)
cpi.index = [x.split(' ')[0] for x in cpi.index]
cpi = cpi.rename(index=idx_dict)


cpi.index = multipliers_all.index

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
total_spend = pd.DataFrame()
for year in years:
    temp = pd.DataFrame(spend[year].sum(axis=0)).T
    temp['year'] = year
    total_spend = total_spend.append(temp)
total_spend = total_spend.set_index('year').T.fillna(0)
total_spend.index = [x.split(' ')[0] for x in total_spend.index]
total_spend = total_spend.rename(index=idx_dict)

# get toal emissions per year per product
total_co2 = pd.DataFrame()
for year in years:
    temp = pd.DataFrame(hhd_ghg[year].apply(lambda x: x*hhd_ghg[year]['no_people']).sum(0).loc[multipliers_all.index]).T
    temp['year'] = year
    total_co2 = total_co2.append(temp)
total_co2 = total_co2.set_index('year').T
total_co2.index = [x.split(' ')[0] for x in total_co2.index]
total_co2 = total_co2.rename(index=idx_dict)


# proportional spend          
prop_spend = total_spend.apply(lambda x: x/x.sum())

# deflated total spend
defl_spend = pd.DataFrame(index = total_spend.index)
for year in years:
    defl_spend[year] = total_spend[year] * (cpi[year] / 100)
defl_tot_spend = defl_spend.sum(axis=0)

# deflated emission 
defl_co2 = pd.DataFrame(total_co2 / defl_spend).astype(float)
defl_co2 = defl_co2.loc[multipliers_all.index]

# population
population = pd.read_csv(data_filepath + 'raw/Population/series-281123.csv', index_col=0, header=7)
population.columns = ['population']

# make SDA variables

sda_vars = {}

for year in years[1:]:
    
    # pop (1x1)
    pop_0 = population.loc[2001, 'population']
    pop_1 = population.loc[year, 'population']
    
    # emission intensities, deflated (1x105)
    ghg_0 = defl_co2.loc[:, 2001]
    ghg_1 = defl_co2.loc[:, year]
    
    # prop spend from yhh (105x1)
    
    
    # total spend from yhh, deflated (1x1)
    

for year in range(2002,2021):
    
    ghg_0 = np.transpose(df(ghg_mults_deflated.loc[:,2001].values))
    ghg_1 = np.transpose(df(ghg_mults_deflated.loc[:,year].values))
    
    prp_0 = df(np.transpose(prop_spend_gen[2001]).values)
    prp_1 = df(np.transpose(prop_spend_gen[year]).values)
    
    exp_0 = df(np.diag(per_cap_spend.loc[:,2001]))
    exp_1 = df(np.diag(per_cap_spend.loc[:,year]))
       
    popp_0 = df(np.diag(pop_prop.loc[:,2001]))
    popp_1 = df(np.diag(pop_prop.loc[:,year]))
    
    tpop_0 = df(np.tile(total_pop.loc[:,2001],(4,1)))
    tpop_1 = df(np.tile(total_pop.loc[:,year],(4,1)))


    foot_0 = df.dot(ghg_0,prp_0)
    foot_0 = df.dot(foot_0,exp_0)
    foot_0 = df.dot(foot_0,popp_0)
    foot_0 = df.dot(foot_0,tpop_0)
    
    foot_1 = df.dot(ghg_1,prp_1)
    foot_1 = df.dot(foot_1,exp_1)
    foot_1 = df.dot(foot_1,popp_1)
    foot_1 = df.dot(foot_1,tpop_1)

    sda_0 = {}
    sda_0[0] = ghg_0
    sda_0[1] = prp_0
    sda_0[2] = exp_0
    sda_0[3] = popp_0
    sda_0[4] = tpop_0
    
    sda_1 = {}
    sda_1[0] = ghg_1
    sda_1[1] = prp_1
    sda_1[2] = exp_1
    sda_1[3] = popp_1
    sda_1[4] = tpop_1
    
    sda2[year] = io.sda(sda_1,sda_0)


# need emission intensities (deflate)
# population data
# proportioned spend


##############
## Save All ##
##############

equ_hhd.to_csv(outputs_filepath + 'basket_2024/equivalised_household.csv')
cm_index.to_csv(outputs_filepath + 'basket_2024/carbon_multiplier_index.csv')
basket_change.to_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change.csv')
