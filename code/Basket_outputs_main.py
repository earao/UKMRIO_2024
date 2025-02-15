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
import seaborn as sns
import matplotlib.pyplot as plt
import Basket_emissions_functions as bem

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
hhd_ghg_dom = pickle.load(open(outputs_filepath + 'results_2024/GHG_by_hhds_dom.p', 'rb')) # emissions by household in survey
hhd_ghg_imp = pickle.load(open(outputs_filepath + 'results_2024/GHG_by_hhds_imp.p', 'rb')) # emissions by household in survey
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
cpi_base_year = 2010
cpi = bof.import_cpi(data_filepath, idx_dict, cpi_base_year)
cpi.loc['4.2.1 Imputed rentals of owner occupiers', :] = cpi.loc['4.1.1 Actual rentals paid by tenants', :]

###########################
## Equivalised Household ##
###########################

years = list(hhd_ghg.keys())

equ_hhd = pd.DataFrame(columns=hhd_ghg[years[0]].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].columns)
equ_hhd_dom = pd.DataFrame(columns=hhd_ghg_dom[years[0]].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].columns)
equ_hhd_imp = pd.DataFrame(columns=hhd_ghg_imp[years[0]].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].columns)

for year in years:
    temp = cp.copy(hhd_ghg[year])
    temp_d = cp.copy(hhd_ghg_dom[year])
    temp_i = cp.copy(hhd_ghg_imp[year])
    if temp['OECD_mod'].sum() == 0:
        print(str(year) + ' missing OECD modifier')
    temp['pop'] = temp['OECD_mod'] * temp['weight']
    temp.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'] = temp.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.']\
        .apply(lambda x:x*temp['weight'])
    temp_d.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'] = temp_d.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.']\
        .apply(lambda x:x*temp['weight'])
    temp_i.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'] = temp_i.loc[:, '1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.']\
        .apply(lambda x:x*temp['weight'])
    temp = temp.sum(0)
    temp_d = temp_d.sum(0)
    temp_i = temp_i.sum(0)
    pop = temp['pop']
    temp = (temp.loc['1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].astype(float) / pop.astype(float)) * 1.5
    temp['total_ghg'] = temp.sum()
    temp['year'] = year
    temp_d = (temp_d.loc['1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].astype(float) / pop.astype(float)) * 1.5
    temp_d['total_ghg'] = temp_d.sum()
    temp_d['year'] = year
    temp_i = (temp_i.loc['1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'].astype(float) / pop.astype(float)) * 1.5
    temp_i['total_ghg'] = temp_i.sum()
    temp_i['year'] = year
    equ_hhd = equ_hhd.append(pd.DataFrame(temp).T)
    equ_hhd_dom = equ_hhd_dom.append(pd.DataFrame(temp_d).T)
    equ_hhd_imp = equ_hhd_imp.append(pd.DataFrame(temp_i).T)
equ_hhd = equ_hhd.set_index('year') 
equ_hhd_dom = equ_hhd_dom.set_index('year') 
equ_hhd_imp = equ_hhd_imp.set_index('year')   

#####################
## Basket of goods ##
#####################

# import basket data
basket = pd.read_csv(data_filepath + 'processed/Basket_data/Basket_spends.csv')
basket['lcfs'] = [x.split(' ')[0] for x in basket['lcfs']]
basket = basket.set_index(['INDEX_DATE', 'lcfs', 'year']).drop(201711).mean(level=['lcfs', 'year'])

# calculate basket emissions
temp = pd.DataFrame(multipliers_all.unstack())
temp['lcfs'] = [x[1].split(' ')[0] for x in temp.index.tolist()]
temp.index.names = ['year', 'full_label']
temp = temp.set_index('lcfs', append=True); temp.columns = ['multiplier']

# get ghg emissions
basket_ghg = basket.join(temp, how='left')
basket_ghg['ghg'] = basket_ghg['spend'] * basket_ghg['multiplier'] 

basket_ghg = basket_ghg[['ghg']].unstack(level=['lcfs', 'full_label']).droplevel(axis=1, level=[0, 1])
order = []
for item in multipliers_all.index.tolist():
    if item in basket_ghg.columns:
        order.append(item)
        
# get mean ghg for each year
basket_ghg = basket_ghg.fillna(basket_ghg.mean())

# sort
basket_ghg = basket_ghg[order]
basket_ghg['TOTAL'] = basket_ghg.sum(1)
basket_ghg = basket_ghg.T

# calculate basket_change compared to 2015
basket_change = basket_ghg.apply(lambda x: (x / basket_ghg[(2015)] * 100).replace(np.inf, np.nan))

# get 3 year average
basket_change_3y = cp.copy(basket_change)
years = basket_change_3y.columns.tolist()

basket_change_3y[years[-1]] = basket_change[[years[-1], years[-1]-1]].mean(1)
for year in years:
    if year == min(years):
        basket_change_3y[year] = basket_change[[year, year+1]].mean(1)
    elif year == max(years):
        basket_change_3y[year] = basket_change[[year-1, year]].mean(1)
    else:
        basket_change_3y[year] = basket_change[[year-1, year, year+1]].mean(1)

basket_change_3y = basket_change_3y.apply(lambda x: (x / basket_change_3y[(2015)] * 100))

#########
## SDA ##
#########

# Population data: https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland

# use only household emission
# Structure: pop (1x1) * emission intensities, deflated (1x105) * prop spend from yhh (105x1) * total spend from yhh, deflated (1x1)

years = list(hhd_ghg.keys())
sda_base_year = 2001
cpi_sda = cpi.apply(lambda x: x / cpi[sda_base_year] * 100)

# population
population = pd.read_csv(data_filepath + 'raw/Population/series-281123.csv', index_col=0, header=7)
population.columns = ['population']

# get total spend to make other data (from yhh)
spend = pickle.load(open(outputs_filepath + 'results_2024/SPEND_yhh.p', 'rb'))
# get toal spend per year by product
total_spend = bof.make_total(spend, years, idx_dict)
# get toal emissions per year per product
temp = {year:hhd_ghg[year].apply(lambda x: pd.to_numeric(x, errors='coerce')*hhd_ghg[year]['weight'])[list(idx_dict.values())] 
        for year in years}
total_co2 = bof.make_total(temp, years, idx_dict)

# deflated total spend
defl_spend = pd.DataFrame(index = total_spend.index)
for year in years:
    defl_spend[year] = total_spend[year] / (cpi_sda[year] / 100)
defl_tot_spend = pd.DataFrame(defl_spend.sum(axis=0))
defl_tot_spend_percapita = defl_tot_spend.join(pd.DataFrame(population))
defl_tot_spend_percapita = defl_tot_spend_percapita[0] / defl_tot_spend_percapita['population']

# proportional spend          
prop_spend = defl_spend.apply(lambda x: x/x.sum())

# deflated emission 
defl_co2 = (total_co2 / (defl_spend)).astype(float)

# make SDA variables
sda_vars = {}
footprint = {}
for year in years:
    temp = {}
    # pop (1x1)
    temp['population'] = np.array(population.loc[year, 'population'])
    # emission intensities, deflated (1x105)
    temp['ghg_intensities_deflated'] = np.array(defl_co2.loc[:, year:year].T.fillna(0))
    # prop spend from yhh (105x1)
    temp['proportional_spend'] = np.array(prop_spend.loc[:, year:year].fillna(0))
    # total spend from yhh, deflated (1x1)
    temp['total_spend_deflated_percapita'] = np.array(defl_tot_spend_percapita[year])
    
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
for year in years:
    sda_0, sda_1 = sda_vars[sda_base_year], sda_vars[year]
    sda[year] = bof.sda(sda_1, sda_0)
    sda[year].columns = ['total'] + sda_order
    
# make summary table

sda_mean = pd.DataFrame()
for year in list(sda.keys()):
    temp = cp.copy(sda[year])
    temp['total'] = temp.loc['SDA_0', 'total']
    temp = pd.DataFrame(temp.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp['year'] = year
    sda_mean = sda_mean.append(temp.fillna(0))


# Production SDA

# Structure - for domestic consumption: UK industrty emission intensities, deflated (1x112) * domestic section of Leontief inverse, deflated (112x112) * UK domestic final demand, deflated (112x1)
# Structure - for reimport : UK industrty emission intensities, deflated (1x112) * domestic section of Leontief inverse, deflated (112x1562) * UK final demand, deflated (1562x1)
# Structure - for export : UK industrty emission intensities, deflated (1x112) * domestic section of Leontief inverse, deflated (112x1680) * row final demand, deflated (112x1)
# Structure - for export : UK industrty emission intensities, deflated (1x112) * domestic section of Leontief inverse, deflated (112x1680) * row final demand, deflated (1562x1)

# load meta data from [UKMRIO]
meta = pickle.load(open(outputs_filepath + 'results_2024/meta.p', "rb" ))
   
ukmrio = {}; #means = {}
for data in ['ghg', 'uk_ghg_direct', 'S', 'U', 'Y']:
    ukmrio[data] = pickle.load(open(outputs_filepath + 'results_2024/' + data + '.p', "rb" ))
   
 # create year lists
years = list(hhd_ghg.keys())

deflator_filepath = wd + 'data/raw/ONS/ONS deflators/'
S_d, U_d, Y_d = bof.import_deflated_MRIO(deflator_filepath,ukmrio,years) 
L_d = {}
x_d = {}
for yr in years:
    Z = bem.make_Z_from_S_U(S_d[yr], U_d[yr])
    bigY = np.zeros(shape = [np.size(Y_d[yr], 0)*2, np.size(Y_d[yr], 1)])
    bigY[np.size(Y_d[yr], 0):np.size(Y_d[yr], 0)*2, 0:] = Y_d[yr] 
    x_d[yr] = bem.make_x(Z, bigY)
    L_d[yr] = bem.make_L(Z, x_d[yr])
    
# make SDA variables
sda_vars1 = {}
sda_vars2 = {}
sda_vars3 = {}
sda_vars4 = {}
footprint = {}
for yr in years:
    temp1 = {}
    temp2 = {}
    temp3 = {}
    temp4 = {}
    # emission intensities, deflated (1x112)
    temp1['ghg_intensity_deflated'] = np.array(np.sum(ukmrio['ghg'][yr].iloc[0:112],1)/x_d[yr][0:112])
    # domestic Leontief, deflated (112x112)
    temp1['domestic_leontief_deflated'] = np.array(L_d[yr][0:112, 1680:1792])
    # domestic final demand, deflated (112x1)
    temp1['domestic_final_demand'] = np.array(np.sum(Y_d[yr].iloc[0:112,0:42],1))
    # emission intensities, deflated (1x112)
    temp2['ghg_intensity_deflated'] = np.array(np.sum(ukmrio['ghg'][yr].iloc[0:112],1)/x_d[yr][0:112])
    # domestic Leontief, deflated (112x112)
    temp2['row_leontief_deflated'] = np.array(L_d[yr][0:112, 1792:3360])
    # domestic final demand, deflated (112x1)
    temp2['domestic_final_demand_of_row'] = np.array(np.sum(Y_d[yr].iloc[112:1680,0:42],1))
    # emission intensities, deflated (1x112)
    temp3['ghg_intensity_deflated'] = np.array(np.sum(ukmrio['ghg'][yr].iloc[0:112],1)/x_d[yr][0:112])
    # domestic Leontief, deflated (112x1680)
    temp3['domestic_leontief_deflated'] = np.array(L_d[yr][0:112, 1680:1792])
    # row final demand, deflated (1680x1)
    temp3['row_final_demand'] = np.array(np.sum(Y_d[yr].iloc[0:112,42:43],1))
    # emission intensities, deflated (1x112)
    temp4['ghg_intensity_deflated'] = np.array(np.sum(ukmrio['ghg'][yr].iloc[0:112],1)/x_d[yr][0:112])
    # domestic Leontief, deflated (112x1680)
    temp4['domestic_leontief_deflated'] = np.array(L_d[yr][0:112, 1792:3350])
    # row final demand, deflated (1680x1)
    temp4['row_final_demand'] = np.array(np.sum(Y_d[yr].iloc[112:1680,42:43],1))
    
    sda_order1 = list(temp1.keys())
    sda_order2 = list(temp2.keys())
    sda_order3 = list(temp3.keys())
    sda_order4 = list(temp4.keys())
    
    # format to match function
    sda_vars1[yr] = {}
    sda_vars2[yr] = {}
    sda_vars3[yr] = {}
    sda_vars4[yr] = {}
    for i in range(len(sda_order1)):
        sda_vars1[yr][i] = temp1[sda_order1[i]]
        sda_vars2[yr][i] = temp2[sda_order2[i]]
        sda_vars3[yr][i] = temp3[sda_order3[i]]
        sda_vars4[yr][i] = temp3[sda_order4[i]]

# Run Analysis   
sda1 = {}
sda2 = {}
sda3 = {}
sda4 = {}
for yr in years:
    sda_0, sda_1 = sda_vars1[2001], sda_vars1[yr]
    sda1[yr] = bof.sda(sda_1, sda_0)
    sda1[yr].columns = ['total'] + sda_order1
    
    sda_2, sda_3 = sda_vars2[2001], sda_vars2[yr]
    sda2[yr] = bof.sda(sda_3, sda_2)
    sda2[yr].columns = ['total'] + sda_order2
    
    sda_4, sda_5 = sda_vars3[2001], sda_vars3[yr]
    sda3[yr] = bof.sda(sda_5, sda_4)
    sda3[yr].columns = ['total'] + sda_order3
    
    sda_6, sda_7 = sda_vars4[2001], sda_vars4[yr]
    sda4[yr] = bof.sda(sda_7, sda_6)
    sda4[yr].columns = ['total'] + sda_order4

# make summary table

sda_mean1 = pd.DataFrame()
sda_mean2 = pd.DataFrame()
sda_mean3 = pd.DataFrame()
sda_mean4 = pd.DataFrame()
for yr in list(sda1.keys()):
    temp1 = cp.copy(sda1[yr])
    temp1['total'] = temp1.loc['SDA_0', 'total']
    temp1 = pd.DataFrame(temp1.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp1['year'] = yr
    sda_mean1 = sda_mean1.append(temp1.fillna(0))
    
    temp2 = cp.copy(sda2[yr])
    temp2['total'] = temp2.loc['SDA_0', 'total']
    temp2 = pd.DataFrame(temp2.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp2['year'] = yr
    sda_mean2 = sda_mean2.append(temp2.fillna(0))
    
    temp3 = cp.copy(sda3[yr])
    temp3['total'] = temp3.loc['SDA_0', 'total']
    temp3 = pd.DataFrame(temp3.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp3['year'] = yr
    sda_mean3 = sda_mean3.append(temp3.fillna(0))
    
    temp4 = cp.copy(sda4[yr])
    temp4['total'] = temp4.loc['SDA_0', 'total']
    temp4 = pd.DataFrame(temp4.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp4['year'] = yr
    sda_mean4 = sda_mean4.append(temp3.fillna(0))

# Comsumption SDA anne edit

# Structure - UK industry emission intensities, deflated (1x112) * d_d section of Leontief inverse, deflated (112x112) * UK domestic final demand, deflated (112x1)
# Structure - UK industry emission intensities, deflated (1x112) * d_i section of Leontief inverse, deflated (112x1562) * UK row demand, deflated (1562x1)
# Structure - RoW industry emission intensities, deflated (1x1568) * i_d  section of Leontief inverse, deflated (1562x112) * UK domestic  final demand, deflated (112x1)
# Structure - RoW industrty emission intensities, deflated (1x1568) * i_i section of Leontief inverse, deflated (1562x1562) * UK row final demand, deflated (1562x1)

# load meta data from [UKMRIO]
meta = pickle.load(open(outputs_filepath + 'results_2024/meta.p', "rb" ))
   
ukmrio = {}; #means = {}
for data in ['ghg', 'uk_ghg_direct', 'S', 'U', 'Y']:
    ukmrio[data] = pickle.load(open(outputs_filepath + 'results_2024/' + data + '.p', "rb" ))
   
 # create year lists
years = list(hhd_ghg.keys())

deflator_filepath = wd + 'data/raw/ONS/ONS deflators/'
S_d, U_d, Y_d = bof.import_deflated_MRIO(deflator_filepath,ukmrio,years) 
eL_d = {}
for yr in years:
    Z = bem.make_Z_from_S_U(S_d[yr], U_d[yr])
    bigY = np.zeros(shape = [np.size(Y_d[yr], 0)*2, np.size(Y_d[yr], 1)])
    bigY[np.size(Y_d[yr], 0):np.size(Y_d[yr], 0)*2, 0:] = Y_d[yr] 
    x_d = bem.make_x(Z, bigY)
    L_d = bem.make_L(Z, x_d)
    bigstressor = np.zeros(shape = [np.size(Y_d[yr], 0)*2, 1])
    bigstressor[0:np.size(Y_d[yr], 0), :] = ukmrio['ghg'][yr]
    e = np.diag(np.sum(bigstressor, 1)/x_d)
    eL_d[yr] = np.dot(e, L_d)
    
# make SDA variables
sda_vars1 = {}
sda_vars2 = {}
sda_vars3 = {}
sda_vars4 = {}
footprint = {}
for yr in years:
    temp1 = {}
    temp2 = {}
    temp3 = {}
    temp4 = {}
    # emission intensities, deflated (1x112)
    temp1['domestic_domestic_emissions_deflated'] = np.sum(np.array(eL_d[yr][0:112, 1680:1792]),0)
    # domestic final demand, deflated (112x1)
    temp1['domestic_final_demand'] = np.array(np.sum(Y_d[yr].iloc[0:112,0:42],1))
    # emission intensities, deflated (1x112)
    temp2['domestic_row_emissions_deflated'] = np.sum(np.array(eL_d[yr][0:112, 1792:3360]),0)
    # domestic final demand, deflated (112x1)
    temp2['domestic_final_demand_of_row'] = np.array(np.sum(Y_d[yr].iloc[112:1680,0:42],1))
    # emission intensities, deflated (1x112)
    temp3['row_domestic_emissions_deflated'] = np.sum(np.array(eL_d[yr][112:1680, 1680:1792]),0)
    # row final demand, deflated (1680x1)
    temp3['domestic_final_demand'] = np.array(np.sum(Y_d[yr].iloc[0:112,0:42],1))
    # emission intensities, deflated (1x112)
    temp4['row_row_emissions_deflated'] = np.sum(np.array(eL_d[yr][112:1680, 1792:3360]),0)
    # row final demand, deflated (1680x1)
    temp4['domestic_final_demand_of_row'] = np.array(np.sum(Y_d[yr].iloc[112:1680,0:42],1))
    
    sda_order1 = list(temp1.keys())
    sda_order2 = list(temp2.keys())
    sda_order3 = list(temp3.keys())
    sda_order4 = list(temp4.keys())
    
    # format to match function
    sda_vars1[yr] = {}
    sda_vars2[yr] = {}
    sda_vars3[yr] = {}
    sda_vars4[yr] = {}
    for i in range(len(sda_order1)):
        sda_vars1[yr][i] = temp1[sda_order1[i]]
        sda_vars2[yr][i] = temp2[sda_order2[i]]
        sda_vars3[yr][i] = temp3[sda_order3[i]]
        sda_vars4[yr][i] = temp4[sda_order4[i]]

# Run Analysis   
sda1 = {}
sda2 = {}
sda3 = {}
sda4 = {}
for yr in years:
    sda_0, sda_1 = sda_vars1[2001], sda_vars1[yr]
    sda1[yr] = bof.sda(sda_1, sda_0)
    sda1[yr].columns = ['total'] + sda_order1
    
    sda_2, sda_3 = sda_vars2[2001], sda_vars2[yr]
    sda2[yr] = bof.sda(sda_3, sda_2)
    sda2[yr].columns = ['total'] + sda_order2
    
    sda_4, sda_5 = sda_vars3[2001], sda_vars3[yr]
    sda3[yr] = bof.sda(sda_5, sda_4)
    sda3[yr].columns = ['total'] + sda_order3
    
    sda_6, sda_7 = sda_vars4[2001], sda_vars4[yr]
    sda4[yr] = bof.sda(sda_7, sda_6)
    sda4[yr].columns = ['total'] + sda_order4

# make summary table

c_sda_mean1 = pd.DataFrame()
c_sda_mean2 = pd.DataFrame()
c_sda_mean3 = pd.DataFrame()
c_sda_mean4 = pd.DataFrame()
for yr in list(sda1.keys()):
    temp1 = cp.copy(sda1[yr])
    temp1['total'] = temp1.loc['SDA_0', 'total']
    temp1 = pd.DataFrame(temp1.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp1['year'] = yr
    c_sda_mean1 = c_sda_mean1.append(temp1.fillna(0))
    
    temp2 = cp.copy(sda2[yr])
    temp2['total'] = temp2.loc['SDA_0', 'total']
    temp2 = pd.DataFrame(temp2.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp2['year'] = yr
    c_sda_mean2 = c_sda_mean2.append(temp2.fillna(0))
    
    temp3 = cp.copy(sda3[yr])
    temp3['total'] = temp3.loc['SDA_0', 'total']
    temp3 = pd.DataFrame(temp3.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp3['year'] = yr
    c_sda_mean3 = c_sda_mean3.append(temp3.fillna(0))
    
    temp4 = cp.copy(sda4[yr])
    temp4['total'] = temp4.loc['SDA_0', 'total']
    temp4 = pd.DataFrame(temp4.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
    temp4['year'] = yr
    c_sda_mean4 = c_sda_mean4.append(temp4.fillna(0))
    
# Comsumption 12 SDA anne edit

# Structure - UK industry emission intensities, deflated (1x112) * d_d section of Leontief inverse, deflated (112x112) * UK domestic final demand, deflated (112x1)
# Structure - UK industry emission intensities, deflated (1x112) * d_i section of Leontief inverse, deflated (112x1562) * UK row demand, deflated (1562x1)
# Structure - RoW industry emission intensities, deflated (1x1568) * i_d  section of Leontief inverse, deflated (1562x112) * UK domestic  final demand, deflated (112x1)
# Structure - RoW industrty emission intensities, deflated (1x1568) * i_i section of Leontief inverse, deflated (1562x1562) * UK row final demand, deflated (1562x1)

# load meta data from [UKMRIO]
meta = pickle.load(open(outputs_filepath + 'results_2024/meta.p', "rb" ))
   
ukmrio = {}; #means = {}
for data in ['ghg', 'uk_ghg_direct', 'S', 'U', 'Y']:
    ukmrio[data] = pickle.load(open(outputs_filepath + 'results_2024/' + data + '.p', "rb" ))
   
 # create year lists
years = list(hhd_ghg.keys())

deflator_filepath = wd + 'data/raw/ONS/ONS deflators/'
S_d, U_d, Y_d = bof.import_deflated_MRIO(deflator_filepath,ukmrio,years) 
eL_d = {}
Y_d_12 = {}
for yr in years:
    temp = np.zeros((1680,12))
    Z = bem.make_Z_from_S_U(S_d[yr], U_d[yr])
    bigY = np.zeros(shape = [np.size(Y_d[yr], 0)*2, np.size(Y_d[yr], 1)])
    bigY[np.size(Y_d[yr], 0):np.size(Y_d[yr], 0)*2, 0:] = Y_d[yr] 
    x_d = bem.make_x(Z, bigY)
    L_d = bem.make_L(Z, x_d)
    bigstressor = np.zeros(shape = [np.size(Y_d[yr], 0)*2, 1])
    bigstressor[0:np.size(Y_d[yr], 0), :] = ukmrio['ghg'][yr]
    e = np.diag(np.sum(bigstressor, 1)/x_d)
    eL_d[yr] = np.dot(e, L_d)    
    temp[:,0]  = np.sum(Y_d[yr].iloc[:,0:2],1) 
    temp[:,1]  = np.sum(Y_d[yr].iloc[:,2:5],1) 
    temp[:,2]  = np.sum(Y_d[yr].iloc[:,5:7],1) 
    temp[:,3]  = np.sum(Y_d[yr].iloc[:,7:12],1) 
    temp[:,4]  = np.sum(Y_d[yr].iloc[:,12:18],1) 
    temp[:,5]  = np.sum(Y_d[yr].iloc[:,18:21],1) 
    temp[:,6]  = np.sum(Y_d[yr].iloc[:,21:24],1) 
    temp[:,7]  = np.sum(Y_d[yr].iloc[:,24:27],1) 
    temp[:,8]  = np.sum(Y_d[yr].iloc[:,27:33],1) 
    temp[:,9]  = np.sum(Y_d[yr].iloc[:,33:34],1) 
    temp[:,10]  = np.sum(Y_d[yr].iloc[:,34:35],1) 
    temp[:,11]  = np.sum(Y_d[yr].iloc[:,35:36],1) 
    Y_d_12[yr] = temp  
    
# make SDA variables
sda_vars1_c = {}
sda_vars2_c = {}
sda_vars3_c = {}
sda_vars4_c = {}
for yr in years:
    sda_vars1 = {}
    sda_vars2 = {}
    sda_vars3 = {}
    sda_vars4 = {}
    for c in range(0,12):   
        temp1 = {}
        temp2 = {}
        temp3 = {}
        temp4 = {}
        # emission intensities, deflated (1x112)
        temp1['domestic_domestic_emissions_deflated'] = np.sum(np.array(eL_d[yr][0:112, 1680:1792]),0)
        # domestic final demand, deflated (112x1)
        temp1['domestic_final_demand'] = np.array(Y_d_12[yr][0:112,c])
        # emission intensities, deflated (1x112)
        temp2['domestic_row_emissions_deflated'] = np.sum(np.array(eL_d[yr][0:112, 1792:3360]),0)
        # domestic final demand, deflated (112x1)
        temp2['domestic_final_demand_of_row'] = np.array(Y_d_12[yr][112:1680,c])
        # emission intensities, deflated (1x112)
        temp3['row_domestic_emissions_deflated'] = np.sum(np.array(eL_d[yr][112:1680, 1680:1792]),0)
        # row final demand, deflated (1680x1)
        temp3['domestic_final_demand'] = np.array(Y_d_12[yr][0:112,c])
        # emission intensities, deflated (1x112)
        temp4['row_row_emissions_deflated'] = np.sum(np.array(eL_d[yr][112:1680, 1792:3360]),0)
        # row final demand, deflated (1680x1)
        temp4['domestic_final_demand_of_row'] = np.array(Y_d_12[yr][112:1680,c])
        
        sda_order1 = list(temp1.keys())
        sda_order2 = list(temp2.keys())
        sda_order3 = list(temp3.keys())
        sda_order4 = list(temp4.keys())
        
        # format to match function
        sda_vars1[c] = {}
        sda_vars2[c] = {}
        sda_vars3[c] = {}
        sda_vars4[c] = {}
        for i in range(len(sda_order1)):
            sda_vars1[c][i] = temp1[sda_order1[i]]
            sda_vars2[c][i] = temp2[sda_order2[i]]
            sda_vars3[c][i] = temp3[sda_order3[i]]
            sda_vars4[c][i] = temp4[sda_order4[i]]
    sda_vars1_c[yr] = sda_vars1
    sda_vars2_c[yr] = sda_vars2
    sda_vars3_c[yr] = sda_vars3
    sda_vars4_c[yr] = sda_vars4
    
# Run Analysis 
sda1_c = {}
sda2_c = {}
sda3_c = {}
sda4_c = {}  
for c in range (0,12):
    sda1 = {}
    sda2 = {}
    sda3 = {}
    sda4 = {}
    for yr in years:
        sda_0, sda_1 = sda_vars1_c[2001][c], sda_vars1_c[yr][c]
        sda1[yr] = bof.sda(sda_1, sda_0)
        sda1[yr].columns = ['total'] + sda_order1
        
        sda_2, sda_3 = sda_vars2_c[2001][c], sda_vars2_c[yr][c]
        sda2[yr] = bof.sda(sda_3, sda_2)
        sda2[yr].columns = ['total'] + sda_order2
        
        sda_4, sda_5 = sda_vars3_c[2001][c], sda_vars3_c[yr][c]
        sda3[yr] = bof.sda(sda_5, sda_4)
        sda3[yr].columns = ['total'] + sda_order3
        
        sda_6, sda_7 = sda_vars4_c[2001][c], sda_vars4_c[yr][c]
        sda4[yr] = bof.sda(sda_7, sda_6)
        sda4[yr].columns = ['total'] + sda_order4
    sda1_c[c] = sda1
    sda2_c[c] = sda2
    sda3_c[c] = sda3
    sda4_c[c] = sda4

# make summary table

c_sda_mean1_c = {}
c_sda_mean2_c = {}
c_sda_mean3_c = {}
c_sda_mean4_c = {}

for c in range(0,12):
    c_sda_mean1 = pd.DataFrame()
    c_sda_mean2 = pd.DataFrame()
    c_sda_mean3 = pd.DataFrame()
    c_sda_mean4 = pd.DataFrame()
    for yr in list(sda1.keys()):
        temp1 = cp.copy(sda1_c[c][yr])
        temp1['total'] = temp1.loc['SDA_0', 'total']
        temp1 = pd.DataFrame(temp1.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
        temp1['year'] = yr
        c_sda_mean1 = c_sda_mean1.append(temp1.fillna(0))
        
        temp2 = cp.copy(sda2_c[c][yr])
        temp2['total'] = temp2.loc['SDA_0', 'total']
        temp2 = pd.DataFrame(temp2.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
        temp2['year'] = yr
        c_sda_mean2 = c_sda_mean2.append(temp2.fillna(0))
        
        temp3 = cp.copy(sda3_c[c][yr])
        temp3['total'] = temp3.loc['SDA_0', 'total']
        temp3 = pd.DataFrame(temp3.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
        temp3['year'] = yr
        c_sda_mean3 = c_sda_mean3.append(temp3.fillna(0))
        
        temp4 = cp.copy(sda4_c[c][yr])
        temp4['total'] = temp4.loc['SDA_0', 'total']
        temp4 = pd.DataFrame(temp4.loc[['mean']].unstack()).T.swaplevel(axis=1).droplevel(axis=1, level=0)
        temp4['year'] = yr
        c_sda_mean4 = c_sda_mean4.append(temp4.fillna(0))
    c_sda_mean1_c[c] = c_sda_mean1
    c_sda_mean2_c[c] = c_sda_mean2
    c_sda_mean3_c[c] = c_sda_mean3
    c_sda_mean4_c[c] = c_sda_mean4

#############################
## Carbon multiplier index ##
#############################

years = list(hhd_ghg.keys())
cm_base_year = 2015
cpi_cm = cpi.apply(lambda x: x / cpi[cm_base_year] * 100)

# get total spend to make other data (from yhh)
spend = pickle.load(open(outputs_filepath + 'results_2024/SPEND_yhh.p', 'rb'))
# get toal spend per year by product
total_spend = bof.make_total(spend, years, idx_dict)
# get toal emissions per year per product
temp = {year:hhd_ghg[year].apply(lambda x: pd.to_numeric(x, errors='coerce')*hhd_ghg[year]['weight'])[list(idx_dict.values())] 
        for year in years}
total_co2 = bof.make_total(temp, years, idx_dict)

# save order
order = total_spend.index.tolist()

# calculate index by getting multiplier and deflating it
defl_cm = total_co2 / total_spend * cpi_cm
defl_cm = defl_cm.drop(2022, axis=1)

# fill in missing values with mean from either side
temp = defl_cm[defl_cm.isna().any(axis=1)]
# fill 2001 with nearest value
for year in years[1:]:
    temp.loc[temp[years[0]].isna() == True, years[0]] = temp[year]
# fill 2021 with nearest value
years.reverse()
for year in years[1:]:
    temp.loc[temp[years[0]].isna() == True, years[0]] = temp[year]
years.reverse()
# fill in other missing values with gradual mean to either side
temp = temp.T
for item in temp.columns:
    fill = []
    
    temp2 = temp[[item]]
    missing = temp2[temp2.isna().any(axis=1)].index.tolist()
    
    for year in years:
        if year in missing:
            for i in list(range(years[0], year)):
                if i not in missing:
                    start = i
            for i in list(reversed(range(year+1, years[-1]+1))):
                if i not in missing:
                    end = i
            fill.append((temp.loc[[start, end], item]).mean())
            print(year, (temp.loc[[start, end], item]).mean())
        else:
            fill.append(temp.loc[year, item])
            print(year, temp.loc[year, item])
    temp[item] = fill

defl_cm = defl_cm.drop(temp.columns.tolist()).append(temp.T).loc[order]

# adjust base year to 100%
cm_index = defl_cm.apply(lambda x: x/defl_cm[cm_base_year]*100)

# plot for presentation
plot_data = pd.DataFrame(multipliers_all.apply(lambda x: x/multipliers_all[cm_base_year]*100).stack())\
    .rename(columns={0:'Original'}).join(pd.DataFrame(cm_index.stack()).rename(columns={0:'Deflated (CM Index)'}))\
        .stack().reset_index().rename(columns={'level_0':'ccp3', 'level_1':'Year', 'level_2':'Type', 0:'Multiplier value rescaled (2015 = 100)'})


ccp_list = ['7.2.2 Fuels and lubricants for personal transport equipment', '4.5.1 Electricity'] 
            #'7.3.1_2 Passenger transport by railway and road', '4.5.2 Gas', 

from matplotlib.ticker import MaxNLocator
for item in ccp_list:
    temp = plot_data.loc[plot_data['ccp3'] == item]
    ax = plt.figure().gca()
    sns.lineplot(data=temp, x='Year', y='Multiplier value rescaled (2015 = 100)', hue='Type')
    plt.title(item)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.axhline(100, c='k', linestyle=':')
    plt.savefig(outputs_filepath + 'basket_2024/plots/cm_index_plot_' + item  + '.png', bbox_inches='tight')
    plt.show()



##############
## Save All ##
##############

equ_hhd.to_csv(outputs_filepath + 'basket_2024/equivalised_household.csv')
equ_hhd_dom.to_csv(outputs_filepath + 'basket_2024/equivalised_household_domestic.csv')
equ_hhd_imp.to_csv(outputs_filepath + 'basket_2024/equivalised_household_imports.csv')
cm_index.to_csv(outputs_filepath + 'basket_2024/carbon_multiplier_index.csv')
basket_change.to_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change.csv')
basket_change_3y.to_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change_3yr_avg.csv')
sda_mean.to_csv(outputs_filepath + 'basket_2024/SDA_mean.csv')

writer = pd.ExcelWriter(outputs_filepath + 'basket_2024/SDA.xlsx')
for year in sda.keys():
    sda[year].to_excel(writer, sheet_name=str(year))
writer.save()
pickle.dump(sda, open(outputs_filepath + 'basket_2024/SDA.p', 'wb' ) )