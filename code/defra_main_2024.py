#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Tue Sep 18 09:56:30 2018

This code builds the demand vector from LCFS data

@author: earao!

Requirements: Run ukmrio_main_2024.py and LCFS_data_main_2024.py and sub_national_footprints_2024_main.py first
'''

import pandas as pd
import numpy as np
import os
df = pd.DataFrame
import pickle
import LCF_functions as lcf
import demand_functions as dm
import defra_functions as defra
import ukmrio_functions as uk
from sys import platform


# set working directory
# make different path depending on operating system
if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/' 

# define filepaths
inputs_filepath = wd + 'data/model_inputs/'
results_filepath = wd + 'outputs/results_2023/'
deflator_filepath = wd + 'data/raw/ONS/ONS deflators/'
census_filepath = wd + 'data/raw/census data/'

years = np.array([int(x) for x in range(2001, 2021)])
allyears =  np.array([int(x) for x in range(1990, 2021)])

S = pickle.load(open(inputs_filepath + 'S.p', 'rb' ))
U = pickle.load(open(inputs_filepath + 'U.p', 'rb' ))
Y = pickle.load(open(inputs_filepath + 'Y.p', 'rb' ))

hhspenddata = pickle.load( open(results_filepath + 'hhspenddata.p', 'rb' ) )
meta = pickle.load( open(results_filepath + 'meta.p', 'rb' ) )

stressor_data = {}
for item in ['ghg', 'uk_ghg_direct', 'co2', 'uk_co2_direct', 'WATblu_cons', 'uk_wat_blu_cons_direct', 
             'WATgrn_cons', 'WATblu_wdrl', 'uk_wat_blu_wdrl_direct', 'mat', 'nrg', 'uk_nrg_direct',
            'bio', 'ffl', 'nmm', 'ore']:
    stressor_data[item.lower()] = pickle.load(open(results_filepath + item + '.p', 'rb' ) )

concs_dict = pd.read_excel(os.path.join(inputs_filepath, 'ONS_to_COICOP_LCF_concs.xlsx'), sheet_name=None, header = (0), index_col=0)

Y2 = lcf.convert43to41(Y, concs_dict, allyears)

Y=Y2

total_Yhh_109 = dm.make_Yhh_109_34(Y, years, meta)

coicop_exp_tot = lcf.make_totals_2023(hhspenddata, years)
coicop_exp_tot2 = lcf.convert_exp_tot_sizes(coicop_exp_tot, concs_dict, years, '456_to_105')
coicop_exp_tot3 = lcf.make_balanced_totals_2023(coicop_exp_tot2, total_Yhh_109, concs_dict, years)

yhh_wide = lcf.make_y_hh_105(Y, coicop_exp_tot3, years, concs_dict, meta)
    
newY = lcf.make_new_Y_105(Y, yhh_wide, years)

(io_deflators, cc_deflators) = uk.get_deflator_data_2023(deflator_filepath)

hhspenddata2 = lcf.convert_hhspend_sizes(hhspenddata, concs_dict, years, '456_to_105')

# import region population data
regions = ['North East', 'North West', 'Yorkshire and The Humber', 'East Midlands', 'West Midlands', 'East', 'London', 'South East', 'South West', 'Scotland', 'Wales']
regoacsyr = pickle.load(open(inputs_filepath + 'regoacsyr.p', 'rb' ))
# import region population data
regpophholdsyr = {}
for yr in years:
    reglaspend = {}
    regpophholds = {}
    for r in regions:
        ladpophholds = regoacsyr[yr][r].groupby(['LAD_code']).sum()
        regpophholds[r] = ladpophholds
    regpophholdsyr[yr] = regpophholds

y_regions = lcf.make_y_regions_2023(wd, hhspenddata2, regions, regpophholdsyr, newY, years) # need to define regions and regpophholdsyr

(S_d, U_d, Y_d, newY_d, y_regions_d) = uk.deflate_io_regions(S, U, Y, newY, y_regions, allyears, years, io_deflators, meta)

v = uk.make_v(U_d, Y_d, allyears, meta)

emptydirect = df(np.zeros((1, len(allyears))), index = ['direct'], columns = allyears)

# make sure index is alphabetical
for item in ['uk_wat_blu_cons_direct', 'uk_wat_blu_wdrl_direct']:
    idx = stressor_data[item].columns.tolist()
    stressor_data[item.replace('wat_blu', 'watblu')] = stressor_data[item][sorted(idx)]

# make lookup for indicator abbreviations
indicator_dict = {item:item for item in ['ghg', 'co2', 'nrg', 'bio', 'ffl', 'nmm', 'ore', 'mat']}
indicator_dict['watblu_cons'] = 'blc'; indicator_dict['watblu_wdrl'] = 'blw'; indicator_dict['watgrn_cons'] = 'grc'

defra_results_uk = {}
for item in list(indicator_dict.keys()):
    # define direct data
    if item in ['ghg', 'co2', 'nrg', 'watblu_cons', 'watblu_wdrl']:
        direct =  stressor_data['uk_' + item + '_direct']
    else:
        direct = emptydirect
    
    indicator = indicator_dict[item]
    print(item)
    defra_results_uk['defra_' + item + '_uk'] = defra.makeukresults2023(wd, S, U, Y, newY, coicop_exp_tot2, meta, stressor_data[item], direct, indicator, allyears, years)

defra.printdefradata('UK', results_filepath, years, list(defra_results_uk.values()), [x.split('_')[1] for x in list(defra_results_uk.keys())])
