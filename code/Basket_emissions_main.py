#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 2023

Calculate emissions for LCFS 2001-2021
"""

import Basket_emissions_functions as estimate_emissions
from sys import platform
import pickle

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

years = list(range(2001, 2022))

# load LFC data
lcfs = pickle.load(open(results_filepath + 'hhspenddata.p', 'rb'))

people = {}; hhdspend={}
for year in years:
    people[year] = lcfs[year].loc[:,:'1.1.1.1.1'].iloc[:,:-1]
    hhdspend[year] = lcfs[year].loc[:,'1.1.1.1.1':'12.7.1.1.6'].astype(float) # already multiplied by weight
  
# calculate emissions
hhd_ghg, multipliers, yhh_wide = estimate_emissions.make_footprint(hhdspend, results_filepath, model_inputs)

# save product names
idx = hhd_ghg[years[0]].columns.tolist()

# calculate emission for individual households
for year in years:
    hhd_ghg[year] = hhd_ghg[year].fillna(0).apply(lambda x: x/people[year]['weight'])
    hhd_ghg[year] = people[year].join(hhd_ghg[year])
    
# save household results
pickle.dump(hhd_ghg, open(results_filepath + 'GHG_by_hhds.p', 'wb'))
pickle.dump(multipliers, open(results_filepath + 'GHG_multipliers.p', 'wb'))
pickle.dump(yhh_wide, open(results_filepath + 'SPEND_yhh.p', 'wb'))
