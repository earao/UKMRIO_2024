#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 13:46:00 2018

This file explores the LCF original survey data from the UKDA

@author: earao
"""

import pandas as pd
import os
import pickle
from sys import platform

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


# import LCFS data
years = list(range(2001,2022))

coicop_file = data_filepath + 'processed/LCFS/LCF_20231117.csv'
coicop_lookup = pd.read_csv(coicop_file, header = 0).fillna(0)

hhspenddata = {}
for year in years:
    yr = str(year)
    
    dvhh_file = coicop_lookup.loc[coicop_lookup['Desc_full'] == 'dvhh', yr].tolist()[0]
    dvper_file = coicop_lookup.loc[coicop_lookup['Desc_full'] == 'dvper', yr].tolist()[0]
    
    dvhh = pd.read_csv(lcf_filepath + dvhh_file, sep='\t')
    dvper = pd.read_csv(lcf_filepath + dvper_file, sep='\t')
    
    dvhh.columns = [x.lower() for x in dvhh.columns]
    dvper.columns = [x.lower() for x in dvper.columns]
    
    dvhh = dvhh.set_index('case')
    dvper = dvper.set_index(['case', 'person'])
    
    # import person data
    dvper_lookup = coicop_lookup.loc[coicop_lookup['Dataset'] == 'dvper'].set_index('Desc_full')

    person_data = pd.DataFrame(index=dvper.index)
    for item in dvper_lookup.index.tolist():
        var = dvper_lookup.loc[item, yr]
        if var == 0:
            person_data[item] = 0
        else:
            person_data[item] = dvper[var]
    person_data['no_people'] = 1
    
    # make OECD modified        
    age = person_data[['age']].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    nan_list = age.loc[age['age'].isna() == True].reset_index()[['case']].drop_duplicates()['case'].tolist()
    
    if len(nan_list) > len(age.index.levels[0]) / 100:
        print('More than 1% of ages missing')
        raise SystemExit
        
    age = age.sort_values('age', ascending = False)
    age['count'] = 0.5
    age.loc[age['age'] < 14, 'count'] = 0.3
    age = age[['count']].unstack(level='person')
    age[('count', 1)] = 1
    age = age.fillna(0).sum(1)
    
    # merge with person data
    person_data = person_data.apply(lambda x: pd.to_numeric(x, errors='coerce')).sum(axis=0, level='case', skipna=True)
    person_data['OECD_mod'] = age
        
    # import househols data
    dvhh_lookup = coicop_lookup.loc[coicop_lookup['Dataset'] == 'dvhh'].set_index('Coicop_full')
    
    useful_data = pd.DataFrame(index=dvhh.index)

    exp_items = [] # collect these to multiply by weight later
    for item in dvhh_lookup.index.tolist():
        if item[0] == '0':
            desc = dvhh_lookup.loc[item, 'Desc_full']
        else:
            desc = item
            exp_items.append(desc)
        
        var = dvhh_lookup.loc[item, yr]
        if var == 0:
            useful_data[desc] = 0
        else:
            useful_data[desc] = dvhh[var]
    
    # multiply expenditure variables by weight to get UK total
    useful_data[exp_items] = useful_data[exp_items].apply(lambda x: x * useful_data['weight'])
    
    # combine person and household data
    useful_data = person_data.join(useful_data)
    
    # fix imputed rent by using home ownership (values 5, 6, 7 indicate home ownership)
    useful_data.loc[useful_data['home_ownership'].isin([5, 6, 7]) == False, 'home_ownership'] = 0
    useful_data.loc[useful_data['home_ownership'].isin([5, 6, 7]) == True, 'home_ownership'] = 1
    
    useful_data['4.2.1.1.1'] = useful_data['4.2.1.1.1'] * useful_data['home_ownership']
    
    # save to dictionary
    hhspenddata[year]  = useful_data
    
    
file = os.path.join(results_filepath, 'hhspenddata.xlsx')
writer = pd.ExcelWriter(file)

for yr in years:
    hhspenddata[yr].to_excel(writer,str(yr))
    print(str(yr) + ' saved to hhspenddata.xlsx')
writer.save()

pickle.dump(hhspenddata, open(results_filepath + "hhspenddata.p", "wb" ))
