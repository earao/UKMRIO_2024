# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:02:36 2023

@author: geolki
"""

import pandas as pd
from sys import platform
import matplotlib.pyplot as plt
import copy as cp
import seaborn as sns

# set working directory
# make different path depending on operating system

if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
outputs_filepath = wd + 'outputs/'
plots_filepath = outputs_filepath + 'basket_2024/plots/'
processed_data_filepath = wd + 'data/processed/Basket_data/'

cpi_base_year = 2010
years = list(range(2001, 2022))

# load data
equ_hhd = pd.read_csv(outputs_filepath + 'basket_2024/equivalised_household.csv', index_col=0).drop(['total_ghg'], axis=1)
equ_hhd_source = pd.read_csv(processed_data_filepath + 'compare_2001_2021_imports_domestic.csv', index_col=0)
equ_hhd_dom = pd.read_csv(outputs_filepath + 'basket_2024/equivalised_household_domestic.csv', index_col=0).drop(['total_ghg'], axis=1)
equ_hhd_imp = pd.read_csv(outputs_filepath + 'basket_2024/equivalised_household_imports.csv', index_col=0).drop(['total_ghg'], axis=1)
cm_index = pd.read_csv(outputs_filepath + 'basket_2024/carbon_multiplier_index.csv', index_col=0)
basket_change = pd.read_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change.csv', index_col=0, header=[0, 1])
basket_change_3y = pd.read_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change_3yr_avg.csv', index_col=0, header=[0, 1])
sda = pd.read_csv(outputs_filepath + 'basket_2024/SDA_mean.csv')

# Change to coicop 1
ccp1_dict = {1: 'Food', 2: 'Alcohol & tobacco',
             3: 'Clothing', 4: 'Housing',
             5: 'Furnishings', 6: 'Health', 7: 'Transport',
             8: 'Info & comms', 9: 'Recreation etc', 10: 'Education',
             11: 'Restaurants & hotels', 12: 'Miscellaneous'}

equ_hhd.columns = [int(x.split('.')[0]) for x in equ_hhd.columns]
equ_hhd_dom.columns = [int(x.split('.')[0]) for x in equ_hhd_dom.columns]
equ_hhd_imp.columns = [int(x.split('.')[0]) for x in equ_hhd_imp.columns]
equ_hhd = equ_hhd.sum(axis=1, level=0).rename(columns=ccp1_dict)
equ_hhd_dom = equ_hhd_dom.sum(axis=1, level=0).rename(columns=ccp1_dict)
equ_hhd_imp = equ_hhd_imp.sum(axis=1, level=0).rename(columns=ccp1_dict)

cm_index.index = pd.MultiIndex.from_arrays([[ccp1_dict[int(x.split('.')[0])] for x in cm_index.index], cm_index.index.tolist()])

# clean equ_hhd_source
equ_hhd_source['Product'] = [x.split(' two-adult')[0] for x in equ_hhd_source.index]
equ_hhd_source['Product'] = equ_hhd_source['Product'].str.replace('Recreation', 'Recreation etc').str.replace('etc etc', 'etc')
equ_hhd_source['Emissions'] = [x.split('two-adult ')[-1].replace('hhold ', '').replace('(', '').replace(')', '') for x in equ_hhd_source.index]
equ_hhd_source = equ_hhd_source.set_index(['Product', 'Emissions', 'Percentage change']).stack().reset_index()\
    .rename(columns={0:'Emissions (tCO2e)', 'level_3':'Year'})
equ_hhd_source.loc[equ_hhd_source['Emissions'].isin(['domestic', 'imports']) == False, 'Emissions'] = 'total'
equ_hhd_source['Year'] = [x.split('in ')[-1].split(' (')[0] for x in equ_hhd_source['Year']]

check = equ_hhd_source.set_index(['Product', 'Emissions', 'Year']).unstack(['Emissions', 'Year'])

equ_hhd_source_pct = equ_hhd_source.set_index(['Product', 'Emissions', 'Year'])[['Percentage change']].unstack('Emissions')
equ_hhd_source_pct = equ_hhd_source_pct.reset_index()


# fix years
equ_hhd.index = [int(x) for x in equ_hhd.index]

###########
## Plots ##
###########

# equivalised household

# total
equ_hhd.plot(kind='bar', stacked=True); plt.legend(bbox_to_anchor=(1,1)); 
#plt.title('Equivalised household emissions'); 
plt.ylabel('tCO2e'); plt.savefig(plots_filepath + 'Equ_hhlds.png', dpi=200, bbox_inches='tight')
plt.show();

# imports vs domestic
fig, axs = plt.subplots(ncols=2, sharey=True, figsize=(10, 5), sharex=True)
for i in range(2):
    co2 = ['domestic', 'imports'][i]
    temp = equ_hhd_source.loc[equ_hhd_source['Emissions'] == co2]
    sns.barplot(ax=axs[i], data=temp, y='Product', x='Emissions (tCO2e)', hue='Year')
    axs[i].set_title(co2.capitalize())
    axs[i].set_ylabel('')
    axs[i].legend(loc='lower right')
plt.tight_layout()
plt.savefig(plots_filepath + 'Equ_hhlds_source.png', dpi=200, bbox_inches='tight')
plt.show()


# carbon multipler index
cm_index.columns = [int(x) for x in cm_index.columns]
avg = cp.copy(cm_index)
temp = cp.copy(avg)
for year in avg.columns:
    year_list = []
    for y in [year-1, year, year+1]:
        if y in avg.columns:
            year_list.append(y)
    avg[year] = temp[year_list].mean(axis=1)
    
cm_index = avg.apply(lambda x: x/avg[2010] * 100)
cm_index = cm_index.drop(('Housing, water, electricity, gas and other fuels', '4.2.1 Imputed rentals of owner occupiers')).loc[:,2010:]
for item in list(ccp1_dict.values()):
    temp = cm_index.loc[item].T
    temp.plot(); plt.legend(bbox_to_anchor=(1,1)); plt.axhline(100, c='k');
    plt.savefig(plots_filepath + 'CM_index_' + item + '.png', dpi=200, bbox_inches='tight')
    plt.show()


# basket change
basket_change = basket_change.drop('TOTAL').mean(axis=1, level=0, skipna=True)
basket_change.columns = [int(x) for x in basket_change.columns]
basket_change = basket_change.stack().reset_index().rename(columns={'level_0':'full_label', 'level_1':'year', 0:'Pct. Change'})
basket_change['Coicop 1'] = [ccp1_dict[int(x.split('.')[0])] for x in basket_change['full_label']]

plt.show()
for item in list(ccp1_dict.values()):
    temp = basket_change.loc[(basket_change['Coicop 1'] == item)]
    sns.lineplot(data=temp, x='year', y='Pct. Change', hue='full_label'); 
    plt.legend(bbox_to_anchor=(1,1)); #plt.title=(item); 
    plt.axhline(100, c='k')
    plt.savefig(plots_filepath + 'Basket_change_' + item + '.png', dpi=200, bbox_inches='tight')
    plt.show()
    
# basket change 3 year average
basket_change_3y = basket_change_3y.drop('TOTAL').mean(axis=1, level=0, skipna=True)
basket_change_3y.columns = [int(x) for x in basket_change_3y.columns]
basket_change_3y = basket_change_3y.stack().reset_index().rename(columns={'level_0':'full_label', 'level_1':'year', 0:'Pct. Change'})
basket_change_3y['Coicop 1'] = [ccp1_dict[int(x.split('.')[0])] for x in basket_change_3y['full_label']]

plt.show()
for item in list(ccp1_dict.values()):
    temp = basket_change_3y.loc[(basket_change_3y['Coicop 1'] == item)]
    sns.lineplot(data=temp, x='year', y='Pct. Change', hue='full_label'); 
    plt.legend(bbox_to_anchor=(1,1)); #plt.title=(item); 
    plt.axhline(100, c='k')
    plt.savefig(plots_filepath + 'basket_change_3yr_avg☺_' + item + '.png', dpi=200, bbox_inches='tight')
    plt.show()


# SDA
sda.set_index('year').drop(['total', 'Unnamed: 0'], axis=1).plot(kind='bar', stacked=True); plt.legend(bbox_to_anchor=(1,1)); plt.axhline(0, c='k'); 
plt.savefig(plots_filepath + 'SDA_mean_bar.png', dpi=200, bbox_inches='tight')

sda.set_index('year').drop(['total', 'Unnamed: 0'], axis=1).plot();plt.legend(bbox_to_anchor=(1,1)); plt.axhline(0, c='k'); 
plt.savefig(plots_filepath + 'SDA_mean_line.png', dpi=200, bbox_inches='tight')
