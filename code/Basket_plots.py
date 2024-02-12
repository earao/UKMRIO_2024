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
import pickle
import numpy as np

# set working directory
# make different path depending on operating system

if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
outputs_filepath = wd + 'outputs/'
plots_filepath = outputs_filepath + 'basket_2024/plots/'

cpi_base_year = 2010
years = list(range(2001, 2022))

# load data
equ_hhd = pd.read_csv(outputs_filepath + 'basket_2024/equivalised_household.csv', index_col=0).drop(['total_ghg'], axis=1)
cm_index = pd.read_csv(outputs_filepath + 'basket_2024/carbon_multiplier_index.csv', index_col=0)
basket_change = pd.read_csv(outputs_filepath + 'basket_2024/basket_items_ghg_change.csv', index_col=0, header=[0, 1])
sda_all = pickle.load(open(outputs_filepath + 'basket_2024/SDA.p', 'rb'))

sda = pd.DataFrame()
for year in list(sda_all.keys()):
    temp = cp.copy(sda_all[year])
    temp['total'] = temp.loc['SDA_0', 'total']
    temp = pd.DataFrame(temp.loc[['mean', 'pct_mean']].unstack()).T.swaplevel(axis=1)
    temp['year'] = year
    sda = sda.append(temp.fillna(0))


# Change to coicop 1
ccp1_dict = {1: 'Food and non-alcoholic beverages', 2: 'Alcoholic beverages, tobacco and narcotics',
             3: 'Clothing and footwear', 4: 'Housing, water, electricity, gas and other fuels',
             5: 'Furnishings, household equipment and routine household maintenance', 6: 'Health', 7: 'Transport',
             8: 'Information and communication', 9: 'Recreation, sport and culture', 10: 'Education services',
             11: 'Restaurants and accommodation services', 12: 'Miscellaneous'}

equ_hhd.columns = [int(x.split('.')[0]) for x in equ_hhd.columns]
equ_hhd = equ_hhd.sum(axis=1, level=0).rename(columns=ccp1_dict)

cm_index.index = pd.MultiIndex.from_arrays([[ccp1_dict[int(x.split('.')[0])] for x in cm_index.index], cm_index.index.tolist()])

# fix years
equ_hhd.index = [int(x) for x in equ_hhd.index]

###########
## Plots ##
###########

# equivalised household
equ_hhd.plot(kind='bar', stacked=True); plt.legend(bbox_to_anchor=(1,1)); 
#plt.title('Equivalised household emissions'); 
plt.ylabel('tCO2e'); plt.savefig(plots_filepath + 'Equ_hhlds.png', dpi=200, bbox_inches='tight')


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

avg = cp.copy(basket_change)
temp = cp.copy(avg)
avg.iloc[:, 0] = avg.iloc[:, :1].mean(axis=1)
avg.iloc[:, -1] = avg.iloc[:, -2:].mean(axis=1)
for year in avg.columns[1:-1]:
    avg[year] = temp[[year-1, year, year+1]].mean(axis=1)
    
basket_change = avg.apply(lambda x: x/avg[2010]*100).unstack().reset_index().rename(columns={0:'Pct. Change'})
basket_change['Coicop 1'] = [ccp1_dict[int(x.split('.')[0])] for x in basket_change['full_label']]

plt.show()
for item in list(ccp1_dict.values()):
    temp = basket_change.loc[(basket_change['Coicop 1'] == item)]
    sns.lineplot(data=temp, x='level_0', y='Pct. Change', hue='full_label'); 
    plt.legend(bbox_to_anchor=(1,1)); #plt.title=(item); 
    plt.axhline(100, c='k')
    plt.savefig(plots_filepath + 'Basket_change_' + item + '.png', dpi=200, bbox_inches='tight')
    plt.show()

# SDA
sda_mean = sda.set_index('year')['mean']
sda_mean.plot(kind='line');plt.legend(bbox_to_anchor=(1,1)); plt.axhline(0, c='k'); 
plt.savefig(plots_filepath + 'SDA_mean_line.png', dpi=200, bbox_inches='tight')

sda_mean.drop('total', axis=1).plot(kind='bar', stacked=True);plt.legend(bbox_to_anchor=(1,1)); plt.axhline(0, c='k'); 
plt.savefig(plots_filepath + 'SDA_mean_bar.png', dpi=200, bbox_inches='tight')


sda_pct = sda_mean.drop('total', axis=1)
temp = cp.copy(sda_pct)
for item in temp.columns:
    temp.loc[temp[item] < 0, item] = 0
sda_pct['total'] = temp.sum(1)
sda_pct = sda_pct.apply(lambda x: x/sda_pct['total'] * 100)

sda_pct.drop('total', axis=1).plot(kind='bar', stacked=True);plt.legend(bbox_to_anchor=(1,1)); plt.axhline(0, c='k'); 
plt.axhline(100, linestyle='--', c='k'); plt.axhline(-100, linestyle='--', c='k'); 
plt.savefig(plots_filepath + 'SDA_pct_bar.png', dpi=200, bbox_inches='tight')
