# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 16:04:10 2021

@author: earao


Requirements: Run ukmrio_main_2023.py and LCFS_data_main_2023.py first
"""

import pandas as pd
import numpy as np
import os
df = pd.DataFrame
import LCF_functions as lcf
import pickle
import defra_functions as defra
import census_functions as cs
import ukmrio_functions as uk
import copy
from sys import platform

# set working directory
# make different path depending on operating system
if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
inputs_filepath = wd + 'data/model_inputs/'
results_filepath = wd + 'outputs/results_2024/'
lcf_filepath =  wd +'data/raw/LCFS/'
deflator_filepath =  wd +'data/raw/ONS/ONS deflators/'
census_filepath =  wd +'data/raw/census data/'
energy_filepath = wd + 'data/processed/uk energy/'


years = np.array([int(x) for x in range(2001,2022)])
allyears =  np.array([int(x) for x in range(1990,2022)])
oac_none_years = np.array([int(x) for x in range(2001,2007)])
oac_2001_years = np.array([int(x) for x in range(2007,2014)])
oac_2011_years = np.array([int(x) for x in range(2014,2022)])

S = pickle.load( open(results_filepath + "S.p", "rb" ) )
U = pickle.load( open(results_filepath + "U.p", "rb" ) )
Y = pickle.load( open(results_filepath + "Y.p", "rb" ) )
hhspenddata = pickle.load( open(results_filepath + "hhspenddata.p", "rb" ) )
meta = pickle.load( open(results_filepath + "meta.p", "rb" ) )

ghg = pickle.load( open(results_filepath + "ghg.p", "rb" ) )
uk_ghg_direct = pickle.load( open(results_filepath + "uk_ghg_direct.p", "rb" ) )

co2 = pickle.load( open(results_filepath + "co2.p", "rb" ) )
uk_co2_direct = pickle.load( open(results_filepath + "uk_co2_direct.p", "rb" ) )

hhspenddata2 = copy.deepcopy(hhspenddata)
hhspenddata2 = lcf.removeoutliers(hhspenddata2,years)
        
concs_dict = pd.read_excel(os.path.join(inputs_filepath, 'ONS_to_COICOP_LCF_concs_2024.xlsx'), sheet_name=None, header = 0, index_col=0)

Y2 = lcf.convert43to41(Y,concs_dict,allyears)

Y=Y2

total_Yhh_109 = lcf.make_Yhh_109_34(Y,years,meta) 

coicop_exp_tot = lcf.make_totals(hhspenddata2,years)

coicop_exp_tot2 = lcf.convert_exp_tot_sizes(coicop_exp_tot,concs_dict,years,'456_to_105')
coicop_exp_tot3 = lcf.make_balanced_totals(coicop_exp_tot2,total_Yhh_109,concs_dict,years)

yhh_wide = lcf.make_y_hh_105(Y,coicop_exp_tot3,years,concs_dict,meta)
newY = lcf.make_new_Y_105(Y,yhh_wide,years)

(io_deflators,cc_deflators) = uk.get_deflator_data(deflator_filepath)
cc_deflators = np.transpose(cc_deflators.loc[years,:])

hhspenddata3 = lcf.convert_hhspend_sizes(hhspenddata2,concs_dict,years,'456_to_105')

regions = ['North East','North West','Yorkshire and The Humber','East Midlands','West Midlands','East','London','South East','South West','Wales','Scotland']
regions_lc = ['north_east','north_west','yorkshire','east_midlands','west_midlands','east','london','south_east','south_west','wales','scotland']
region_dict = dict(zip(regions_lc, list(range(1,13))))

oa_lookup01 = cs.make_01_lookup(census_filepath) 
(oa_lookup11,s11) = cs.make_11_lookup(census_filepath) 
oa_lookup21 = cs.make_21_lookup(census_filepath,s11)

regoacsyr = lcf.make_pop_hhold_by_oac_region_year(oa_lookup01,oa_lookup11,oa_lookup21,census_filepath,regions) # save?
pickle.dump(regoacsyr, open(results_filepath + "regoacsyr.p", "wb" ))

# turn OAC classes into string variable
for yr in list(hhspenddata3.keys()):
    hhspenddata3[yr]['OA class 1'] = hhspenddata3[yr]['OA class 1'].astype(str)
(oacyrspends,oacyrmeta) = lcf.makeoacspends(hhspenddata3,oac_none_years,oac_2001_years,oac_2011_years)

(reglaspendyr,regpophholdsyr) = lcf.make_la_spends_pop_by_region_year(regoacsyr,oacyrspends,regions,years) 
reglaspendyr = lcf.correct_reglaspendyr_zero(reglaspendyr,regions,years)                 
reglaspropyr = lcf.make_la_spend_props_by_region_year(reglaspendyr,regions,years)    
reglaspropyr = lcf.correct_la_spend_props_gas_elec(reglaspropyr,regions,regions_lc,years,wd)
reglaspropyr = lcf.add_pop_factors(reglaspropyr,regpophholdsyr,newY,regions,years)

y_regions = lcf.make_y_regions(wd,hhspenddata3,regions,regpophholdsyr,newY,years)

defra_ghg_uk = defra.makeukresults(wd,S,U,Y,newY,coicop_exp_tot2,meta,ghg,uk_ghg_direct,'ghg',allyears,years)
defra_co2_uk = defra.makeukresults(wd,S,U,Y,newY,coicop_exp_tot2,meta,co2,uk_co2_direct,'co2',allyears,years)

y_regional = {}
for region in regions_lc:
    y_regional[region] = lcf.make_y_regional(region, region_dict, y_regions, years)
        
(y_England,y_N_Ireland,y_Scotland,y_Wales) = lcf.make_y_countries(y_regions,years)

defra_ghg_reg = {}
for country in ['England', 'Scotland', 'Wales', 'Northern Ireland']:
    y_country = eval('y_' + country.replace('Northern ', 'N_'))
    defra_ghg_reg[country] = defra.makeregionresults(S,U,Y,newY,y_country,meta,ghg,defra_ghg_uk['direct'],'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,country,cc_deflators)

reg_england = ['North East','North West','Yorkshire and The Humber','East Midlands','West Midlands','East','London','South East','South West']
for reg in reg_england:
    reg_lower = reg.replace(' and The Humber', '').replace(' ', '_').lower()
    defra_ghg_reg[reg] = defra.makeregionresults(S,U,Y,newY,y_regional[reg_lower],meta,ghg,defra_ghg_uk['direct'],
                                                           'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,reg,cc_deflators)
    
defra_ghg_reg['Wales'] = defra.makeregionresults(S,U,Y,newY,y_regional['wales'],meta,ghg,defra_ghg_uk['direct'],
                                                       'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'Wales',cc_deflators)
defra_ghg_reg['Northern Ireland'] = defra.makeregionresults(S,U,Y,newY,y_N_Ireland,meta,ghg,defra_ghg_uk['direct'],
                                                       'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'Northern Ireland',cc_deflators)
defra_ghg_reg['Scotland'] = defra.makeregionresults(S,U,Y,newY,y_regional['scotland'],meta,ghg,defra_ghg_uk['direct'],
                                                       'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'Scotland',cc_deflators)

defra_co2_reg = {}
for country in ['England', 'Scotland', 'Wales', 'Northern Ireland']:
    y_country = eval('y_' + country.replace('Northern ', 'N_'))
    defra_co2_reg[country] = defra.makeregionresults(S,U,Y,newY,y_country,meta,co2,defra_co2_uk['direct'],'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,country,cc_deflators)

reg_england = ['North East','North West','Yorkshire and The Humber','East Midlands','West Midlands','East','London','South East','South West']
for reg in reg_england:
    reg_lower = reg.replace(' and The Humber', '').replace(' ', '_').lower()
    defra_co2_reg[reg] = defra.makeregionresults(S,U,Y,newY,y_regional[reg_lower],meta,co2,defra_co2_uk['direct'],
                                                           'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,reg,cc_deflators)
    
defra_co2_reg['Wales'] = defra.makeregionresults(S,U,Y,newY,y_regional['wales'],meta,co2,defra_co2_uk['direct'],
                                                       'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,'Wales',cc_deflators)
defra_co2_reg['Northern Ireland'] = defra.makeregionresults(S,U,Y,newY,y_N_Ireland,meta,co2,defra_co2_uk['direct'],
                                                       'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,'Northern Ireland',cc_deflators)
defra_co2_reg['Scotland'] = defra.makeregionresults(S,U,Y,newY,y_regional['scotland'],meta,co2,defra_ghg_uk['direct'],
                                                       'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,'Scotland',cc_deflators)

file = os.path.join(census_filepath, 'laregionlookup2022.xlsx')   
LAcodesnames = pd.read_excel(file, header = 4, index_col=0)
oac01_names = pd.read_excel(os.path.join(census_filepath, 'oac2001.xlsx'), sheet_name=None, index_col = 0)
oac11_names = pd.read_excel(os.path.join(census_filepath, 'oac2011.xlsx'), sheet_name=None, index_col = 0)

ghg_reg_las = {}; oac_lookup = {}
for reg in reg_england + ['Scotland', 'Wales']:
    ghg_reg_las[reg], oac_lookup[reg] = defra.makelaresults(defra_ghg_reg[reg],reglaspropyr,regpophholdsyr,LAcodesnames,regoacsyr,oacyrmeta,oacyrspends,oac01_names,oac11_names,reg,years,concs_dict,cc_deflators)
    
for reg in reg_england:
    defra.printdefradata(reg,results_filepath,years,[defra_ghg_reg[reg]],['ghg'])

for reg in list(ghg_reg_las.keys()):
    for la in ghg_reg_las[reg]:
        defra.printdefradata2(la,results_filepath,years,ghg_reg_las[reg][la],'ghg')
        
co2_reg_las = {}; oac_lookup = {}
for reg in reg_england + ['Scotland', 'Wales']:
    co2_reg_las[reg], oac_lookup[reg] = defra.makelaresults(defra_co2_reg[reg],reglaspropyr,regpophholdsyr,LAcodesnames,regoacsyr,oacyrmeta,oacyrspends,oac01_names,oac11_names,reg,years,concs_dict,cc_deflators)
    
for reg in reg_england:
    defra.printdefradata(reg,results_filepath,years,[defra_co2_reg[reg]],['co2'])

for reg in list(ghg_reg_las.keys()):
    for la in co2_reg_las[reg]:
        defra.printdefradata2(la,results_filepath,years,co2_reg_las[reg][la],'co2')

# formula = 'defra.regioncheck2023('
# for reg in reg_england + ['Scotland', 'Wales', 'Northern Ireland']:
#     formula += 'defra_ghg_reg["' + reg + '"], '
# formula += 'years)'

# region_results_check = eval(formula)

region_results_check = defra.regioncheck(defra_ghg_reg,years)
region_results_check2 = defra.regioncheck(defra_co2_reg,years)

defrawriter = os.path.join(results_filepath, 'Region Mastersheet.xlsx')
writer = pd.ExcelWriter(defrawriter)
for item in region_results_check:
    print(item)
    sheetlabel = item
    region_results_check[item].to_excel(writer,str(sheetlabel))    
writer.save()

la_results_check = {}
for reg in reg_england + ['Scotland', 'Wales']:
    la_results_check[reg] = defra.lacheck(defra_ghg_reg[reg],ghg_reg_las[reg],years)
    
la_results_check2 = {}
for reg in reg_england + ['Scotland', 'Wales']:
    la_results_check2[reg] = defra.lacheck(defra_co2_reg[reg],co2_reg_las[reg],years)
    
    
for reg in reg_england + ['Scotland', 'Wales']:
    defrawriter = os.path.join(results_filepath, reg + ' Mastersheet.xlsx')
    writer = pd.ExcelWriter(defrawriter)
    for item in la_results_check[reg]:
        print(item)
        sheetlabel = item
        la_results_check[reg][item].to_excel(writer,str(sheetlabel))
    writer.save()

'''
# Below not yet working - need to fix
la_results_check = {}
for reg in reg_england + ['Scotland', 'Wales']:
    la_results_check[reg] = defra.lacheck(defra_ghg_reg[reg],ghg_reg_las[reg],years)


   


for reg in reg_england + ['Scotland', 'Wales']:
    defrawriter = os.path.join(results_filepath, reg + ' Mastersheet.xlsx')
    writer = pd.ExcelWriter(defrawriter)
    for item in la_results_check[reg]:
        print(item)
        sheetlabel = item
        la_results_check[reg][item].to_excel(writer,str(sheetlabel))
    writer.save()


footdata = lcf.processdataforfoots(defra_ghg_uk,hhspenddata3,years)

hh_age_foots = lcf.hhagefoots(footdata,years)
pop_age_foots = lcf.popagefoots(footdata,years)
hh_income_foots = lcf.hhincfoots(footdata,years)
pop_income_foots = lcf.popincfoots(footdata,years)
hh_income2_foots = lcf.hhinc2foots(footdata,years)
pop_income2_foots = lcf.popinc2foots(footdata,years)
hh_reg_foots = lcf.hhregfoots(footdata,years)
pop_reg_foots = lcf.popregfoots(footdata,years)
hh_oac_foots = lcf.hhoacfoots(footdata,years)
pop_oac_foots = lcf.popoacfoots(footdata,years)
'''