# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:32:01 2023

This code makes the spreadsheets for all footprint types for the UK, England, Scotland, Wales and NI

@authors: Anne Owen & Lena Kilian  
"""

import pandas as pd
import numpy as np
import os
df = pd.DataFrame
import LCF_functions as lcf
import pickle
import defra_functions as defra
import ukmrio_functions as uk
import census_functions as cs
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
lcf_filepath = wd +'data/raw/LCFS/'
deflator_filepath = wd +'data/raw/ONS/ONS deflators/'
census_filepath = wd +'data/raw/census data/'

years = np.array([int(x) for x in range(2001,2021)])
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
watblu_cons = pickle.load( open(results_filepath + "watblu_cons.p", "rb" ) )
uk_watbluc_direct = pickle.load( open(results_filepath + "uk_wat_blu_cons_direct.p", "rb" ) )
watgrn_cons = pickle.load( open(results_filepath + "watgrn_cons.p", "rb" ) )
watblu_wdrl = pickle.load( open(results_filepath + "watblu_wdrl.p", "rb" ) )
uk_watbluw_direct = pickle.load( open(results_filepath + "uk_wat_blu_wdrl_direct.p", "rb" ) )
mat = pickle.load( open(results_filepath + "mat.p", "rb" ) )
bio = pickle.load( open(results_filepath + "bio.p", "rb" ) )
ore = pickle.load( open(results_filepath + "ore.p", "rb" ) )
nmm = pickle.load( open(results_filepath + "nmm.p", "rb" ) )
ffl = pickle.load( open(results_filepath + "ffl.p", "rb" ) )
nrg = pickle.load( open(results_filepath + "nrg.p", "rb" ) )
uk_nrg_direct = pickle.load( open(results_filepath + "uk_nrg_direct.p", "rb" ) )

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

hhspenddata3 = lcf.convert_hhspend_sizes(hhspenddata2,concs_dict,years,'456_to_105')

regions = ['North East','North West','Yorkshire and The Humber','East Midlands','West Midlands','East','London','South East','South West','Scotland','Wales']
regions_lc = ['north_east','north_west','yorkshire','east_midlands','west_midlands','east','london','south_east','south_west','scotland','wales']
region_dict = dict(zip(regions_lc, list(range(1,13))))

oa_lookup01 = cs.make_01_lookup(census_filepath) 
(oa_lookup11,s11) = cs.make_11_lookup(census_filepath) 
oa_lookup21 = cs.make_21_lookup(census_filepath,s11)

regoacsyr = pickle.load(open(results_filepath + "regoacsyr.p", "rb" ))

# regoacsyr = lcf.make_pop_hhold_by_oac_region_year_2023(oa_lookup01,oa_lookup11,census_filepath,regions)

hhspenddata3 = lcf.makestrings(hhspenddata3,oac_2001_years,oac_2011_years)
(oacyrspends,oacyrmeta) = lcf.makeoacspends_2023(hhspenddata3,oac_none_years,oac_2001_years,oac_2011_years)
oacyrspends = lcf.correct_index(oacyrspends,regions)

(reglaspendyr,regpophholdsyr) = lcf.make_la_spends_pop_by_region_year_2023(regoacsyr,oacyrspends,regions,years) 
reglaspendyr = lcf.correct_reglaspendyr_zero(reglaspendyr,regions,years)                 

regpophholdsyr = {}
for yr in years:
    reglaspend = {}
    regpophholds = {}
    for r in regions:
        ladpophholds = regoacsyr[yr][r].groupby(['LAD_code']).sum()
        regpophholds[r] = ladpophholds
    regpophholdsyr[yr] = regpophholds

y_regions = lcf.make_y_regions_2023(wd, hhspenddata3, regions, regpophholdsyr, newY, years) # need to define regions and regpophholdsyr

(S_d, U_d, Y_d, newY_d, y_regions_d) = uk.deflate_io_regions(S, U, Y, newY, y_regions, allyears, years, io_deflators, meta)

v = uk.make_v(U_d, Y_d, allyears, meta)

emptydirect = df(np.zeros((1, len(allyears))), index = ['direct'], columns = allyears)

uk_watbluc_direct = np.transpose(uk_watbluc_direct)
uk_watbluw_direct = np.transpose(uk_watbluw_direct)

temp = np.zeros((1, 31))
temp[:, 0] = uk_watbluc_direct.iloc[:, 0]
temp[:, 1] = uk_watbluc_direct.iloc[:, 0]
temp[:, 2] = uk_watbluc_direct.iloc[:, 0]
temp[:, 3] = uk_watbluc_direct.iloc[:, 0]
temp[:, 4] = uk_watbluc_direct.iloc[:, 0]
temp[:, 5:31] = uk_watbluc_direct.iloc[:, 0:26]
uk_watbluc_direct2 = df(temp, columns = allyears)

temp = np.zeros((1, 31))
temp[:, 0] = uk_watbluw_direct.iloc[:, 0]
temp[:, 1] = uk_watbluw_direct.iloc[:, 0]
temp[:, 2] = uk_watbluw_direct.iloc[:, 0]
temp[:, 3] = uk_watbluw_direct.iloc[:, 0]
temp[:, 4] = uk_watbluw_direct.iloc[:, 0]
temp[:, 5:31] = uk_watbluw_direct.iloc[:, 0:26]
uk_watbluw_direct2 = df(temp, columns = allyears)

y_regional = {}
for region in regions_lc:
    y_regional[region] = lcf.make_y_regional(region, region_dict, y_regions, years)
        
(y_England,y_N_Ireland,y_Scotland,y_Wales) = lcf.make_y_countries_2023(y_regions,years)
(y_England_d,y_N_Ireland_d,y_Scotland_d,y_Wales_d) = lcf.make_y_countries_2023(y_regions_d,years)

defra_ghg_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, ghg, uk_ghg_direct, 'ghg', allyears, years)
defra_co2_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, co2, uk_co2_direct, 'co2', allyears, years)
defra_nrg_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, nrg, uk_nrg_direct, 'nrg', allyears, years)
defra_mat_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, mat, emptydirect, 'mat', allyears, years)
defra_bio_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, bio, emptydirect, 'bio', allyears, years)
defra_ore_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, ore, emptydirect, 'ore', allyears, years)
defra_ffl_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, ffl, emptydirect, 'ffl', allyears, years)
defra_nmm_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, nmm, emptydirect, 'nmm', allyears, years)
defra_blc_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, watblu_cons, uk_watbluc_direct2, 'blc', allyears, years)
defra_blw_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, watblu_wdrl, uk_watbluw_direct2, 'blw', allyears, years)
defra_grc_uk = defra.makeukresults2023(S, U, Y, newY, coicop_exp_tot2, meta, watgrn_cons, emptydirect, 'grc', allyears, years)
defra_gva_uk = defra.makeukresults2023(S_d, U_d, Y_d, newY_d, coicop_exp_tot2, meta, v, emptydirect, 'gva', allyears, years)

defra.printdefradata('UK', results_filepath, years, 
                     [defra_ghg_uk, defra_co2_uk, defra_nrg_uk, defra_mat_uk, defra_bio_uk, defra_ore_uk, defra_ffl_uk, defra_nmm_uk, 
                      defra_blc_uk, defra_blw_uk, defra_grc_uk, defra_gva_uk], 
                     ['ghg', 'co2', 'nrg', 'mat', 'bio', 'ore', 'ffl', 'nmm', 'blc', 'blw', 'grc', 'gva'])

defra_ghg_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,ghg,defra_ghg_uk['direct'],'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'England')
defra_co2_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,co2,defra_co2_uk['direct'],'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,'England')
defra_nrg_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,nrg,defra_nrg_uk['direct'],'nrg',years,concs_dict,defra_nrg_uk['coicop_mult'],defra_nrg_uk['sic_mult'],regpophholdsyr,'England')
defra_mat_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,mat,defra_mat_uk['direct'],'mat',years,concs_dict,defra_mat_uk['coicop_mult'],defra_mat_uk['sic_mult'],regpophholdsyr,'England')
defra_bio_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,bio,defra_bio_uk['direct'],'bio',years,concs_dict,defra_bio_uk['coicop_mult'],defra_bio_uk['sic_mult'],regpophholdsyr,'England')
defra_ore_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,ore,defra_ore_uk['direct'],'ore',years,concs_dict,defra_ore_uk['coicop_mult'],defra_ore_uk['sic_mult'],regpophholdsyr,'England')
defra_ffl_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,ffl,defra_ffl_uk['direct'],'ffl',years,concs_dict,defra_ffl_uk['coicop_mult'],defra_ffl_uk['sic_mult'],regpophholdsyr,'England')
defra_nmm_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,nmm,defra_nmm_uk['direct'],'nmm',years,concs_dict,defra_nmm_uk['coicop_mult'],defra_nmm_uk['sic_mult'],regpophholdsyr,'England')
defra_blc_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,watblu_cons,defra_blc_uk['direct'],'blc',years,concs_dict,defra_blc_uk['coicop_mult'],defra_blc_uk['sic_mult'],regpophholdsyr,'England')
defra_blw_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,watblu_wdrl,defra_blw_uk['direct'],'blw',years,concs_dict,defra_blw_uk['coicop_mult'],defra_blw_uk['sic_mult'],regpophholdsyr,'England')
defra_grc_eng = defra.makeregionresults2023(S,U,Y,newY,y_England,meta,watgrn_cons,defra_grc_uk['direct'],'grc',years,concs_dict,defra_grc_uk['coicop_mult'],defra_grc_uk['sic_mult'],regpophholdsyr,'England')
defra_gva_eng = defra.makeregionresults2023(S_d,U_d,Y_d,newY_d,y_England_d,meta,v,defra_gva_uk['direct'],'gva',years,concs_dict,defra_gva_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'England')

defra.printdefradata('England', results_filepath, years, 
                     [defra_ghg_eng, defra_co2_eng, defra_nrg_eng, defra_mat_eng, defra_bio_eng, defra_ore_eng, defra_ffl_eng, defra_nmm_eng, 
                      defra_blc_eng, defra_blw_eng, defra_grc_eng, defra_gva_eng], 
                     ['ghg', 'co2', 'nrg', 'mat', 'bio', 'ore', 'ffl', 'nmm', 'blc', 'blw', 'grc', 'gva'])

defra_ghg_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,ghg,defra_ghg_uk['direct'],'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_co2_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,co2,defra_co2_uk['direct'],'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_nrg_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,nrg,defra_nrg_uk['direct'],'nrg',years,concs_dict,defra_nrg_uk['coicop_mult'],defra_nrg_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_mat_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,mat,defra_mat_uk['direct'],'mat',years,concs_dict,defra_mat_uk['coicop_mult'],defra_mat_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_bio_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,bio,defra_bio_uk['direct'],'bio',years,concs_dict,defra_bio_uk['coicop_mult'],defra_bio_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_ore_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,ore,defra_ore_uk['direct'],'ore',years,concs_dict,defra_ore_uk['coicop_mult'],defra_ore_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_ffl_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,ffl,defra_ffl_uk['direct'],'ffl',years,concs_dict,defra_ffl_uk['coicop_mult'],defra_ffl_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_nmm_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,nmm,defra_nmm_uk['direct'],'nmm',years,concs_dict,defra_nmm_uk['coicop_mult'],defra_nmm_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_blc_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,watblu_cons,defra_blc_uk['direct'],'blc',years,concs_dict,defra_blc_uk['coicop_mult'],defra_blc_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_blw_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,watblu_wdrl,defra_blw_uk['direct'],'blw',years,concs_dict,defra_blw_uk['coicop_mult'],defra_blw_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_grc_sco = defra.makeregionresults2023(S,U,Y,newY,y_Scotland,meta,watgrn_cons,defra_grc_uk['direct'],'grc',years,concs_dict,defra_grc_uk['coicop_mult'],defra_grc_uk['sic_mult'],regpophholdsyr,'Scotland')
defra_gva_sco = defra.makeregionresults2023(S_d,U_d,Y_d,newY_d,y_Scotland_d,meta,v,defra_gva_uk['direct'],'gva',years,concs_dict,defra_gva_uk['coicop_mult'],defra_gva_uk['sic_mult'],regpophholdsyr,'Scotland')

defra.printdefradata('Scotland', results_filepath, years, 
                     [defra_ghg_sco, defra_co2_sco, defra_nrg_sco, defra_mat_sco, defra_bio_sco, defra_ore_sco, defra_ffl_sco, defra_nmm_sco, 
                      defra_blc_sco, defra_blw_sco, defra_grc_sco, defra_gva_sco], 
                     ['ghg', 'co2', 'nrg', 'mat', 'bio', 'ore', 'ffl', 'nmm', 'blc', 'blw', 'grc', 'gva'])

defra_ghg_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,ghg,defra_ghg_uk['direct'],'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'Wales')
defra_co2_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,co2,defra_co2_uk['direct'],'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,'Wales')
defra_nrg_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,nrg,defra_nrg_uk['direct'],'nrg',years,concs_dict,defra_nrg_uk['coicop_mult'],defra_nrg_uk['sic_mult'],regpophholdsyr,'Wales')
defra_mat_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,mat,defra_mat_uk['direct'],'mat',years,concs_dict,defra_mat_uk['coicop_mult'],defra_mat_uk['sic_mult'],regpophholdsyr,'Wales')
defra_bio_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,bio,defra_bio_uk['direct'],'bio',years,concs_dict,defra_bio_uk['coicop_mult'],defra_bio_uk['sic_mult'],regpophholdsyr,'Wales')
defra_ore_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,ore,defra_ore_uk['direct'],'ore',years,concs_dict,defra_ore_uk['coicop_mult'],defra_ore_uk['sic_mult'],regpophholdsyr,'Wales')
defra_ffl_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,ffl,defra_ffl_uk['direct'],'ffl',years,concs_dict,defra_ffl_uk['coicop_mult'],defra_ffl_uk['sic_mult'],regpophholdsyr,'Wales')
defra_nmm_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,nmm,defra_nmm_uk['direct'],'nmm',years,concs_dict,defra_nmm_uk['coicop_mult'],defra_nmm_uk['sic_mult'],regpophholdsyr,'Wales')
defra_blc_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,watblu_cons,defra_blc_uk['direct'],'blc',years,concs_dict,defra_blc_uk['coicop_mult'],defra_blc_uk['sic_mult'],regpophholdsyr,'Wales')
defra_blw_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,watblu_wdrl,defra_blw_uk['direct'],'blw',years,concs_dict,defra_blw_uk['coicop_mult'],defra_blw_uk['sic_mult'],regpophholdsyr,'Wales')
defra_grc_wal = defra.makeregionresults2023(S,U,Y,newY,y_Wales,meta,watgrn_cons,defra_grc_uk['direct'],'grc',years,concs_dict,defra_grc_uk['coicop_mult'],defra_grc_uk['sic_mult'],regpophholdsyr,'Wales')
defra_gva_wal = defra.makeregionresults2023(S_d,U_d,Y_d,newY_d,y_Wales_d,meta,v,defra_gva_uk['direct'],'gva',years,concs_dict,defra_gva_uk['coicop_mult'],defra_gva_uk['sic_mult'],regpophholdsyr,'Wales')

defra.printdefradata('Wales', results_filepath, years, 
                     [defra_ghg_wal, defra_co2_wal, defra_nrg_wal, defra_mat_wal, defra_bio_wal, defra_ore_wal, defra_ffl_wal, defra_nmm_wal, 
                      defra_blc_wal, defra_blw_wal, defra_grc_wal, defra_gva_wal], 
                     ['ghg', 'co2', 'nrg', 'mat', 'bio', 'ore', 'ffl', 'nmm', 'blc', 'blw', 'grc', 'gva'])

defra_ghg_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,ghg,defra_ghg_uk['direct'],'ghg',years,concs_dict,defra_ghg_uk['coicop_mult'],defra_ghg_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_co2_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,co2,defra_co2_uk['direct'],'co2',years,concs_dict,defra_co2_uk['coicop_mult'],defra_co2_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_nrg_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,nrg,defra_nrg_uk['direct'],'nrg',years,concs_dict,defra_nrg_uk['coicop_mult'],defra_nrg_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_mat_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,mat,defra_mat_uk['direct'],'mat',years,concs_dict,defra_mat_uk['coicop_mult'],defra_mat_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_bio_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,bio,defra_bio_uk['direct'],'bio',years,concs_dict,defra_bio_uk['coicop_mult'],defra_bio_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_ore_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,ore,defra_ore_uk['direct'],'ore',years,concs_dict,defra_ore_uk['coicop_mult'],defra_ore_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_ffl_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,ffl,defra_ffl_uk['direct'],'ffl',years,concs_dict,defra_ffl_uk['coicop_mult'],defra_ffl_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_nmm_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,nmm,defra_nmm_uk['direct'],'nmm',years,concs_dict,defra_nmm_uk['coicop_mult'],defra_nmm_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_blc_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,watblu_cons,defra_blc_uk['direct'],'blc',years,concs_dict,defra_blc_uk['coicop_mult'],defra_blc_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_blw_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,watblu_wdrl,defra_blw_uk['direct'],'blw',years,concs_dict,defra_blw_uk['coicop_mult'],defra_blw_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_grc_ni = defra.makeregionresults2023(S,U,Y,newY,y_N_Ireland,meta,watgrn_cons,defra_grc_uk['direct'],'grc',years,concs_dict,defra_grc_uk['coicop_mult'],defra_grc_uk['sic_mult'],regpophholdsyr,'Northern Ireland')
defra_gva_ni = defra.makeregionresults2023(S_d,U_d,Y_d,newY_d,y_N_Ireland_d,meta,v,defra_gva_uk['direct'],'gva',years,concs_dict,defra_gva_uk['coicop_mult'],defra_gva_uk['sic_mult'],regpophholdsyr,'Northern Ireland')

defra.printdefradata('Northern Ireland', results_filepath, years, 
                     [defra_ghg_ni, defra_co2_ni, defra_nrg_ni, defra_mat_ni, defra_bio_ni, defra_ore_ni, defra_ffl_ni, defra_nmm_ni, 
                      defra_blc_ni, defra_blw_ni, defra_grc_ni, defra_gva_ni], 
                     ['ghg', 'co2', 'nrg', 'mat', 'bio', 'ore', 'ffl', 'nmm', 'blc', 'blw', 'grc', 'gva'])