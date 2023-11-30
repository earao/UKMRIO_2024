#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 2021

CO2 emissions for MSOAs or LSOAs combining 2 years at a time, IO part adapted from code by Anne Owen

@author: lenakilian
"""

import pandas as pd
import pickle
import numpy as np

df = pd.DataFrame

################
# IO functions #
################

def make_Z_from_S_U(S,U):
    
    Z = np.zeros(shape = (np.size(S,0)+np.size(U,0),np.size(S,1)+np.size(U,1)))
    
    Z[np.size(S,0):,0:np.size(U,1)] = U
    Z[0:np.size(S,0),np.size(U,1):] = S
        
    return Z


def make_x(Z,Y):
    
    x = np.sum(Z,1)+np.sum(Y,1)
    x[x == 0] = 0.000000001
    
    return x


def make_L(Z,x):
    
    bigX = np.zeros(shape = (len(Z)))    
    bigX = np.tile(np.transpose(x),(len(Z),1))
    A = np.divide(Z,bigX)    
    L = np.linalg.inv(np.identity(len(Z))-A)

    return L

####################
# demand functions #
####################

def make_Yhh_106(Y_d,years,meta):
    
    total_Yhh_106 = {}
    col =  Y_d[years[0]].columns[0:36]
    idx =  Y_d[years[0]].index[0:106]
    for yr in years:
        temp = np.zeros(shape = [106,36])
        
        for r in range(0,meta['reg']['len']):
            temp  = temp + Y_d[yr].iloc[r*106:(r+1)*106,0:36].values
            
        total_Yhh_106[yr] = df(temp, index =idx, columns =col)
    
    return total_Yhh_106


def make_Yhh_109_34(Y_d,years,meta):
    
    total_Yhh_109 = {}
    col =  Y_d[years[0]].columns[0:34]
    idx =  Y_d[years[0]].index[0:109]
    for yr in years:
        temp = np.zeros(shape = [109,34])
        
        for r in range(0,meta['reg']['len']):
            temp  = temp + Y_d[yr].iloc[r*109:(r+1)*109,0:34].values
            
        total_Yhh_109[yr] = df(temp, index =idx, columns =col)
    
    return total_Yhh_109

##################
# LCFS functions #
##################

def convert43to41(Y,concs_dict,years):
    Y2 = {}
    for yr in years:
      temp = np.dot(Y[yr],concs_dict['C43_to_C41'])
      Y2[yr] = df(temp, index = Y[yr].index, columns = concs_dict['C43_to_C41'].columns)
    
    return Y2

def make_totals_2023(hhdspend,years):
  
  coicop_exp_tot = {}
  
  for yr in years:
      coicop_exp_tot[yr] = np.sum(hhdspend[yr],0)
  return coicop_exp_tot

def convert_exp_tot_sizes(coicop_exp_tot,concs_dict,years,size_str):

  coicop_exp_tot2 = {}
  for yr in years:
    temp = np.sum(np.dot(np.diag(coicop_exp_tot[yr]),concs_dict[size_str]),0)
    coicop_exp_tot2[yr] = df(temp, index = concs_dict[size_str].columns)
    
  return coicop_exp_tot2

def make_balanced_totals_2023(coicop_exp_tot2,total_Yhh_112,concs_dict,years):
  
  coicop_exp_tot3 = {}
  
  for yr in years:
    corrector = np.zeros(shape = 105)
    countstart = 0
    countend = 0
    for numb in range(0,34):
      conc = concs_dict[str(numb)+'a']
      countend = np.sum(np.sum(conc))+countstart
      lcf_subtotal = np.sum(np.dot(np.transpose(coicop_exp_tot2[yr]),conc)) #*52/1000)
      required_subtotal = np.sum(total_Yhh_112[yr].iloc[:,numb])
      correction_factor = required_subtotal/lcf_subtotal
      for c in range(countstart,countend):
        corrector[c] = correction_factor
      countstart = countend
    coicop_exp_tot3[yr] = np.dot(np.transpose(coicop_exp_tot2[yr]),np.diag(corrector))
  
  return coicop_exp_tot3

def make_y_hh_105(Y,coicop_exp_tot3,years,concs_dict,meta):
  
  yhh_wide = {}
  
  for yr in years:
    temp = np.zeros(shape = [meta['fd']['len_idx'], 105])

    countstart = 0
    countend = 0
    col = []
    for a in range(0,34):
      conc = np.tile(concs_dict[str(a)],(meta['reg']['len'],1))
      countend = np.sum(np.sum(concs_dict[str(a)+'a']))+countstart
      category_total = np.dot(coicop_exp_tot3[yr],concs_dict[str(a)+'a'])
      test1 = np.dot(conc,np.diagflat(category_total))
      test2 = np.tile(np.dot(conc,np.transpose(category_total)),(1,np.size(conc,1)))
      test3 = test1/test2
      test3 = np.nan_to_num(test3, copy=True)
      test4 = np.dot(np.diag(Y[yr].iloc[:,a]),test3)
      temp[:,countstart:countend] = test4
      col[countstart:countend] = concs_dict[str(a) + 'a'].columns
      countstart = countend
    yhh_wide[yr] = df(temp, index = Y[yr].index, columns = col)
      
  return yhh_wide


def make_new_Y_105(Y,yhh_wide,years):
  newY = {}
  col = []
  
  for yr in years:
    temp = np.zeros(shape=[len(Y[yr]),112])
    temp[:,0:105] = yhh_wide[yr]
    temp[:,105:112] = Y[yr].iloc[:,34:41]
    col[0:105] = yhh_wide[yr].columns
    col[105:112] = Y[yr].iloc[:,34:41].columns
    newY[yr] = df(temp, index = Y[yr].index, columns = col)
  
  return newY


def make_y_hh_prop(Y,total_Yhh_106,meta,years):
    
    yhh_prop = {}
    
    for yr in years:
        temp = np.zeros(shape=(len(Y[yr])))
    
        for r in range(0,meta['reg']['len']):
            temp[r*106:(r+1)*106] = np.divide(np.sum(Y[yr].iloc[r*106:(r+1)*106,0:36],1),np.sum(total_Yhh_106[yr],1))
            np.nan_to_num(temp, copy = False)
        
        yhh_prop[yr] = temp

        
    return yhh_prop


def make_new_Y(Y,yhh_wide,meta,years):
    
    newY = {}
    col = []
    
    for yr in years:
        temp = np.zeros(shape=[len(Y[yr]),314])
        temp[:,0:307] = yhh_wide[yr]
        temp[:,307:314] = Y[yr].iloc[:,33:40]
        col[0:307] = yhh_wide[yr].columns
        col[307:314] = Y[yr].iloc[:,33:40].columns
        newY[yr] = df(temp, index = Y[yr].index, columns = col)
    
    return newY

def make_ylcf_props(hhdspend,years, concs_dict, size_str):
    
    ylcf_props = {}; ylcfs_total = {}
    
    var_dict = df(concs_dict[size_str].stack())
    var_dict = var_dict.loc[var_dict[0] > 0].reset_index()
    var_dict = dict(zip(var_dict['level_0'], var_dict['level_1']))
    
    for yr in years:
        temp = hhdspend[yr].rename(columns=var_dict).sum(axis=1, level=0)
        ylcf_props[yr] = df(temp.apply(lambda x: x/x.sum()))
        ylcfs_total[yr] = temp.sum()
        
    return ylcf_props, ylcfs_total


def makefoot(S,U,Y,stressor,years):
    footbyCOICOP = {}
    for yr in years:
        temp = np.zeros(shape = np.size(Y[yr],1))
        Z = make_Z_from_S_U(S[yr], U[yr]) 
        bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
        bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:] = Y[yr]     
        x = make_x(Z,bigY)
        L = make_L(Z,x)
        bigstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
        bigstressor[0:np.size(Y[yr],0),:] = stressor[yr]
        e = np.sum(bigstressor,1)/x
        eL = np.dot(e,L)
        for a in range(np.size(Y[yr],1)):
            temp[a] = np.dot(eL,bigY[:,a])
        footbyCOICOP[yr] = temp  

    return footbyCOICOP

###########
# Run all #
###########

def make_footprint(hhdspend, wd, concs_filepath):
    
    """
    Calculate consumption-based household GHG emissions for MSOAs or LSOAs from the LCFS (emissios calculated in LCFS_aggregation_combined_years.py) and the UKMRIO 2020
    """

#############
# load data #
#############

    # load meta data from [UKMRIO]
    meta = pickle.load(open(wd + 'meta.p', "rb" ))
   
    # create year lists
    years = list(hhdspend.keys())

    # load and clean up concs to make it usable
    # these translate IO data sectors to LCFS products/services
    concs_dict = pd.read_excel(concs_filepath + 'ONS_to_COICOP_LCF_concs.xlsx', sheet_name=None, index_col=0)

#######################
# aggregate emissions #
#######################

    # get mean from 2 years
    # calculate differnece between years in household data to calculate means for other vairables
    
    # Load UKMRIO and calculate means for UKMRIO data
    ukmrio = {}; #means = {}
    for data in ['ghg', 'uk_ghg_direct', 'S', 'U', 'Y']:
        ukmrio[data] = pickle.load(open(wd + data + '.p', "rb" ))
        
    ukmrio['Y'] = convert43to41(ukmrio['Y'], concs_dict, years)
    total_Yhh_109 = make_Yhh_109_34(ukmrio['Y'], years, meta)
    
    coicop_exp_tot = make_totals_2023(hhdspend, years)
    coicop_exp_tot2 = convert_exp_tot_sizes(coicop_exp_tot, concs_dict, years, '456_to_105')
    coicop_exp_tot3 = make_balanced_totals_2023(coicop_exp_tot2, total_Yhh_109, concs_dict, years) # cannot multiply because it is 109 instead of 112 columns

    yhh_wide = make_y_hh_105(ukmrio['Y'], coicop_exp_tot3, years, concs_dict, meta)
    
    newY = make_new_Y_105(ukmrio['Y'], yhh_wide, years)

    ylcf_props, ylcfs_total = make_ylcf_props(hhdspend, list(hhdspend.keys()), concs_dict, '456_to_105')

    COICOP_ghg = makefoot(ukmrio['S'], ukmrio['U'], newY, ukmrio['ghg'], list(hhdspend.keys()))
    
    Total_ghg = {}; multipliers = {}
    for year in list(hhdspend.keys()):
        # add index
        new_index = concs_dict['456_to_105'].columns.tolist() + ukmrio['Y'][year].loc[:,'13 Non-profit instns serving households':].columns.tolist()
        COICOP_ghg[year] = df(COICOP_ghg[year], index=new_index, columns=['total_ghg'])
        # add direct emissions
        for item in new_index:
            if '4.5.2' in item:
                gas_direct = item
            if '7.2.2' in item:
                travel_direct = item
        COICOP_ghg[year].loc[gas_direct, 'total_ghg'] += ukmrio['uk_ghg_direct'][year]['Consumer expenditure - not travel']
        COICOP_ghg[year].loc[travel_direct, 'total_ghg'] += ukmrio['uk_ghg_direct'][year]['Consumer expenditure - travel']
        
        # multipliers tCO2e/GBP 
        multipliers[year] = COICOP_ghg[year].join(df(ylcfs_total[year], columns=['total_spend']), how='right')
        multipliers[year]['multipliers'] = multipliers[year]['total_ghg'] / multipliers[year]['total_spend']
    
        # this gives GHG emissions for the groups, break down to per capita emissions
        temp = ylcf_props[year].T.apply(lambda x: x*COICOP_ghg[year]['total_ghg'])
        Total_ghg[year] = temp.T[ylcf_props[year].columns.tolist()]
        
        print('Emissions calculated for ' + str(year))
    
    return(Total_ghg, multipliers, yhh_wide)