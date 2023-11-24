#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 13:46:00 2018
This file explores the LCF original survey data from the UKDA
@author: earao
"""

import numpy as np
import pandas as pd
import os
df = pd.DataFrame
import copy as cp

###########################
# used in defra_main_2023 #
###########################

def make_Yhh_109_34(Y_d,years,meta):
    
    total_Yhh_109 = {}
    col = Y_d[2016].columns[0:34]
    idx = Y_d[2016].index[0:109]
    for yr in years:
        temp = np.zeros(shape = [109,34])
        
        for r in range(0,meta['reg']['len']):
            temp  = temp + Y_d[yr].iloc[r*109:(r+1)*109,0:34].values
            
        total_Yhh_109[yr] = df(temp, index =idx, columns =col)
    
    return total_Yhh_109



def convert_hhspend_sizes(hhspenddata,concs_dict,years,size_str):

    concs_dict[size_str] = concs_dict[size_str].apply(lambda x: pd.to_numeric(x, errors='coerce')).dropna(axis=0, how='all').dropna(axis=1, how='all')
  
    hhspenddata2 = {}
    for yr in years:
        loc = hhspenddata[yr].columns.get_loc('1.1.1.1.1')
        beginning = hhspenddata[yr].iloc[:,:loc]
        end = hhspenddata[yr].iloc[:,loc:]
        temp = np.array(end.apply(lambda x: pd.to_numeric(x, errors='coerce')).fillna(0))
        temp2 = np.dot(temp, concs_dict[size_str])
        end = df(temp2, columns = concs_dict[size_str].columns, index = hhspenddata[yr].index)
 
        hhspenddata2[yr] = beginning.join(end)
        
    return hhspenddata2
        
def convert_exp_tot_sizes(coicop_exp_tot,concs_dict,years,size_str):

  coicop_exp_tot2 = {}
  for yr in years:
    temp = np.sum(np.dot(np.diag(coicop_exp_tot[yr]),concs_dict[size_str]),0)
    coicop_exp_tot2[yr] = df(temp, index = concs_dict[size_str].columns)
    
  return coicop_exp_tot2

def convert43to41(Y,concs_dict,years):
  
  Y2 = {}
  
  for yr in years:
    temp = np.dot(Y[yr],concs_dict['C43_to_C41'])
    Y2[yr] = df(temp, index = Y[yr].index, columns = concs_dict['C43_to_C41'].columns)
  
  return Y2

def make_balanced_totals(coicop_exp_tot2,total_Yhh_112,concs_dict,years):
  
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

def make_totals(hhspenddata,years):
  
  coicop_exp_tot = {}
  
  for yr in years:
      hhspenddata[yr].columns
      coicop_exp_tot[yr] = np.sum(hhspenddata[yr].loc[:, '1.1.1.1.1':],0)
  return coicop_exp_tot

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

def make_y_regions_2023(wd, hhspenddata3, regions, regpophholdsyr, newY, years):
  y_regions={}
  ytemp = np.zeros((np.size(newY[2001],0),112))
  filepath = wd + 'data/raw/BEIS energy/'
  file = os.path.join(filepath, 'Sub-national_energy_consumption_statistics_2005-2020.xlsx')
  gas_prop = pd.read_excel(file, sheet_name = 'reg_prop_gas', index_col = 0)
  elc_prop = pd.read_excel(file, sheet_name = 'reg_prop_elec', index_col = 0)
  nipop = [1686000,1693000,1701000,1709000,1721000,1735000,1752000,1770000,1786000,1799000,1810000,1819000,1827000,1835000,1846000,1857000,1867000,1876000,1885000,1896000]
  for y, yr in enumerate(years):
#    totalhholds = np.sum(hhspenddata3[yr].iloc[:,0],0)
    totalhholdspend = np.sum(hhspenddata3[yr].loc[:,'1.1.1 Bread and cereals':],0)
    temp = hhspenddata3[yr].drop('OA class 1',1)
    temp = temp.drop('OA class 2',1)
    temp = temp.groupby(['GOR']).sum()
    spendbyregion = temp.loc[:,'1.1.1 Bread and cereals':]
#    regpop = temp.loc[:,'weight']
#    for reg in range(1,13):
    totalpop=0
    for r, reg in enumerate(regions):
      totalpop = totalpop + np.sum(regpophholdsyr[yr][reg]['pop'])
    totalpop = totalpop+nipop[y]
    regions_eng = regions
    for r, reg in enumerate(regions_eng):
      prop_spend = df(spendbyregion.loc[r+1,:]/totalhholdspend,index=spendbyregion.columns)
      prop_spend.iloc[30] = elc_prop.loc[r+1,yr]
      prop_spend.iloc[31] = gas_prop.loc[r+1,yr]
      
      prop_pop = np.sum(regpophholdsyr[yr][reg]['pop'])/totalpop
      ytemp[:,0:105] = np.multiply(newY[yr].iloc[:,0:105],np.tile(np.transpose(prop_spend),(np.size(newY[2001],0),1)))
      ytemp[:,105:112] = np.multiply(newY[yr].iloc[:,105:112],np.tile(prop_pop,(np.size(newY[2001],0),7)))
      y_regions[str(yr)+'_'+str(r+1)] = df(ytemp,index=newY[yr].index,columns=newY[yr].columns)
      y_regions[str(yr)+'_'+str(r+1)] = y_regions[str(yr)+'_'+str(r+1)].fillna(0)
    
    r=10
    
    prop_spend = df(spendbyregion.loc[r,:]/totalhholdspend,index=spendbyregion.columns)
    prop_spend.iloc[30] = elc_prop.loc[r,yr]
    prop_spend.iloc[31] = gas_prop.loc[r,yr]
      
    prop_pop = np.sum(regpophholdsyr[yr]['Wales']['pop'])/totalpop
    ytemp[:,0:105] = np.multiply(newY[yr].iloc[:,0:105],np.tile(np.transpose(prop_spend),(np.size(newY[2001],0),1)))
    ytemp[:,105:112] = np.multiply(newY[yr].iloc[:,105:112],np.tile(prop_pop,(np.size(newY[2001],0),7)))
    y_regions[str(yr)+'_'+str(r)] = df(ytemp,index=newY[yr].index,columns=newY[yr].columns)
    y_regions[str(yr)+'_'+str(r)] = y_regions[str(yr)+'_'+str(r)].fillna(0)
    
    r=11
    
    prop_spend = df(spendbyregion.loc[r,:]/totalhholdspend,index=spendbyregion.columns)
    prop_spend.iloc[30] = elc_prop.loc[r,yr]
    prop_spend.iloc[31] = gas_prop.loc[r,yr]
      
    prop_pop = np.sum(regpophholdsyr[yr]['Scotland']['pop'])/totalpop
    ytemp[:,0:105] = np.multiply(newY[yr].iloc[:,0:105],np.tile(np.transpose(prop_spend),(np.size(newY[2001],0),1)))
    ytemp[:,105:112] = np.multiply(newY[yr].iloc[:,105:112],np.tile(prop_pop,(np.size(newY[2001],0),7)))
    y_regions[str(yr)+'_'+str(r)] = df(ytemp,index=newY[yr].index,columns=newY[yr].columns)
    y_regions[str(yr)+'_'+str(r)] = y_regions[str(yr)+'_'+str(r)].fillna(0)
    
    r=12
    
    prop_spend = df(spendbyregion.loc[r,:]/totalhholdspend,index=spendbyregion.columns)
    prop_spend.iloc[30] = elc_prop.loc[r,yr]
    prop_spend.iloc[31] = gas_prop.loc[r,yr]
      
    prop_pop = nipop[y]/totalpop
    ytemp[:,0:105] = np.multiply(newY[yr].iloc[:,0:105],np.tile(np.transpose(prop_spend),(np.size(newY[2001],0),1)))
    ytemp[:,105:112] = np.multiply(newY[yr].iloc[:,105:112],np.tile(prop_pop,(np.size(newY[2001],0),7)))
    y_regions[str(yr)+'_'+str(r)] = df(ytemp,index=newY[yr].index,columns=newY[yr].columns)
    y_regions[str(yr)+'_'+str(r)] = y_regions[str(yr)+'_'+str(r)].fillna(0)
    
  return y_regions

##########################################
# used in defra_uk_devolved_regions_2023 #
##########################################

# Also used in defra_main_2023 (see defra_main_2023):
# convert_hhspend_sizes
# convert_exp_tot_sizes
# convert43to41
# make_balanced_totals
# make_new_Y_109
# make_totals
# make_y_hh_109
# make_y_regions_2023

def make_y_countries_2023(y_regions,years):
  y_England={}
  y_N_Ireland={}
  y_Scotland={}
  y_Wales={}
  for yr in years:
    temp = np.zeros((np.size(y_regions['2001_1'],0),112))
    for reg in range(1,10):
      temp = temp+y_regions[str(yr)+'_'+str(reg)].values
    y_England[yr] = df(temp,index=y_regions[str(yr)+'_1'].index,columns=y_regions[str(yr)+'_1'].columns)
    y_N_Ireland[yr] = df(y_regions[str(yr)+'_12'].values,index=y_regions[str(yr)+'_12'].index,columns=y_regions[str(yr)+'_12'].columns)
    y_Scotland[yr] = df(y_regions[str(yr)+'_11'].values,index=y_regions[str(yr)+'_11'].index,columns=y_regions[str(yr)+'_11'].columns)
    y_Wales[yr] = df(y_regions[str(yr)+'_10'].values,index=y_regions[str(yr)+'_10'].index,columns=y_regions[str(yr)+'_10'].columns)
  
  return (y_England,y_N_Ireland,y_Scotland,y_Wales)

def removeoutliers(hhspenddata2,years):
# replaces spends that are larger than 4 Standard Deviations from mean with a value exactly 4 standard deviations from mean  
  for yr in years:
    
    for i in range(18,np.size(hhspenddata2[yr],1)):
      stdev4 = np.std(hhspenddata2[yr].iloc[:,i][hhspenddata2[yr].iloc[:,i]!=0])*4
      hhspenddata2[yr].iloc[:,i] = np.where(hhspenddata2[yr].iloc[:,i] > stdev4,stdev4,hhspenddata2[yr].iloc[:,i])
    
  return hhspenddata2


#####################################
# used in scotlanddeciles_2023_main #
#####################################

# Also used in defra_main_2023 (see defra_main_2023):
# convert43to41
# convert_hhspend_sizes
# convert_exp_tot_sizes
# make_balanced_totals
# make_new_Y_109
# make_totals
# make_y_hh_109
# removeoutliers

def processdataforscotdecile(ghg_mults,hhspenddata3,years):  
  footdata={}
  
  agg_filepath = 'C:/Users/earao/OneDrive - University of Leeds/Projects/milena/'
  agg = pd.read_excel(os.path.join(agg_filepath, 'aggregation.xlsx'), sheet_name='Sheet1', header = 0, index_col=0)
 
  
  for yr in years:
    temp = hhspenddata3[yr].loc[hhspenddata3[yr]['GOR']==11]
    
    temp = temp.reset_index()
      
    
    temp = temp.sort_values(by = "Income anonymised")
    
    temp = temp.reset_index()
    
    for i in range(0,len(temp)):
      if i == 0:
        temp.loc[i,'cum rank'] = temp.loc[i,'weight']
      else:
        temp.loc[i,'cum rank'] = temp.loc[i-1,'cum rank']+temp.loc[i,'weight']
        
    decile = np.sum(temp['weight'])/20
    
    for i in range(0,len(temp)):
      if(temp.loc[i,'cum rank'] ) < decile:
        temp.loc[i,'income decile'] = 1
      elif(temp.loc[i,'cum rank'] ) < decile*2:
        temp.loc[i,'income decile'] = 2
      elif(temp.loc[i,'cum rank'] ) < decile*3:
        temp.loc[i,'income decile'] = 3
      elif(temp.loc[i,'cum rank'] ) < decile*4:
        temp.loc[i,'income decile'] = 4
      elif(temp.loc[i,'cum rank'] ) < decile*5:
        temp.loc[i,'income decile'] = 5
      elif(temp.loc[i,'cum rank'] ) < decile*6:
        temp.loc[i,'income decile'] = 6
      elif(temp.loc[i,'cum rank'] ) < decile*7:
        temp.loc[i,'income decile'] = 7
      elif(temp.loc[i,'cum rank'] ) < decile*8:
        temp.loc[i,'income decile'] = 8
      elif(temp.loc[i,'cum rank'] ) < decile*9:
        temp.loc[i,'income decile'] = 9
      elif(temp.loc[i,'cum rank'] ) < decile*10:
        temp.loc[i,'income decile'] = 10
      elif(temp.loc[i,'cum rank'] ) < decile*11:
        temp.loc[i,'income decile'] = 11
      elif(temp.loc[i,'cum rank'] ) < decile*12:
        temp.loc[i,'income decile'] = 12
      elif(temp.loc[i,'cum rank'] ) < decile*13:
        temp.loc[i,'income decile'] = 13
      elif(temp.loc[i,'cum rank'] ) < decile*14:
        temp.loc[i,'income decile'] = 14
      elif(temp.loc[i,'cum rank'] ) < decile*15:
        temp.loc[i,'income decile'] = 15
      elif(temp.loc[i,'cum rank'] ) < decile*16:
        temp.loc[i,'income decile'] = 16
      elif(temp.loc[i,'cum rank'] ) < decile*17:
        temp.loc[i,'income decile'] = 17
      elif(temp.loc[i,'cum rank'] ) < decile*18:
        temp.loc[i,'income decile'] = 18
      elif(temp.loc[i,'cum rank'] ) < decile*19:
        temp.loc[i,'income decile'] = 19
      else:
        temp.loc[i,'income decile'] = 20
        
    temp['income equiv'] = temp['Income anonymised']/temp['OECD scale']*2
    
    temp = temp.sort_values(by = "income equiv")
    
    temp = temp.reset_index()
    
    for i in range(0,len(temp)):
      if i == 0:
        temp.loc[i,'cum rank2'] = temp.loc[i,'weight']
      else:
        temp.loc[i,'cum rank2'] = temp.loc[i-1,'cum rank2']+temp.loc[i,'weight']
        
    decile = np.sum(temp['weight'])/20
    
    for i in range(0,len(temp)):
      if(temp.loc[i,'cum rank2'] ) < decile:
        temp.loc[i,'income decile2'] = 1
      elif(temp.loc[i,'cum rank2'] ) < decile*2:
        temp.loc[i,'income decile2'] = 2
      elif(temp.loc[i,'cum rank2'] ) < decile*3:
        temp.loc[i,'income decile2'] = 3
      elif(temp.loc[i,'cum rank2'] ) < decile*4:
        temp.loc[i,'income decile2'] = 4
      elif(temp.loc[i,'cum rank2'] ) < decile*5:
        temp.loc[i,'income decile2'] = 5
      elif(temp.loc[i,'cum rank2'] ) < decile*6:
        temp.loc[i,'income decile2'] = 6
      elif(temp.loc[i,'cum rank2'] ) < decile*7:
        temp.loc[i,'income decile2'] = 7
      elif(temp.loc[i,'cum rank2'] ) < decile*8:
        temp.loc[i,'income decile2'] = 8
      elif(temp.loc[i,'cum rank2'] ) < decile*9:
        temp.loc[i,'income decile2'] = 9
      elif(temp.loc[i,'cum rank2'] ) < decile*10:
        temp.loc[i,'income decile2'] = 10
      elif(temp.loc[i,'cum rank2'] ) < decile*11:
        temp.loc[i,'income decile2'] = 11
      elif(temp.loc[i,'cum rank2'] ) < decile*12:
        temp.loc[i,'income decile2'] = 12
      elif(temp.loc[i,'cum rank2'] ) < decile*13:
        temp.loc[i,'income decile2'] = 13
      elif(temp.loc[i,'cum rank2'] ) < decile*14:
        temp.loc[i,'income decile2'] = 14
      elif(temp.loc[i,'cum rank2'] ) < decile*15:
        temp.loc[i,'income decile2'] = 15
      elif(temp.loc[i,'cum rank2'] ) < decile*16:
        temp.loc[i,'income decile2'] = 16
      elif(temp.loc[i,'cum rank2'] ) < decile*17:
        temp.loc[i,'income decile2'] = 17
      elif(temp.loc[i,'cum rank2'] ) < decile*18:
        temp.loc[i,'income decile2'] = 18
      elif(temp.loc[i,'cum rank2'] ) < decile*19:
        temp.loc[i,'income decile2'] = 19
      else:
        temp.loc[i,'income decile2'] = 20
        
         
    for i in range(0,len(ghg_mults)):
      temp[ghg_mults.index[i]+' ghg'] = temp.iloc[:,23+i]*ghg_mults.loc[ghg_mults.index[i],yr]*52/1000
    
                                     
    data = np.dot(temp.loc[:,'1.1.1 Bread and cereals ghg':'12.7.1 Other services n.e.c. ghg'].values,np.transpose(agg.values))
    data = df(data, index = temp.index, columns = agg.index)
    data['case'] = temp['case']
    data['weight'] = temp['weight']
    data['pop x weight'] = temp['weight']*temp['no_people']
    data['sex HRP'] = temp['sex hrp']
    data['income decile'] = temp['income decile']
    data['income decile2'] = temp['income decile2']
    data['GOR'] = temp['GOR']
    data['OA class 1'] = temp['OA class 1']
    data['OECD scale'] = temp['OECD scale']
    data['income'] = temp['Income anonymised']
    data['income equiv'] = temp['income equiv']
    
    footdata[yr] = data  

  return footdata

#################################
# used in generations_2023_main #
#################################

# Also used in defra_main_2023 (see defra_main_2023):
# convert43to41
# convert_hhspend_sizes
# convert_exp_tot_sizes
# make_balanced_totals
# make_new_Y_109
# make_totals
# make_y_hh_109
# removeoutliers

def processdataforfoots(ghg_mults,hhspenddata3,years):  
  footdata={}
  
  agg_filepath = 'C:/Users/earao/OneDrive - University of Leeds/Projects/milena/'
  agg = pd.read_excel(os.path.join(agg_filepath, 'aggregation.xlsx'), sheet_name='Sheet1', header = 0, index_col=0)
 
  
  for yr in years:
    temp = hhspenddata3[yr]
    
    temp = temp.reset_index()
    
    for i in range(0,len(temp)):
      if(temp.loc[i,'age hrp']<30):
         temp.loc[i,'age code'] = 1
      elif(temp.loc[i,'age hrp']<45):
         temp.loc[i,'age code'] = 2
      elif(temp.loc[i,'age hrp']<60):
         temp.loc[i,'age code'] = 3
      elif(temp.loc[i,'age hrp']<75):
         temp.loc[i,'age code'] = 4
      else:
         temp.loc[i,'age code'] = 5
    
    
    temp = temp.sort_values(by = "Income anonymised")
    
    temp = temp.reset_index()
    
    for i in range(0,len(temp)):
      if i == 0:
        temp.loc[i,'cum rank'] = temp.loc[i,'weight']
      else:
        temp.loc[i,'cum rank'] = temp.loc[i-1,'cum rank']+temp.loc[i,'weight']
        
    decile = np.sum(temp['weight'])/20
    
    for i in range(0,len(temp)):
      if(temp.loc[i,'cum rank'] ) < decile:
        temp.loc[i,'income decile'] = 1
      elif(temp.loc[i,'cum rank'] ) < decile*2:
        temp.loc[i,'income decile'] = 2
      elif(temp.loc[i,'cum rank'] ) < decile*3:
        temp.loc[i,'income decile'] = 3
      elif(temp.loc[i,'cum rank'] ) < decile*4:
        temp.loc[i,'income decile'] = 4
      elif(temp.loc[i,'cum rank'] ) < decile*5:
        temp.loc[i,'income decile'] = 5
      elif(temp.loc[i,'cum rank'] ) < decile*6:
        temp.loc[i,'income decile'] = 6
      elif(temp.loc[i,'cum rank'] ) < decile*7:
        temp.loc[i,'income decile'] = 7
      elif(temp.loc[i,'cum rank'] ) < decile*8:
        temp.loc[i,'income decile'] = 8
      elif(temp.loc[i,'cum rank'] ) < decile*9:
        temp.loc[i,'income decile'] = 9
      elif(temp.loc[i,'cum rank'] ) < decile*10:
        temp.loc[i,'income decile'] = 10
      elif(temp.loc[i,'cum rank'] ) < decile*11:
        temp.loc[i,'income decile'] = 11
      elif(temp.loc[i,'cum rank'] ) < decile*12:
        temp.loc[i,'income decile'] = 12
      elif(temp.loc[i,'cum rank'] ) < decile*13:
        temp.loc[i,'income decile'] = 13
      elif(temp.loc[i,'cum rank'] ) < decile*14:
        temp.loc[i,'income decile'] = 14
      elif(temp.loc[i,'cum rank'] ) < decile*15:
        temp.loc[i,'income decile'] = 15
      elif(temp.loc[i,'cum rank'] ) < decile*16:
        temp.loc[i,'income decile'] = 16
      elif(temp.loc[i,'cum rank'] ) < decile*17:
        temp.loc[i,'income decile'] = 17
      elif(temp.loc[i,'cum rank'] ) < decile*18:
        temp.loc[i,'income decile'] = 18
      elif(temp.loc[i,'cum rank'] ) < decile*19:
        temp.loc[i,'income decile'] = 19
      else:
        temp.loc[i,'income decile'] = 20
        
    temp['income equiv'] = temp['Income anonymised']/temp['OECD scale']*2
    
    temp = temp.sort_values(by = "income equiv")
    
    temp = temp.reset_index()
    
    for i in range(0,len(temp)):
      if i == 0:
        temp.loc[i,'cum rank2'] = temp.loc[i,'weight']
      else:
        temp.loc[i,'cum rank2'] = temp.loc[i-1,'cum rank2']+temp.loc[i,'weight']
        
    decile = np.sum(temp['weight'])/20
    
    for i in range(0,len(temp)):
      if(temp.loc[i,'cum rank2'] ) < decile:
        temp.loc[i,'income decile2'] = 1
      elif(temp.loc[i,'cum rank2'] ) < decile*2:
        temp.loc[i,'income decile2'] = 2
      elif(temp.loc[i,'cum rank2'] ) < decile*3:
        temp.loc[i,'income decile2'] = 3
      elif(temp.loc[i,'cum rank2'] ) < decile*4:
        temp.loc[i,'income decile2'] = 4
      elif(temp.loc[i,'cum rank2'] ) < decile*5:
        temp.loc[i,'income decile2'] = 5
      elif(temp.loc[i,'cum rank2'] ) < decile*6:
        temp.loc[i,'income decile2'] = 6
      elif(temp.loc[i,'cum rank2'] ) < decile*7:
        temp.loc[i,'income decile2'] = 7
      elif(temp.loc[i,'cum rank2'] ) < decile*8:
        temp.loc[i,'income decile2'] = 8
      elif(temp.loc[i,'cum rank2'] ) < decile*9:
        temp.loc[i,'income decile2'] = 9
      elif(temp.loc[i,'cum rank2'] ) < decile*10:
        temp.loc[i,'income decile2'] = 10
      elif(temp.loc[i,'cum rank2'] ) < decile*11:
        temp.loc[i,'income decile2'] = 11
      elif(temp.loc[i,'cum rank2'] ) < decile*12:
        temp.loc[i,'income decile2'] = 12
      elif(temp.loc[i,'cum rank2'] ) < decile*13:
        temp.loc[i,'income decile2'] = 13
      elif(temp.loc[i,'cum rank2'] ) < decile*14:
        temp.loc[i,'income decile2'] = 14
      elif(temp.loc[i,'cum rank2'] ) < decile*15:
        temp.loc[i,'income decile2'] = 15
      elif(temp.loc[i,'cum rank2'] ) < decile*16:
        temp.loc[i,'income decile2'] = 16
      elif(temp.loc[i,'cum rank2'] ) < decile*17:
        temp.loc[i,'income decile2'] = 17
      elif(temp.loc[i,'cum rank2'] ) < decile*18:
        temp.loc[i,'income decile2'] = 18
      elif(temp.loc[i,'cum rank2'] ) < decile*19:
        temp.loc[i,'income decile2'] = 19
      else:
        temp.loc[i,'income decile2'] = 20
        
    for i in range(0,len(temp)):
      temp.loc[i,'birth year'] = yr - temp.loc[i,'age hrp']
      
      # greatest gen 1901-27 = 1, silent gen 1928-45 = 2, babyboomer 1946-64 = 3
      # gen x 1965-80 = 4, millenial 1981-96 = 5, gen z 1997-2012 = 6

    for i in range(0,len(temp)):
      if(temp.loc[i,'birth year']>1981):
        temp.loc[i,'gen code'] = 3
      elif(temp.loc[i,'birth year']>1965):
        temp.loc[i,'gen code'] = 2
      elif(temp.loc[i,'birth year']>1946):
        temp.loc[i,'gen code'] = 1
      else:
        temp.loc[i,'gen code'] = 0
       
    for i in range(0,len(ghg_mults)):
      temp[ghg_mults.index[i]+' ghg'] = temp.iloc[:,23+i]*ghg_mults.loc[ghg_mults.index[i],yr]*52/1000
    
                                     
    data = np.dot(temp.loc[:,'1.1.1 Bread and cereals ghg':'12.7.1 Other services n.e.c. ghg'].values,np.transpose(agg.values))
    data = df(data, index = temp.index, columns = agg.index)
    data['case'] = temp['case']
    data['weight'] = temp['weight']
    data['pop x weight'] = temp['weight']*temp['no_people']
    data['age code'] = temp['age code']
    data['gen code'] = temp['gen code']
    data['sex HRP'] = temp['sex hrp']
    data['income decile'] = temp['income decile']
    data['income decile2'] = temp['income decile2']
    data['GOR'] = temp['GOR']
    data['OA class 1'] = temp['OA class 1']
    data['OECD scale'] = temp['OECD scale']
    
    footdata[yr] = data  

  return footdata

def totalbygen(footdata,years):  
  totalbygen = {} 
  for yr in years:  
    temp = footdata[yr].iloc[:,0:47].groupby(by='gen code').sum()    
    totalbygen[yr] = temp    
  return totalbygen

def hhagecohortfoots(footdata,years):  
  hh_agecohort_foots = {} 
  for yr in years:  
    temp = footdata[yr].iloc[:,0:47].groupby(by='gen code').sum()    
    hh_agecohort_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['weight'],[46,1])))      
  return hh_agecohort_foots

def popagecohortfoots(footdata,years):  
  pop_agecohort_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:46].groupby(by='gen code').sum()    
    pop_agecohort_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['pop x weight'],[45,1])))
  return pop_agecohort_foots

def hhagefoots(footdata,years):  
  hh_age_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:44].groupby(by='age code').sum()    
    hh_age_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['weight'],[43,1])))      
  return hh_age_foots

def popagefoots(footdata,years):  
  pop_age_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:44].groupby(by='age code').sum()    
    pop_age_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['pop x weight'],[43,1])))
  return pop_age_foots

def hhincfoots(footdata,years):  
  hh_income_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:44].groupby(by='income decile').sum()    
    hh_income_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['weight'],[43,1])))
  return hh_income_foots

def popincfoots(footdata,years):  
  pop_income_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:44].groupby(by='income decile').sum()    
    pop_income_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['pop x weight'],[43,1])))
  return pop_income_foots

def hhinc2foots(footdata,years):  
  hh_income2_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:44].groupby(by='income decile2').sum()    
    hh_income2_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['weight'],[43,1])))
  return hh_income2_foots

def popinc2foots(footdata,years):  
  pop_income2_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:44].groupby(by='income decile2').sum()    
    pop_income2_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['pop x weight'],[43,1])))
  return pop_income2_foots

def hhregfoots(footdata,years):  
  hh_reg_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:45].groupby(by='GOR').sum()    
    hh_reg_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['weight'],[44,1])))
  return hh_reg_foots

def popregfoots(footdata,years):  
  pop_reg_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:45].groupby(by='GOR').sum()    
    pop_reg_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['pop x weight'],[44,1])))
  return pop_reg_foots

def hhoacfoots(footdata,years):  
  hh_oac_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:46].groupby(by='OA class 1').sum()    
    hh_oac_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['weight'],[45,1])))
  return hh_oac_foots

def popoacfoots(footdata,years):  
  pop_oac_foots = {}  
  for yr in years:    
    temp = footdata[yr].iloc[:,0:46].groupby(by='OA class 1').sum()    
    pop_oac_foots[yr] = np.divide(temp,np.transpose(np.tile(temp['pop x weight'],[45,1])))
  return pop_oac_foots

#############################################
# used in sub_national_footprints_2023_main #
#############################################

# Also used in defra_main_2023 (see defra_main_2023):
# convert43to41
# convert_hhspend_sizes
# convert_exp_tot_sizes
# make_balanced_totals
# make_new_Y_109
# make_totals
# make_y_hh_109
# removeoutliers
# make_y_countries_2023	
# make_y_regions_2023	
    
# Also used in generations_2023_main (see generations_2023_main):
# hhagefoots	
# hhincfoots	
# hhinc2foots	
# hhoacfoots	
# hhregfoots
# popagefoots	
# popincfoots	
# popinc2foots	
# popoacfoots	
# popregfoots	
# processdataforfoots	

def add_pop_factors(reglaspropyr,regpophholdsyr,newY,regions,years):
 
 for yr in years:
  
  for r in regions:
   
   pop_reg = np.sum(regpophholdsyr[yr][r].loc[:,'pop'])
   
   popprop = np.zeros((len(regpophholdsyr[yr][r].loc[:,'pop'])))
   
   for i, la in enumerate(regpophholdsyr[yr][r].index):
    pop_la = regpophholdsyr[yr][r].loc[la,'pop']
    
    popprop[i] = pop_la/pop_reg
   
   reglaspropyr[yr][r][newY[2001].columns[106]] = popprop
   reglaspropyr[yr][r][newY[2001].columns[107]] = popprop
   reglaspropyr[yr][r][newY[2001].columns[108]] = popprop
   reglaspropyr[yr][r][newY[2001].columns[109]] = popprop
   reglaspropyr[yr][r][newY[2001].columns[110]] = popprop
   reglaspropyr[yr][r][newY[2001].columns[111]] = popprop
 
 
 return reglaspropyr

def correct_la_spend_props_gas_elec_2023(reglaspropyr,regions,regions_lc,years,wd):
  
  filepath = wd + 'data/raw/BEIS energy/'
  file = os.path.join(filepath, 'Sub-national_energy_consumption_statistics_2005-2020.xlsx')
  for i,r in enumerate(regions_lc):
    gas_prop = pd.read_excel(file, sheet_name = [r+'_prop_gas'], index_col = 0)
    elc_prop = pd.read_excel(file, sheet_name = [r+'_prop_elec'], index_col = 0)
    
    for yr in years:
      reglaspropyr[yr][regions[i]].iloc[:,30] = elc_prop[r+'_prop_elec'][yr]
      reglaspropyr[yr][regions[i]].iloc[:,31] = gas_prop[r+'_prop_gas'][yr]
      
  return reglaspropyr

def correct_reglaspendyr_zero(reglaspendyr,regions,years):
  
  for yr in years:
    for r in regions:
      total = np.sum(reglaspendyr[yr][r].values,1)
      for s in reglaspendyr[yr][r].columns:
        if np.sum(reglaspendyr[yr][r][s])==0:
          reglaspendyr[yr][r][s] = total 
      
  return reglaspendyr

def make_la_spend_props_by_region_year(reglaspendyr,regions,years):
  
  reglaspropyr = {}

  for yr in years:
    regpropspend = {}
  
    for r in regions:
      total = np.sum(reglaspendyr[yr][r])
      prop = np.divide(reglaspendyr[yr][r], np.tile(total,[len(reglaspendyr[yr][r]),1]))
      prop=prop.fillna(0)
    
      regpropspend[r] = prop
  
    reglaspropyr[yr] = regpropspend

  return reglaspropyr

def make_la_spends_pop_by_region_year(regoacsyr,oacyrspends,regions,years):
  
  reglaspendyr = {}
  regpophholdsyr = {}

  for yr in years:
    reglaspend = {}
    regpophholds = {}
    
    for r in regions:
      ladspend = df(np.zeros(shape=(len(regoacsyr[yr][r]),105)),index=regoacsyr[yr][r].index,columns=oacyrspends[yr][r].columns)
      for n in range(0,len(regoacsyr[yr][r])):
        ladspend.iloc[n] = np.transpose(regoacsyr[yr][r]['hholds'][n]*oacyrspends[yr][r].loc[regoacsyr[yr][r]['OAC_code'][n].lower()])
      ladspend['LAD_code'] = ladspend.index.get_level_values(0)
      ladspend = ladspend.groupby(['LAD_code']).sum()
      ladpophholds = regoacsyr[yr][r].groupby(['LAD_code']).sum()
      
      reglaspend[r] = ladspend
      regpophholds[r] = ladpophholds
      
    reglaspendyr[yr] = reglaspend
    regpophholdsyr[yr] = regpophholds

  return (reglaspendyr,regpophholdsyr) 
	
def make_pop_hhold_by_oac_region_year(oa_lookup01,oa_lookup11,oa_lookup21,census_filepath,regions):
  regoacsyr = {}
  
  file = os.path.join(census_filepath, '2001OAC.csv') 
  oac01 = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
  file = os.path.join(census_filepath, '2001OACs.csv') 
  oac02 = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
  file = os.path.join(census_filepath, '2001OACw.csv') 
  oac03 = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
  file = os.path.join(census_filepath, '2001OACn.csv') 
  oac04 = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
  oac=pd.concat([oac01,oac02,oac03,oac04])
  oac2001pophh = oac.join(oa_lookup01)
    
  for yr in range(2001,2014):
    regoacs = {}
    for r in regions:
      if yr == 2001:
        temp = oac2001pophh[oac2001pophh['REG'] == r].groupby(['LAD','AC Subgroup']).sum()
        temp = temp.rename(columns={'uv0010001': 'pop','ks0200001':'hholds'})
        temp['LAD_code'] = temp.index.get_level_values(0)
        temp['OAC_code'] = temp.index.get_level_values(1)
        regoacs[r] = temp
      else:
        temp = oac2001pophh[oac2001pophh['REG'] == r].groupby(['LAD','AC Subgroup']).sum()
        temp = temp.rename(columns={'uv0010001': 'pop','ks0200001':'hholds'})
        temp['LAD_code'] = temp.index.get_level_values(0)
        temp['OAC_code'] = temp.index.get_level_values(1)
      
        file = os.path.join(census_filepath, 'MYEB3_summary_components_of_change_series_UK_(2018_geog19).csv') 
        data = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
        data['popprop'] = data['population_' +str(yr)]/data['population_2001']
        data['LAD_code'] = data.index
        fred = temp.merge(data,how='left')
        fred.index = temp.index
        temp['pop'] = temp['pop']*fred['popprop']
        temp['hholds'] = temp['hholds']*fred['popprop']
        
        regoacs[r] = temp
    regoacsyr[yr] = regoacs
      
  regions2 = ['north-east','north-west','yorkshire-and-the-humber','east-midlands','west-midlands','east','london','south-east','south-west','scotland','wales']
  
  file = os.path.join(census_filepath, 'oa2011pop.csv') 
  pop_e = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
  file = os.path.join(census_filepath, 'oa2011pops.csv') 
  pop_s = pd.read_csv(file, index_col = 0) 
  pop = pd.concat([pop_e,pop_s])
  file = os.path.join(census_filepath, 'oa2011hhold.csv') 
  hhold_e = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
  file = os.path.join(census_filepath, 'oa2011hholds.csv') 
  hhold_s = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
  hhold = df(pd.concat([hhold_e,hhold_s]))
  file = os.path.join(census_filepath, '2011 OAC Clusters and Names csv v2.csv') 
  oac11 = pd.read_csv(file, encoding="iso8859_15", index_col = 0) 
   
  oac2011pophh = oac11.join(pop)
  oac2011pophh = oac2011pophh.join(hhold)
  oac2011pophh = oac2011pophh.join(oa_lookup11) 
        
  for yr in range(2014,2021):
  
    regoacs = {}
    print(yr)
          
    for i,r in enumerate(regions):
       print(r)
       temp = oac2011pophh[oac2011pophh['Region/Country Name'] == r]
       
       if r == 'Scotland':
         
         f='sape-' +str(yr) +'-persons.xlsx'
         file = os.path.join(census_filepath, f)
         sn = 'Table 1a Persons ('+str(yr) +')'
         data = pd.read_excel(file,sheet_name = sn, skiprows = 5, index_col=0) 
         newpop = df(data['Unnamed: 3'])
         LSOApop = temp.groupby('LSOA').sum()
         newpop2 = newpop.join(LSOApop)
         popprop = df(newpop2['Unnamed: 3']/newpop2['pop'],columns = ['popprop'])
         fred = temp.join(popprop,on='LSOA')
         temp['pop'] = temp['pop']*fred['popprop']
         temp['hholds'] = temp['hholds']*fred['popprop']
         
         temp2 = temp.groupby(['LAD','Subgroup Code']).sum()
         temp2['LAD_code'] = temp2.index.get_level_values(0)
         temp2['OAC_code'] = temp2.index.get_level_values(1)
         temp2 = temp2.rename(columns={'Supergroup Code': 'AC supergroup'})
          
       else:
         
         f='mid-' +str(yr) +'-coa-unformatted-syoa-estimates-' +regions2[i] +'.xlsx'
         file = os.path.join(census_filepath, f) 
         sn = 'Mid-'+str(yr) +' Persons'
         data = pd.read_excel(file,sheet_name = sn,index_col=0,header=4) 
         newpop = df(data['All Ages'])
         
         fred = temp.join(newpop)
         fred = fred.replace(0, 1)
         fred['popprop'] = fred['All Ages']/fred['pop']
         temp['pop'] = temp['pop']*fred['popprop']
         temp['hholds'] = temp['hholds']*fred['popprop']
         
         temp2 = temp.groupby(['LAD','Subgroup Code']).sum()
         temp2['LAD_code'] = temp2.index.get_level_values(0)
         temp2['OAC_code'] = temp2.index.get_level_values(1)
         temp2 = temp2.rename(columns={'Supergroup Code': 'AC supergroup'})
         
       regoacs[r] = temp2
         
    regoacsyr[yr] = regoacs
   
  yr = 2021
  
  f = 'OA_(2011)_to_OA_(2021)_to_Local_Authority_District_(2022)_for_England_and_Wales_Lookup_(Version_2).csv'
  file = os.path.join(census_filepath, f)
  lookup = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
  lookup.set_index(['OA11CD'],inplace=True)
  popfile = os.path.join(census_filepath, 'oa2021pop.csv')
  popew = pd.read_csv(popfile, encoding="iso8859_15", index_col = 0)
  hhfile = os.path.join(census_filepath, 'oa2021hhold.csv')
  hhew = pd.read_csv(hhfile,index_col = 0)
  
  oac2021pophh =  oac2011pophh.join(lookup)
  oac2021pophh = oac2021pophh.drop(['pop','hholds'],axis=1) 
  oac2021pophh =  oac2021pophh.merge(popew, left_on = 'OA21CD', right_index=True)
  oac2021pophh =  oac2021pophh.merge(hhew, left_on = 'OA21CD', right_index=True)
  
  for i,r in enumerate(regions):
     print(r)
     temp = oac2021pophh[oac2021pophh['Region/Country Name'] == r]
     
     if r == 'Scotland':
                
       temp = oac2011pophh[oac2011pophh['Region/Country Name'] == r]
             
       f='sape-' +str(yr) +'-persons.xlsx'
       file = os.path.join(census_filepath, f)
       sn = 'Table 1a Persons ('+str(yr) +')'
       data = pd.read_excel(file,sheet_name = sn, skiprows = 5, index_col=0) 
       newpop = df(data['Unnamed: 3'])
       LSOApop = temp.groupby('LSOA').sum()
       newpop2 = newpop.join(LSOApop)
       popprop = df(newpop2['Unnamed: 3']/newpop2['pop'],columns = ['popprop'])
       fred = temp.join(popprop,on='LSOA')
       temp['pop'] = temp['pop']*fred['popprop']
       temp['hholds'] = temp['hholds']*fred['popprop']
       
       temp2 = temp.groupby(['LAD','Subgroup Code']).sum()
       temp2['LAD_code'] = temp2.index.get_level_values(0)
       temp2['OAC_code'] = temp2.index.get_level_values(1)
       temp2 = temp2.rename(columns={'Supergroup Code': 'AC supergroup'})
        
     else:
          
       temp2 = temp.groupby(['LAD','Subgroup Code']).sum()
       temp2['LAD_code'] = temp2.index.get_level_values(0)
       temp2['OAC_code'] = temp2.index.get_level_values(1)
       temp2 = temp2.rename(columns={'Supergroup Code': 'AC supergroup'})
           
       
     regoacs[r] = temp2
    
  regoacsyr[yr] = regoacs
    

  return regoacsyr
	
def make_y_regional(region_str, region_dict, y_regions, years):
  y_regional = {}
  reg_no = str(region_dict[region_str])
  for yr in years:
    y_regional[yr] = df(y_regions[str(yr) + '_' + reg_no].values,
              index=y_regions[str(yr) + '_' + reg_no].index,
              columns=y_regions[str(yr) + '_' + reg_no].columns)
  return y_regional

def makeoacspends(hhspenddata3,oac_none_years,oac_2001_years,oac_2011_years):
  
  oacyrspends = {}
  oacyrmeta = {}
  
  regions = ['North East','North West','Yorkshire and The Humber','East Midlands','West Midlands','East',
             'London','South East','South West','Scotland','Wales']
  
  for yr in oac_2001_years.tolist() + oac_2011_years.tolist():
    regspends = {}
    regmeta = {}
    
    # create OAC lookup from scratch to ensure all categories are included
    lookup = df()
    if yr in oac_2001_years.tolist():
        lookup['OA class 3'] = ['1a1','1a2','1a3','1b1','1b2','1c1','1c2','1c3',
                                '2a1','2a2','2b1','2b2', 
                                '3a1','3a2','3b1','3b2','3c1','3c2',
                                '4a1','4a2','4b1','4b2','4b3','4b4','4c1','4c2','4c3','4d1','4d2',
                                '5a1','5a2','5b1','5b2','5b3','5b4','5c1','5c2','5c3',
                                '6a1','6a2','6b1','6b2','6b3','6c1','6c2','6d1','6d2',
                                '7a1','7a2','7a3','7b1','7b2']
    else:
        lookup['OA class 3'] = ['1a1', '1a2', '1a3', '1a4', '1b1', '1b2', '1b3', '1c1', '1c2', '1c3',
                                '2a1', '2a2', '2a3', '2b1', '2b2', '2c1', '2c2', '2c3', '2d1', '2d2', '2d3',
                                '3a1', '3a2', '3b1', '3b2', '3b3', '3c1', '3c2', '3d1', '3d2', '3d3',
                                '4a1', '4a2', '4a3', '4b1', '4b2', '4c1', '4c2', '4c3',
                                '5a1', '5a2', '5a3', '5b1', '5b2', '5b3',
                                '6a1', '6a2', '6a3', '6a4', '6b1', '6b2', '6b3', '6b4',
                                '7a1', '7a2', '7a3', '7b1', '7b2', '7b3', '7c1', '7c2', '7c3', '7d1', '7d2', '7d3', '7d4',
                                '8a1', '8a2', '8b1', '8b2', '8c1', '8c2', '8c3', '8d1', '8d2', '8d3']
    lookup['OA class 2'] = lookup['OA class 3'].str[:2]
    lookup['OA class 1'] = lookup['OA class 3'].str[:1]
    
    hhspenddata3[yr][['OA class 1', 'OA class 2', 'OA class 3']] = hhspenddata3[yr][['OA class 1', 'OA class 2', 'OA class 3']].apply(lambda x:x.str.lower())                       
    
    for r in range(len(regions)):
        reg = r+1
      
        regavspend = df(np.sum(hhspenddata3[yr][hhspenddata3[yr]["GOR"] == reg].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'])/np.sum(hhspenddata3[yr][hhspenddata3[yr]["GOR"] == reg].loc[:,'weight']))
        
        count = {}; spendbyoacreg = {}
        lookup_temp = cp.copy(lookup)
        for i in [1, 2, 3]:
            
            count[i] = df(hhspenddata3[yr].loc[hhspenddata3[yr]["GOR"] == reg].groupby(['OA class '  + str(i)]).count().iloc[:,0])
            count[i].columns = ['oac' + str(i) + '_count']
            
            temp = hhspenddata3[yr].loc[hhspenddata3[yr]["GOR"] == reg].groupby(['OA class ' + str(i)]).sum()
            spendbyoacreg[i] = temp.loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.']/np.transpose(np.tile(temp.loc[:,'weight'].values,(105,1)))

            lookup_temp = lookup_temp.merge(count[i].reset_index(), on='OA class ' + str(i), how='left').fillna(0).rename(columns={'OA class ' + str(i):'oac' + str(i)})
        
        lookup_temp = lookup_temp.set_index(['oac1', 'oac2', 'oac3'])
        temp = df(columns=spendbyoacreg[1].columns, index=lookup_temp.index).fillna(0)
      
        for o in temp.index:
            if lookup_temp.loc[o,'oac3_count']>19:
                temp.loc[o] = spendbyoacreg[3].loc[o[2]]
            elif lookup_temp.loc[o,'oac2_count']>19:
                temp.loc[o] = spendbyoacreg[2].loc[o[1]]
            elif lookup_temp.loc[o,'oac1_count']>19:
                temp.loc[o] = spendbyoacreg[1].loc[o[0]]
            else:
                temp.loc[o] = regavspend.iloc[:,0]
        
        temp = temp.droplevel(axis=0, level=['oac1', 'oac2'])
        regspends[regions[r]] = temp
        regmeta[regions[r]] = lookup_temp
        
    oacyrspends[yr] = regspends
    oacyrmeta[yr] = regmeta
    
  for yr in oac_none_years:
    
    oacyrspends[yr] = oacyrspends[2007]
    oacyrmeta[yr] = oacyrmeta[2007]
    
  return(oacyrspends,oacyrmeta)
