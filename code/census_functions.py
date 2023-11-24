#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:04:07 2018
This function is used with household_demand_main.py
@author: earao
"""
import pandas as pd
import numpy as np
import os
df = pd.DataFrame

#############################################
# used in sub_national_footprints_2023_main #
#############################################

def make_11_lookup(census_filepath):    
    
    file = os.path.join(census_filepath, 'Output_Area_to_Lower_Layer_Super_Output_Area_to_Middle_Layer_Super_Output_Area_to_Local_Authority_District_2020_Lookup_in_Great_Britain__Classification_Version_2.csv')  
    temp = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    oa_lookup11=temp.drop(['OAC11CD','OAC11NM','SOAC11CD','SOAC11NM','LSOA11NM','MSOA11NM','LACCD','LACNM', 'RGN11CD', 'CTRY11CD','CTRY11NM','FID'], axis=1)
    oa_lookup11=oa_lookup11.rename(index=str, columns={'LSOA11CD': 'LSOA', "MSOA11CD": "MSOA", 'LAD20NM':'LAD_nm',"LAD20CD": "LAD", "RGN11NM":"REG"})

    return oa_lookup11
  
def make_01_lookup(census_filepath,oa_lookup11):
    
    file = os.path.join(census_filepath, 'OA01_LSOA01_MSOA01_EW_LU.csv')  
    ew01 = pd.read_csv(file, encoding="iso8859_15", index_col = 0)    
    file = os.path.join(census_filepath, 'OA01_OA11_LAD20_EW_LU.csv')  
    lookupLAD = pd.read_csv(file, encoding="iso8859_15", index_col = 1)
    file = os.path.join(census_filepath, 'laregionlookup2020.xls')  
    lookupREG = pd.read_excel(file, header = 4, index_col=0)
    ew01 = ew01.join(lookupLAD)
    ew01 = ew01.join(lookupREG, on='LAD20CD')
    ew01=ew01.rename(index=str, columns={'LSOA01CD': 'LSOA', "MSOA01CD": "MSOA", 'LAD20NM':'LAD_nm',"LAD20CD": "LAD", "region_name":"REG"})
    ew01=ew01.drop(['LSOA01NM','MSOA01NM','OA01CD','OA11CD','CHGIND','la_name'], axis=1)

    file = os.path.join(census_filepath, 'scot2001_lookup.csv') 
    s01 = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    file = os.path.join(census_filepath, 'scotlandLAD.csv')  
    lookupLAD = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    s01 = s01.join(lookupLAD, on='council_area')
    s01['REG'] = 'Scotland'
    s01=s01.drop(['council_area'], axis=1)
    s01=s01.rename(index=str, columns={"data_zone": "LSOA", "inter_zone": "MSOA", "Council Area NAME": "LAD_nm",'council_area': 'LAD', })
        
    oa_lookup01 = pd.concat([ew01,s01], sort=True)
    
    return oa_lookup01

##################
# Not sorted yet #
##################

def get_data():
    oa_foot = {}
    oa_pop = {}
    
    oa_foot['2005'] = pd.read_csv('oa_foot_2005.csv', encoding="iso8859_15", index_col = 0)
    oa_foot['2007'] = pd.read_csv('oa_foot_2007.csv', encoding="iso8859_15", index_col = 0)
    oa_foot['2008'] = pd.read_csv('oa_foot_2008.csv', encoding="iso8859_15", index_col = 0)
    oa_foot['2012'] = pd.read_csv('oa_foot_2012.csv', encoding="iso8859_15", index_col = 0)
    oa_foot['2014'] = pd.read_csv('oa_foot_2014.csv', encoding="iso8859_15", index_col = 0)
    oa_foot['2016'] = pd.read_csv('oa_foot_2016.csv', encoding="iso8859_15", index_col = 0)
    
    oa_pop['2005'] = pd.read_csv('oa_pop_2005.csv', encoding="iso8859_15", index_col = 0)
    oa_pop['2007'] = pd.read_csv('oa_pop_2007.csv', encoding="iso8859_15", index_col = 0)
    oa_pop['2008'] = pd.read_csv('oa_pop_2008.csv', encoding="iso8859_15", index_col = 0)
    oa_pop['2012'] = pd.read_csv('oa_pop_2012.csv', encoding="iso8859_15", index_col = 0)
    oa_pop['2014'] = pd.read_csv('oa_pop_2014.csv', encoding="iso8859_15", index_col = 0)
    oa_pop['2016'] = pd.read_csv('oa_pop_2016.csv', encoding="iso8859_15", index_col = 0)

    return(oa_foot,oa_pop)


def make_01_11_lookup(census_lookup_path,oa_lookup11):

    file = os.path.join(census_lookup_path, 'oa2001_oldtonew_lookup.xlsx')  
    scot01_11a = pd.read_excel(file, header=0, index_col = None)
    file = os.path.join(census_lookup_path, 'geog-2011-cen-supp-info-2001-2011-master-postcode.xls')  
    scot01_11b = pd.read_excel(file, header=2, index_col = None)
    scot01_11 = pd.merge(scot01_11a,scot01_11b, how='left')
    scot01_11 = scot01_11.drop(columns=['MasterPostcode2001','OutputArea2001Code'])
    scot01_11 = scot01_11.set_index('OutputArea2011Code')
    scot01_11 = scot01_11.join(oa_lookup11)
    scot01_11 = scot01_11.rename(index=str, columns={"NRSoldOutputArea2001Code": 'OutputArea2001Code'})
    
    file = os.path.join(census_lookup_path, 'OA01_OA11_LAD11_EW_LU.csv')  
    ew01_11 = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew01_11 = ew01_11.drop(columns=['OA11CD', 'CHGIND', 'LAD11CD', 'LAD11NM', 'LAD11NMW'])
    ew01_11 = ew01_11.join(oa_lookup11, how='inner')
    ew01_11 = ew01_11.rename(index=str, columns={"OA01CDO": 'OutputArea2001Code'})
       
    oa01_11 = pd.concat([ew01_11,scot01_11], sort=True)
    
    return oa01_11

def get_census_pop(censusfilepath01,censusfilepath11):
    pop2001 = df
    pop2011 = df
    
    scotpop01 = df
    file = os.path.join(censusfilepath01, 'CAS001.csv')  
    temp = pd.read_csv(file, encoding="iso8859_15", index_col = None, header = 4)
    temp2 = np.zeros(shape = (42604,1))
    temp3 = []
    for a in range(0,42604):
        temp2[a,0] = temp.iloc[a*45,2]
        temp3.append(temp.iloc[a*45,0])
    scotpop01 = df(temp2, index = temp3)
    scotpop01.columns = ['Pop']
   
    file = os.path.join(censusfilepath11, 'Output Area std/KS101SC.csv')  
    scotpop11 = pd.read_csv(file, encoding="iso8859_15", thousands=',', index_col=0, header=4)
    scotpop11 = df(scotpop11.iloc[2:,0], index = scotpop11[2:].index)
    scotpop11.columns = ["Pop"]
    
    file = os.path.join(censusfilepath01, 'KS01.xlsx')  
    nirepop01 = pd.read_excel(file, index_col = 0)
    nirepop01 = df(nirepop01.iloc[2:,0])
    nirepop01.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'Key Statistics Tables (statistical geographies)/SMALL AREAS/KS101NIDATA0.csv')  
    nirepop11 = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    nirepop11 =  df(nirepop11.iloc[:,0])
    nirepop11.columns = ['Pop']
    
    file = os.path.join(censusfilepath01, '217934412.csv')  
    ew_pop01 = pd.read_csv(file, encoding="iso8859_15", index_col = 0, header = 5)
    ew_pop01 =  df(ew_pop01.iloc[:,1])
    ew_pop01.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_A.csv')  
    ew_Apop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Apop =  df(ew_Apop.iloc[:,0])
    ew_Apop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_B.csv')  
    ew_Bpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Bpop =  df(ew_Bpop.iloc[:,0])
    ew_Bpop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_D.csv')  
    ew_Dpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Dpop =  df(ew_Dpop.iloc[:,0])
    ew_Dpop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_E.csv')  
    ew_Epop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Epop =  df(ew_Epop.iloc[:,0])
    ew_Epop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_F.csv')  
    ew_Fpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Fpop =  df(ew_Fpop.iloc[:,0])
    ew_Fpop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_G.csv')  
    ew_Gpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Gpop =  df(ew_Gpop.iloc[:,0])
    ew_Gpop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_H.csv')  
    ew_Hpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Hpop =  df(ew_Hpop.iloc[:,0])
    ew_Hpop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_J.csv')  
    ew_Jpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Jpop =  df(ew_Jpop.iloc[:,0])
    ew_Jpop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_K.csv')  
    ew_Kpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Kpop =  df(ew_Kpop.iloc[:,0])
    ew_Kpop.columns = ['Pop']
    
    file = os.path.join(censusfilepath11, 'qs101ew_2011_oa/QS101EW_2011STATH_NAT_OA_REL_1.2.2/QS101EWDATA01_W.csv')  
    ew_Wpop = pd.read_csv(file, encoding="iso8859_15", index_col = 0)
    ew_Wpop =  df(ew_Wpop.iloc[:,0])
    ew_Wpop.columns = ['Pop']
    
    pop2001=pd.concat([scotpop01,nirepop01,ew_pop01],sort = True)
    pop2011=pd.concat([scotpop11,nirepop11,ew_Apop,ew_Bpop,ew_Dpop,ew_Epop,ew_Fpop,ew_Gpop,ew_Hpop,ew_Jpop,ew_Kpop,ew_Wpop],sort = True)
  
    return (pop2001,pop2011)
    

def remove_duplicates(oa_lookup01):
    
    temp = df(oa_lookup01.index.duplicated(keep = 'first'), index = oa_lookup01.index)
    temp = temp.rename(columns = {0:"boo"})
    temp2 = pd.concat([oa_lookup01,temp],axis=1,sort=False)
    temp3 = temp2[temp2['boo'] == False]
    temp3=temp3.drop(['boo'], axis=1)
    oa_lookup01 = df(temp3)
    
    return oa_lookup01


def make_regional_foots(region,oa_foot,oa_pop,oa_lookup11,oa_lookup01):
    pop_by_type = {}
    foot_by_type_tot = {}
    per_cap = {}
        
    pop_by_type['2016'] = oa_pop['2016'].join(oa_lookup11).groupby([region]).sum()
    pop_by_type['2014'] = oa_pop['2014'].join(oa_lookup11).groupby([region]).sum()
    pop_by_type['2012'] = oa_pop['2012'].join(oa_lookup11).groupby([region]).sum()
    pop_by_type['2008'] = oa_pop['2008'].join(oa_lookup01).groupby([region]).sum()
    pop_by_type['2007'] = oa_pop['2007'].join(oa_lookup01).groupby([region]).sum()
    pop_by_type['2005'] = oa_pop['2005'].join(oa_lookup01).groupby([region]).sum()
    foot_by_type_tot['2016'] = oa_foot['2016'].join(oa_lookup11).groupby([region]).sum()    
    foot_by_type_tot['2014'] = oa_foot['2014'].join(oa_lookup11).groupby([region]).sum()    
    foot_by_type_tot['2012'] = oa_foot['2012'].join(oa_lookup11).groupby([region]).sum()    
    foot_by_type_tot['2008'] = oa_foot['2008'].join(oa_lookup01).groupby([region]).sum()    
    foot_by_type_tot['2007'] = oa_foot['2007'].join(oa_lookup01).groupby([region]).sum()    
    foot_by_type_tot['2005'] = oa_foot['2005'].join(oa_lookup01).groupby([region]).sum()    
    per_cap['2016'] = np.divide(foot_by_type_tot['2016'],pop_by_type['2016'])
    per_cap['2014'] = np.divide(foot_by_type_tot['2014'],pop_by_type['2014'])
    per_cap['2012'] = np.divide(foot_by_type_tot['2012'],pop_by_type['2012'])
    per_cap['2008'] = np.divide(foot_by_type_tot['2008'],pop_by_type['2008'])
    per_cap['2007'] = np.divide(foot_by_type_tot['2007'],pop_by_type['2007'])
    per_cap['2005'] = np.divide(foot_by_type_tot['2005'],pop_by_type['2005'])
        
    return per_cap