# -*- coding: utf-8 -*-
"""
Spyder Editor

This function is used with ukmio_main.py to build the emissions extension vector

@authors: Anne Owen & Lena Kilian
"""
import numpy as np
import pandas as pd
import os
df = pd.DataFrame

############################
# used in ukmrio_main_2024 #
############################

def make_exio382_stressor(use,exioyrs,exiobase_filepath,meta,c_conc,i_conc): # used in ukmrio_main_2024

    exioCO2 = {}
    exioGHG = {}
       
    exio_i_index = []
    for c in range (0,len(c_conc.index)):
        for i in range (0,len(i_conc.index)):            
            exio_i_index.append("".join([i_conc.index[i], c_conc.index[c]]))
                     
    eind_slc = {}     
    for c in range(0,49):
        eind_slc[c] = slice(c*163,(c+1)*163)
      
    uind_slc = {}    
    for c in range(0,meta['reg']['len']):
         uind_slc[c] = slice(c*meta['use_dd']['len_col'],(c+1)*meta['use_dd']['len_col'])
          
    idata = np.zeros(shape = (len(exio_i_index),meta['use']['len_col']))
   
    for i in range (0,len(c_conc.index)):
        for c in range (0,len(c_conc.columns)):
            if c_conc.loc[c_conc.index[i]][c_conc.columns[c]] == 1:
                idata[i*len(i_conc):(i+1)*len(i_conc),uind_slc[c]] = i_conc.values
                
    exioUKconci = df(idata, index = exio_i_index, columns = meta['use']['col'])    
       
    for a in range(0,np.size(exioyrs)):
        
        filepath = exiobase_filepath + "3.8.2/MRSUT_{}/".format(str(exioyrs[a]))
        stressor = pd.read_csv(os.path.join(filepath, 'F.txt'), sep='\t', header = [0,1], index_col = 0)
        
        tempco2 = stressor.loc['Carbon dioxide (CO2) IPCC categories 1 to 4 and 6 to 7 (excl land use, land use change and forestry)']
        tempghg = stressor.loc['GHG emissions AR5 (GWP100) | GWP100 (IPCC, 2010)']
        
        weightedEXIOUKconci = np.zeros(shape=(7987,meta['use']['len_col']));
         
        uk_output_i = df.sum(use[exioyrs[a]].iloc[meta['v_d']['rng'],:], 1)
                             
        for m in range (0,49):
            for n in range (0,meta['reg']['len']):       
                if c_conc.iloc[m,n] == 1:
                    num = np.transpose(np.dot(np.diag(uk_output_i),np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])))
                    den = np.transpose(np.tile(np.dot(uk_output_i,np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])),(meta['use_dd']['len_col'],1)))                                      
                    weightedEXIOUKconci[eind_slc[m],uind_slc[n]] = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
        
        exioCO2[exioyrs[a]] = df(np.dot(tempco2,weightedEXIOUKconci), index = meta['v']['col'])
        exioGHG[exioyrs[a]] = df(np.dot(tempghg,weightedEXIOUKconci), index = meta['v']['col'])
        
    return(exioCO2,exioGHG)

def make_old_exio_stressor_382(exioCO2,exioGHG,edgar_filepath,regions): # used in ukmrio_main_2024
        
    file = os.path.join(edgar_filepath, 'v432_CO2_excl_short-cycle_org_C_1970_2012.xlsx')
    oldCO2prop = pd.read_excel(file, 'CO2_props',index_col = 0)
    
    oldGHGprop = oldCO2prop
    
    exioCO2[1994] = exioCO2[1995]*oldCO2prop.loc[:,1994].values
    exioGHG[1994] = exioGHG[1995]*oldGHGprop.loc[:,1994].values
    exioCO2[1993] = exioCO2[1994]*oldCO2prop.loc[:,1993].values
    exioGHG[1993] = exioGHG[1994]*oldGHGprop.loc[:,1993].values
    exioCO2[1992] = exioCO2[1993]*oldCO2prop.loc[:,1992].values
    exioGHG[1992] = exioGHG[1993]*oldGHGprop.loc[:,1992].values
    exioCO2[1991] = exioCO2[1992]*oldCO2prop.loc[:,1991].values
    exioGHG[1991] = exioGHG[1992]*oldGHGprop.loc[:,1991].values
    exioCO2[1990] = exioCO2[1991]*oldCO2prop.loc[:,1990].values
    exioGHG[1990] = exioGHG[1991]*oldGHGprop.loc[:,1990].values
        
    return (exioCO2,exioGHG)

def make_UK_emissions(ons_filepath,yrs): # used in ukmrio_main_2024
    
    emissionsyrs = np.array([int(x) for x in range(1990,2022)])
    
    file = os.path.join(ons_filepath, 'ONS environmental accounts/2024/atmoshpericemissionsghg.xlsx') # not my spelling!
    
    uk_ghg_sectors = pd.read_excel(file, 'GHG total ', usecols='C:AI', index_col=0, header=0, nrows=131, skiprows = 29)
    uk_co2_sectors = pd.read_excel(file, 'CO2', usecols='C:AI', index_col=0, header=0, nrows=131, skiprows = 29)
    uk_ghg_sectors.columns = emissionsyrs
    uk_co2_sectors.columns = emissionsyrs
    uk_ghg_direct = uk_ghg_sectors.iloc[129:131,0:32]
    uk_co2_direct = uk_co2_sectors.iloc[129:131,0:32]
         
    return (uk_ghg_sectors,uk_co2_sectors,uk_ghg_direct,uk_co2_direct) 

def make_co2_382(exioCO2,uk_co2_sectors,ons_filepath,yrs,meta): # used in ukmrio_main_2024
    
    co2 = {}
        
    file = os.path.join(ons_filepath, 'Analytical tables/sectorconc_112.xlsx')
    uk_emissions_conc = pd.read_excel(file, 'emissions',index_col = 0)
    
    for a in range(0,np.size(yrs)):
        temp = np.zeros(shape=(meta['v']['len']))
        exio_co2 = (exioCO2[yrs[a]])
        
        temp[meta['v_d']['rng']] = np.dot(np.transpose(uk_co2_sectors.loc[:,yrs[a]]),uk_emissions_conc)
        temp[meta['v_i']['rng']] = np.sum(exio_co2[meta['v_i']['rng']],1)
        
        co2[yrs[a]] = df(temp, index = meta['use']['col'])
    
    return(co2)
 
def make_ghg_382(exioGHG,uk_ghg_sectors,ons_filepath,yrs,meta): # used in ukmrio_main_2024
    
    ghg = {}
        
    file = os.path.join(ons_filepath, 'Analytical tables/sectorconc_112.xlsx')
    uk_emissions_conc = pd.read_excel(file, 'emissions',index_col = 0)
    
    for a in range(0,np.size(yrs)):
        temp = np.zeros(shape=(meta['v']['len']))
        exio_ghg = (exioGHG[yrs[a]])/1000000
        
        temp[meta['v_d']['rng']] = np.dot(np.transpose(uk_ghg_sectors.loc[:,yrs[a]]),uk_emissions_conc)
        temp[meta['v_i']['rng']] = np.sum(exio_ghg[meta['v_i']['rng']],1)  
    
        ghg[yrs[a]] = df(temp, index = meta['use']['col'])
        
    return(ghg)