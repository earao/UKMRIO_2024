# -*- coding: utf-8 -*-
"""
Spyder Editor

This function is used with ukmio_main.py to build the energy extension vector from IEA data

@author: earao
"""
import numpy as np
import pandas as pd
import os
df = pd.DataFrame

############################
# used in ukmrio_main_2024 #
############################

def make_IEA_data(iea_filepath): # used in ukmrio_main_2023
     
    filepath = iea_filepath + "Extended energy balances/2020/"
    
    aviation = pd.read_csv(os.path.join(filepath, 'aviation.csv'),encoding='latin-1', header = 0, index_col = 0)
    shipping = pd.read_csv(os.path.join(filepath, 'shipping.csv'),encoding='latin-1', header = 0, index_col = 0)
    iea_fe_data = pd.read_csv(os.path.join(filepath, 'iea_fe_1990_2020.csv'),encoding='latin-1', header = 0, index_col = [0,1])
    iea_fe_data=iea_fe_data.replace('x',0)
    iea_fe_data=iea_fe_data.replace('..',0)
    iea_fe_data=iea_fe_data.replace(' ',0)
              
    return (iea_fe_data,aviation,shipping)

def iea_to_full_exio382(inputs_filepath,exiobase_filepath,iea_fe_data,aviation,shipping): # used in ukmrio_main_2023
    
    fullexioNRG = {}
    
    c_conc_iea = pd.read_excel(os.path.join(inputs_filepath, 'EXIOBASE_data.xlsx'), sheet_name = 'iea_c_conc', index_col = 0, header = 0)
    s_conc_iea = pd.read_excel(os.path.join(inputs_filepath, 'EXIOBASE_data.xlsx'), sheet_name = 'iea_FE', index_col = [0,1], header = 0)
    
    exio_regions = c_conc_iea.index
    iea_regions = c_conc_iea.columns
    big_conc = np.zeros(shape = (163*49,27*147))
    
    for r in range(0,np.size(iea_regions)):
        for e in range(0,np.size(exio_regions)):
            if c_conc_iea.loc[exio_regions[e],iea_regions[r]] == 1:
                big_conc[e*163:(e+1)*163,r*27:(r+1)*27] = s_conc_iea
    
    eind_slc = {}
    
    for c in range(0,49):
         eind_slc[c] = slice(c*163,(c+1)*163)
    
    iind_slc = {}
    
    for c in range(0,170):
         iind_slc[c] = slice(c*27,(c+1)*27)
         
    road_prop = pd.read_excel(os.path.join(inputs_filepath,'Transport energy split 0816.xlsx'), sheet_name = 'for use', index_col = 0, header = 0)
      
    i_spend_on_e = np.zeros(shape=(26,7987)) 
    spend_on_a = np.zeros(shape=(26,49))    
    spend_on_s = np.zeros(shape=(26,4,49))           
    for yr in range(1995,2021):
        print(yr)
        filepath = exiobase_filepath + "3.8.2/MRSUT_{}/".format(str(yr))               
        exio_u = pd.read_csv(os.path.join(filepath, 'use.csv'), sep='\t', header = [0,1], index_col = [0,1])
        i_spend = np.zeros(shape=(200,7987))
        for e in range(0,np.size(exio_regions)): 
            i_spend = i_spend+exio_u.iloc[e*200:(e+1)*200,:].values
            spend_on_a[yr-1995,e] = np.sum(exio_u.iloc[:,e*163+124].values)
            spend_on_s[yr-1995,:,e] = np.sum(exio_u.iloc[:,e*163+121:e*163+125].values,0)
            
        i_spend_on_e = np.sum(i_spend[19:27,:],0)+np.sum(i_spend[63:85,:],0)+np.sum(i_spend[127:147,:],0)
        
        weightedconc = np.zeros(shape = (163*49,27*147))
               
        for r in range(0,np.size(iea_regions)):
            for e in range(0,np.size(exio_regions)):
                if c_conc_iea.loc[exio_regions[e],iea_regions[r]] == 1:
                    num = np.dot(np.diag(i_spend_on_e[eind_slc[e]]),(big_conc[eind_slc[e],iind_slc[r]]))
                    den = np.tile(np.dot(i_spend_on_e[eind_slc[e]],(big_conc[eind_slc[e],iind_slc[r]])),(163,1))                                     
                    weightedconc[eind_slc[e],iind_slc[r]] = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
    
        fullexioNRG[yr] = df(np.dot(weightedconc,iea_fe_data[str(yr)]))
        
        for e in range(0,np.size(exio_regions)):
            fullexioNRG[yr].iloc[e*163+120] = fullexioNRG[yr].iloc[e*163+120]*(1-road_prop.iloc[e,yr-1995])
            fullexioNRG[yr].iloc[e*163+124] = fullexioNRG[yr].iloc[e*163+124].values + aviation[str(yr)]*(spend_on_a[yr-1995,e]/np.sum(spend_on_a[yr-1995,:]))
            extra_shipping = shipping[str(yr)].values*(spend_on_s[yr-1995,:,e]/np.sum(spend_on_s[yr-1995,:,:]))
            for s in range(0,4):
                fullexioNRG[yr].iloc[e*163+121+s] = fullexioNRG[yr].iloc[e*163+121+s] + extra_shipping[s]
    
    for yr in range(1990,1996):       
        print(yr)
        filepath = exiobase_filepath + "3.8.2/MRSUT_1995/".format(str(yr))               
        exio_u = pd.read_csv(os.path.join(filepath, 'use.csv'), sep='\t', header = [0,1], index_col = [0,1])
        i_spend = np.zeros(shape=(200,7987))
        for e in range(0,np.size(exio_regions)): 
            i_spend = i_spend+exio_u.iloc[e*200:(e+1)*200,:].values
            spend_on_a[yr-1995,e] = np.sum(exio_u.iloc[:,e*163+124].values)
            spend_on_s[yr-1995,:,e] = np.sum(exio_u.iloc[:,e*163+121:e*163+125].values,0)
            
        i_spend_on_e = np.sum(i_spend[19:27,:],0)+np.sum(i_spend[63:85,:],0)+np.sum(i_spend[127:147,:],0)
        
        weightedconc = np.zeros(shape = (163*49,27*147))
                 
        for r in range(0,np.size(iea_regions)):
            for e in range(0,np.size(exio_regions)):
                if c_conc_iea.loc[exio_regions[e],iea_regions[r]] == 1:
                    num = np.dot(np.diag(i_spend_on_e[eind_slc[e]]),(big_conc[eind_slc[e],iind_slc[r]]))
                    den = np.tile(np.dot(i_spend_on_e[eind_slc[e]],(big_conc[eind_slc[e],iind_slc[r]])),(163,1))                                     
                    weightedconc[eind_slc[e],iind_slc[r]] = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
    
        fullexioNRG[yr] = df(np.dot(weightedconc,iea_fe_data[str(yr)]))
        
        for e in range(0,np.size(exio_regions)):
            fullexioNRG[yr].iloc[e*163+120] = fullexioNRG[yr].iloc[e*163+120]*(1-road_prop.iloc[e,yr-1995])
            fullexioNRG[yr].iloc[e*163+124] = fullexioNRG[yr].iloc[e*163+124].values + aviation[str(yr)]*(spend_on_a[yr-1995,e]/np.sum(spend_on_a[yr-1995,:]))
            extra_shipping = shipping[str(yr)].values*(spend_on_s[yr-1995,:,e]/np.sum(spend_on_s[yr-1995,:,:]))
            for s in range(0,4):
                fullexioNRG[yr].iloc[e*163+121+s] = fullexioNRG[yr].iloc[e*163+121+s] + extra_shipping[s]
       
    return (fullexioNRG)

def uk_exio_nrg(fullexioNRG,use,yrs,meta,c_conc,i_conc): # used in ukmrio_main_2023
    exioNRG = {}

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
    
    for a in range(0,np.size(yrs)):
        weightedEXIOUKconci = np.zeros(shape=(7987,meta['use']['len_col']));
         
        uk_output_i = df.sum(use[yrs[a]].iloc[meta['v_d']['rng'],:], 1)
                             
        for m in range (0,49):
            for n in range (0,meta['reg']['len']):       
                if c_conc.iloc[m,n] == 1:
                    num = np.transpose(np.dot(np.diag(uk_output_i),np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])))
                    den = np.transpose(np.tile(np.dot(uk_output_i,np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])),(meta['use_dd']['len_col'],1)))                                      
                    weightedEXIOUKconci[eind_slc[m],uind_slc[n]] = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
        
        exioNRG[yrs[a]] = df(np.dot(np.sum(fullexioNRG[yrs[a]],1),weightedEXIOUKconci), index = meta['v']['col'])  
      
    return exioNRG

def make_nrg_2023(uk_energy_filepath,exioNRG,S,yrs,meta): # used in ukmrio_main_2023
    
    ukd = pd.read_excel(os.path.join(uk_energy_filepath,'UKenergy2023.xlsx'), sheet_name = 'data', skiprows=40, usecols='B:AJ')
    ukc = pd.read_excel(os.path.join(uk_energy_filepath,'UKenergy2023.xlsx'), sheet_name = 'conc', index_col=0,header = [0,1])
    
    erg = {}
     
    for a in range(0,np.size(yrs)):
        temp = np.zeros(shape=(meta['v']['len']))
         
        weightedNRGconc = np.zeros(shape=(len(ukc),meta['use']['len_col']))
         
        uk_output_i = df.sum(S[yrs[a]].iloc[meta['v_d']['rng'],:], 1)
                                    
        num = np.transpose(np.dot(np.diag(uk_output_i),np.transpose(ukc)))
        den = np.transpose(np.tile(np.dot(uk_output_i,np.transpose(ukc)),(meta['use_dd']['len_col'],1)))                                      
        weightedNRGconc = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
              
        temp[meta['v_d']['rng']] = np.dot(ukd.iloc[a,:],weightedNRGconc)
        temp[meta['v_i']['rng']] = np.sum(exioNRG[yrs[a]].iloc[meta['v_i']['rng']],1)
        
        erg[yrs[a]] = df(temp, index = meta['use']['col'])
             
    return erg