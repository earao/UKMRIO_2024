# -*- coding: utf-8 -*-
"""
Spyder Editor

This function is used with ukmio_main.py to build the materials extension vector from WU data

@author: earao
"""
import numpy as np
import pandas as pd
import os
df = pd.DataFrame

############################
# used in ukmrio_main_2024 #
############################

def make_exio_stressor_382(mat_filepath,yrs): # used in ukmrio_main_2024
    
    tempdata = pd.read_excel(os.path.join(mat_filepath, 'mat_data.xlsx'), '90-21', usecols = 'A:AL')
    conc_mat = pd.read_excel(os.path.join(mat_filepath, 'mat_data.xlsx'), 'mat_EXIO',index_col=0)
    conc_bio = pd.read_excel(os.path.join(mat_filepath, 'mat_data.xlsx'), 'bio_EXIO',index_col=0)
    conc_ore = pd.read_excel(os.path.join(mat_filepath, 'mat_data.xlsx'), 'ore_EXIO',index_col=0)
    conc_nmm = pd.read_excel(os.path.join(mat_filepath, 'mat_data.xlsx'), 'nmm_EXIO',index_col=0)
    conc_ffl = pd.read_excel(os.path.join(mat_filepath, 'mat_data.xlsx'), 'ffl_EXIO',index_col=0)   
    conc_c = pd.read_excel(os.path.join(mat_filepath, 'mat_data.xlsx'), 'countries_2_EXIO', index_col=0)
    
    data = np.zeros(shape=[62*219,32])
    
    c_list = conc_c.index
    m_list = conc_mat.index
      
    for i in range(0,len(tempdata)):
        for c,country in enumerate(c_list):
            for m,material in enumerate(m_list):
                if tempdata.iloc[i,0]==country:
                    if tempdata.iloc[i,2] == material:
                        data[c*62+m,:]=tempdata.iloc[i,6:].values             
     
    bigconc_mat = np.zeros(shape=[62*219,163*49])
    bigconc_bio = np.zeros(shape=[62*219,163*49])
    bigconc_ore = np.zeros(shape=[62*219,163*49])
    bigconc_nmm = np.zeros(shape=[62*219,163*49])
    bigconc_ffl = np.zeros(shape=[62*219,163*49])
    
    for c in range(0,219):
        for d in range(0,49):
            if conc_c.iloc[c,d] == 1:
                bigconc_mat[c*62:(c+1)*62,d*163:(d+1)*163]=conc_mat.values
                bigconc_bio[c*62:(c+1)*62,d*163:(d+1)*163]=conc_bio.values
                bigconc_ore[c*62:(c+1)*62,d*163:(d+1)*163]=conc_ore.values
                bigconc_nmm[c*62:(c+1)*62,d*163:(d+1)*163]=conc_nmm.values
                bigconc_ffl[c*62:(c+1)*62,d*163:(d+1)*163]=conc_ffl.values

    exioMAT = {} 
    
    for y, yr in enumerate(yrs):
        print(yr)
        temp = data[:,y]
        exioMAT[yr] =  df(np.dot(temp,bigconc_mat), columns = {'mat'})  
        exioMAT[yr]['bio'] = df(np.dot(temp,bigconc_bio))  
        exioMAT[yr]['ore'] =  df(np.dot(temp,bigconc_ore))  
        exioMAT[yr]['nmm'] =  df(np.dot(temp,bigconc_nmm))  
        exioMAT[yr]['ffl'] =  df(np.dot(temp,bigconc_ffl))  
        
    return exioMAT

def make_UK_exioMAT(exioMAT,use,yrs,meta,c_conc,i_conc): # used in ukmrio_main_2024
    
    exioMAT2 = {}
          
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
        temp = exioMAT[yrs[a]]/1000
        exioMAT2[yrs[a]] = df(np.dot(np.transpose(weightedEXIOUKconci),temp), index = meta['v']['col'], columns = exioMAT[yrs[a]].columns)  
      
    return exioMAT2

def make_uk_stressor(ons_filepath,yrs): # used in ukmrio_main_2024
 
    data = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'Domestic extraction', skiprows = 4, usecols = 'B:AG')
    conc_crop = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'crop',index_col=0)
    conc_crpr = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'crop_res',index_col=0)
    conc_wood = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'wood',index_col=0)
    conc_wild = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'wild_catch',index_col=0)
    conc_iron = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'iron',index_col=0)
    conc_nfer = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'non_fer',index_col=0)
    conc_lmgp = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'lime_gyp',index_col=0)
    conc_clay = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'clay',index_col=0)
    conc_sand = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'sand',index_col=0)
    conc_othr = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'other',index_col=0)
    conc_coal = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'coal',index_col=0)
    conc__oil = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'oil',index_col=0)
    conc_ngas = pd.read_excel(os.path.join(ons_filepath, 'ONS environmental accounts/2024/materialflows2023final.xlsx'), 'gas',index_col=0)

    data.fillna(value=0, inplace=True)
     
    uk_MAT_sectors = {} 
    count = 0
    
    for yr in yrs:
        print(yr)
        temp = data.iloc[:,count].values
        uk_MAT_sectors[yr] =  df(np.dot(temp,conc_crop), columns = {'crop'})  
        uk_MAT_sectors[yr]['crop_residue'] = df(np.dot(temp,conc_crpr))  
        uk_MAT_sectors[yr]['wood'] =  df(np.dot(temp,conc_wood))
        uk_MAT_sectors[yr]['wild_catch'] =  df(np.dot(temp,conc_wild))   
        uk_MAT_sectors[yr]['iron'] =  df(np.dot(temp,conc_iron))
        uk_MAT_sectors[yr]['non_ferrous'] =  df(np.dot(temp,conc_nfer)) 
        uk_MAT_sectors[yr]['limestone_gypsum'] =  df(np.dot(temp,conc_lmgp))
        uk_MAT_sectors[yr]['clay'] =  df(np.dot(temp,conc_clay))
        uk_MAT_sectors[yr]['sand'] =  df(np.dot(temp,conc_sand))
        uk_MAT_sectors[yr]['other'] =  df(np.dot(temp,conc_othr))  
        uk_MAT_sectors[yr]['coal'] =  df(np.dot(temp,conc_coal))
        uk_MAT_sectors[yr]['oil'] =  df(np.dot(temp,conc__oil))
        uk_MAT_sectors[yr]['natural_gas'] =  df(np.dot(temp,conc_ngas)) 
        count = count+1
        
    return uk_MAT_sectors

def make_mat(uk_MAT_sectors,exioMAT,S,yrs,meta): # used in ukmrio_main_2024
    
    mat = {}
    bio = {}
    ore = {}
    nmm = {}
    ffl = {}
    
    for a in range(0,np.size(yrs)):
        temp_mat = np.zeros(shape=(meta['v']['len']))       
        temp_mat[meta['v_d']['rng']] = np.sum(uk_MAT_sectors[yrs[a]].values,1)
        temp_mat[meta['v_i']['rng']] = exioMAT[yrs[a]].iloc[meta['v_i']['rng'],0]
        mat[yrs[a]] = df(temp_mat, index = meta['use']['col'])
        
        temp_bio = np.zeros(shape=(meta['v']['len']))       
        temp_bio[meta['v_d']['rng']] = np.sum(uk_MAT_sectors[yrs[a]].iloc[:,0:4].values,1)
        temp_bio[meta['v_i']['rng']] = exioMAT[yrs[a]].iloc[meta['v_i']['rng'],1]
        bio[yrs[a]] = df(temp_bio, index = meta['use']['col'])
        
        temp_ore = np.zeros(shape=(meta['v']['len']))       
        temp_ore[meta['v_d']['rng']] = np.sum(uk_MAT_sectors[yrs[a]].iloc[:,4:6].values,1)
        temp_ore[meta['v_i']['rng']] =exioMAT[yrs[a]].iloc[meta['v_i']['rng'],2]
        ore[yrs[a]] = df(temp_ore, index = meta['use']['col'])
        
        temp_nmm = np.zeros(shape=(meta['v']['len']))       
        temp_nmm[meta['v_d']['rng']] = np.sum(uk_MAT_sectors[yrs[a]].iloc[:,6:10].values,1)
        temp_nmm[meta['v_i']['rng']] =exioMAT[yrs[a]].iloc[meta['v_i']['rng'],3]
        nmm[yrs[a]] = df(temp_nmm, index = meta['use']['col'])
        
        temp_ffl = np.zeros(shape=(meta['v']['len']))       
        temp_ffl[meta['v_d']['rng']] = np.sum(uk_MAT_sectors[yrs[a]].iloc[:,10:13].values,1)
        temp_ffl[meta['v_i']['rng']] =exioMAT[yrs[a]].iloc[meta['v_i']['rng'],4]
        ffl[yrs[a]] = df(temp_ffl, index = meta['use']['col'])
    
    return (mat,bio,ore,nmm,ffl)
