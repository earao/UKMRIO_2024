# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 11:27:33 2019
@author: earao
"""
import io_functions as io
import numpy as np
import pandas as pd
df = pd.DataFrame
import os

###########################
# used in defra_main_2024 #
###########################

def makeukresults(wd, S,U,Y,newY,coicop_exp_tot2,meta,stressor,direct,indicator,allyears,years):

    defra_foot = {}
    tempregagg = np.zeros((meta['fd']['len_idx'],meta['reg']['len'])) 
    tempsecagg = np.zeros((meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    for r in range(0,meta['reg']['len']):
        tempregagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],r] = np.transpose(np.ones((meta['sup_dd']['len_idx'],1)))
        tempsecagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],:] = np.identity(meta['sup_dd']['len_idx'])
    regagg = np.zeros((2*meta['fd']['len_idx'],meta['reg']['len'])) 
    secagg = np.zeros((2*meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    prdagg = np.zeros((2*meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    regagg[0:meta['fd']['len_idx'],:] = tempregagg
    secagg[0:meta['fd']['len_idx'],:] = tempsecagg
    prdagg[meta['fd']['len_idx']:,:] = tempsecagg
    coicop = np.zeros((np.size(newY[2012],1)-1,len(years)))
    coicop2 = np.zeros((np.size(newY[2012],1)-1,len(years)))
    if indicator == 'ghg':
        drct = np.zeros((4,len(allyears)))
    elif indicator == 'co2':
        drct = np.zeros((4,len(allyears)))
    elif indicator == 'nrg':
        drct = np.zeros((6,len(allyears)))
    else:
        drct = np.zeros((np.size(direct,0),len(allyears)))
    sic = np.zeros((meta['sup_dd']['len_idx'],len(allyears)))
    
    energy_filepath = wd + 'data/processed/uk energy/'
    ghg_props = pd.read_excel(os.path.join(energy_filepath, 'UKenergy2024.xlsx'), sheet_name='ghg_props', header = 0, index_col=0)
    co2_props = pd.read_excel(os.path.join(energy_filepath, 'UKenergy2024.xlsx'), sheet_name='co2_props', header = 0, index_col=0)
        
    for i,yr in enumerate(allyears):
        print(yr)
        Z = io.make_Z_from_S_U(S[yr],U[yr])
        bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
        bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
        x = io.make_x(Z,bigY)
        L = io.make_L(Z,x)
        bigstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
        bigstressor[0:np.size(Y[yr],0),:] = stressor[yr]
        e = np.sum(bigstressor,1)/x 
        eL = np.dot(np.diag(e),L)
        reg = np.zeros((meta['reg']['len'],40))
        prd = np.zeros((meta['sup_dd']['len_idx'],40))
        ccc = np.zeros((len(U[yr]),meta['sup_dd']['len_idx']))
        ygg = np.dot(np.sum(bigY[:,0:np.shape(bigY)[1]-1],1),prdagg)
        
        for ysec in range (0,40):
            reg[:,ysec] = np.dot(np.dot(eL,bigY[:,ysec]),regagg)
            prd[:,ysec] = np.dot(np.dot(np.sum(eL,0),np.diag(bigY[:,ysec])),prdagg)
        ccc = np.dot(np.transpose(secagg),np.dot(np.dot(eL,np.diag(np.sum(bigY[:,0:-1],1))),prdagg))
        sic[:,i] = np.sum(prd,1)/ygg
       
        if indicator == 'ghg':
            drct[0,i] = direct.loc['Consumer expenditure - not travel',yr]*ghg_props.loc[yr,'Gas prop']
            drct[1,i] = direct.loc['Consumer expenditure - not travel',yr]*ghg_props.loc[yr,'Liquid fuel prop']
            drct[2,i] = direct.loc['Consumer expenditure - not travel',yr]*ghg_props.loc[yr,'Solid fuel prop']
            drct[3,i] = direct.loc['Consumer expenditure - travel',yr]
        if indicator == 'co2':
           drct[0,i] = direct.loc['Consumer expenditure - not travel',yr]*co2_props.loc[yr,'Gas prop']
           drct[1,i] = direct.loc['Consumer expenditure - not travel',yr]*co2_props.loc[yr,'Liquid fuel prop']
           drct[2,i] = direct.loc['Consumer expenditure - not travel',yr]*co2_props.loc[yr,'Solid fuel prop']
           drct[3,i] = direct.loc['Consumer expenditure - travel',yr]
        if indicator == 'nrg':
            drct[:,i] = direct.loc[:,yr]
        if indicator == 'blc':
            drct[:,i] = direct.loc[yr,:]
        if indicator == 'blw':
            drct[:,i] = direct.loc[yr,:]
        
        defra_foot[str(yr)+'_sic'] = df(ccc, index = meta['sectors']['ind'], columns = meta['sectors']['prd'])
        defra_foot[str(yr)+'_reg'] = df(reg, index = meta['reg']['idx'], columns = Y[yr].columns[0:np.size(Y[yr],1)-1])
       
    if indicator == 'ghg':
        defra_foot['direct'] = df(drct, index = ['Consumer expenditure gas','Consumer expenditure liquid fuel', 'Consumer expenditure solid fuel', 'Consumer expenditure travel'], columns = allyears)
    elif indicator == 'co2':
        defra_foot['direct'] = df(drct, index = ['Consumer expenditure gas','Consumer expenditure liquid fuel', 'Consumer expenditure solid fuel', 'Consumer expenditure travel'], columns = allyears)
    elif indicator == 'nrg':
        defra_foot['direct'] = df(drct, index = ['Consumer expenditure electricity','Consumer expenditure gas','Consumer expenditure liquid fuel', 'Consumer expenditure solid fuel', 'Consumer expenditure heat fuel', 'Consumer expenditure travel'], columns = allyears)
    else:
        defra_foot['direct'] = df(drct, index = direct.index, columns = allyears)
        
    for i,yr in enumerate(years): 
        print(yr)
        Z = io.make_Z_from_S_U(S[yr],U[yr])
        bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
        bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
        newbigY = np.zeros(shape = [np.size(newY[yr],0)*2,np.size(newY[yr],1)])
        newbigY[np.size(Y[yr],0):np.size(newY[yr],0)*2,0:np.size(newY[yr],1)] = newY[yr]
        x = io.make_x(Z,bigY)
        L = io.make_L(Z,x)
        bigstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
        bigstressor[0:np.size(Y[yr],0),:] = stressor[yr]
        e = np.sum(bigstressor,1)/x 
        eL_vector = np.dot(e,L)
      
        for ysec2 in range (0,np.size(newY[yr],1)-1):
            coicop[ysec2,i] = np.dot(eL_vector,newbigY[:,ysec2])
        if indicator == 'ghg':
            coicop[31,i] = coicop[31,i] + defra_foot['direct'].loc['Consumer expenditure gas',yr]
            coicop[32,i] = coicop[32,i] + defra_foot['direct'].loc['Consumer expenditure liquid fuel',yr]
            coicop[33,i] = coicop[33,i] + defra_foot['direct'].loc['Consumer expenditure solid fuel',yr]
            coicop[59,i] = coicop[59,i] + defra_foot['direct'].loc['Consumer expenditure travel',yr]
        if indicator == 'co2':
            coicop[31,i] = coicop[31,i] + defra_foot['direct'].loc['Consumer expenditure gas',yr]
            coicop[32,i] = coicop[32,i] + defra_foot['direct'].loc['Consumer expenditure liquid fuel',yr]
            coicop[33,i] = coicop[33,i] + defra_foot['direct'].loc['Consumer expenditure solid fuel',yr]
            coicop[59,i] = coicop[59,i] + defra_foot['direct'].loc['Consumer expenditure travel',yr]
        if indicator == 'nrg':
            coicop[30,i] = coicop[30,i] + defra_foot['direct'].loc['Consumer expenditure electricity',yr]
            coicop[31,i] = coicop[31,i] + defra_foot['direct'].loc['Consumer expenditure gas',yr]
            coicop[32,i] = coicop[32,i] + defra_foot['direct'].loc['Consumer expenditure liquid fuel',yr]
            coicop[33,i] = coicop[33,i] + defra_foot['direct'].loc['Consumer expenditure solid fuel',yr]
            coicop[34,i] = coicop[34,i] + defra_foot['direct'].loc['Consumer expenditure heat fuel',yr]
            coicop[59,i] = coicop[59,i] + defra_foot['direct'].loc['Consumer expenditure travel',yr]
        if indicator == 'wat':
            coicop[26,i] = coicop[26,i] + direct[yr]
                  
        coicop2[0:105,i] = np.sum(df((coicop[0:105,i]*1000000))/df((coicop_exp_tot2[yr].values*52*1000)),1)
        coicop2[105,i] =  (coicop[105,i]*1000000)/(np.sum(Y[yr].loc[:,('13 Non-profit instns serving households')],0)*1000000)
        coicop2[106,i] =  (coicop[106,i]*1000000)/(np.sum(Y[yr].loc[:,('14 Central_x000D_ government')],0)*1000000)
        coicop2[107,i] =  (coicop[107,i]*1000000)/(np.sum(Y[yr].loc[:,('15 Local_x000D_ Authorities')],0)*1000000)
        coicop2[108,i] =  (coicop[108,i]*1000000)/(np.sum(Y[yr].loc[:,('16 Gross fixed_x000D_ capital_x000D_ formation')],0)*1000000)
        coicop2[109,i] =  (coicop[109,i]*1000000)/(np.sum(Y[yr].loc[:,('17 Valuables')],0)*1000000)
        coicop2[110,i] =  (coicop[110,i]*1000000)/(np.sum(Y[yr].loc[:,('18 Changes in inventories')],0)*1000000)
    
    defra_foot['coicop'] = df(coicop, index = newY[yr].columns[0:np.size(newY[yr],1)-1], columns = years)
    defra_foot['coicop_mult'] = df(coicop2, index = newY[yr].columns[0:np.size(newY[yr],1)-1], columns = years)
    defra_foot['sic_mult'] = df(sic, index =  meta['sectors']['prd'], columns = allyears)
        
    return defra_foot

def printdefradata(region,results_filepath,years,results,indicator):
    
    defrawriter = os.path.join(results_filepath, 'Defra_results_'+region+'.xlsx')
    writer = pd.ExcelWriter(defrawriter)
   
    for item in range(0,len(results)):
        for item2 in results[item]:
            sheetlabel = [indicator[item],item2]
            df(results[item][item2]).to_excel(writer,('_'.join(sheetlabel)))
    print(defrawriter + ' saved')
    writer.save()
    
    return

##########################################
# used in defra_uk_devolved_regions_2023 #
##########################################

# Also used in defra_main_2023 (see defra_main_2023):
# makeukresults2023
# printdefradata

def makeregionresults(S,U,Y,newY,Yregion,meta,stressor,direct,indicator,years,concs_dict,coicop_mult,sic_mult,regpophholdsyr,Region_Name):
    defra_foot = {}
    tempregagg = np.zeros((meta['fd']['len_idx'],meta['reg']['len'])) 
    tempsecagg = np.zeros((meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    for r in range(0,meta['reg']['len']):
        tempregagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],r] = np.transpose(np.ones((meta['sup_dd']['len_idx'],1)))
        tempsecagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],:] = np.identity(meta['sup_dd']['len_idx'])
    regagg = np.zeros((2*meta['fd']['len_idx'],meta['reg']['len'])) 
    secagg = np.zeros((2*meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    prdagg = np.zeros((2*meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    regagg[0:meta['fd']['len_idx'],:] = tempregagg
    secagg[0:meta['fd']['len_idx'],:] = tempsecagg
    prdagg[meta['fd']['len_idx']:,:] = tempsecagg
    coicop = np.zeros((np.size(newY[2012],1)-1,len(years)))
    if indicator == 'ghg':
        drct = np.zeros((4,len(years)))
    elif indicator == 'co2':
        drct = np.zeros((4,len(years)))
    elif indicator == 'nrg':
        drct = np.zeros((6,len(years)))
    else:
        drct = np.zeros((np.size(direct,0),len(years)))
    pop = np.zeros((2,len(years)))
    nipop = [1686000,1693000,1701000,1709000,1721000,1735000,1752000,1770000,1786000,1799000,1810000,1819000,1827000,1835000,1846000,1857000,1867000,1876000,1885000,1896000]
  
    for i, yr in enumerate(years):
        print(yr)
        Z = io.make_Z_from_S_U(S[yr],U[yr])
        bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
        bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
        bigYregion = np.zeros(shape = [np.size(Yregion[yr],0)*2,np.size(Yregion[yr],1)])
        bigYregion[np.size(Yregion[yr],0):np.size(Yregion[yr],0)*2,0:np.size(Yregion[yr],1)] = Yregion[yr]       
        x = io.make_x(Z,bigY)
        L = io.make_L(Z,x)
        bigstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
        bigstressor[0:np.size(Y[yr],0),:] = stressor[yr]
        e = np.sum(bigstressor,1)/x
        eL_vector = np.dot(e,L)
        eL = np.dot(np.diag(e),L)
        reg = np.zeros((meta['reg']['len'],111))
        ind = np.zeros((meta['sup_dd']['len_idx'],111))
        
        for ysec in range (0,np.size(Yregion[yr],1)-1):
            reg[:,ysec] = np.dot(np.dot(eL,bigYregion[:,ysec]),regagg)
            ind[:,ysec] = np.dot(np.dot(eL,bigYregion[:,ysec]),secagg)
            coicop[ysec,i] = np.dot(eL_vector,bigYregion[:,ysec])
        
        ccc = np.dot(eL,np.diag(np.sum(bigYregion[:,0:-1],1)))
        ccc = np.dot(np.transpose(secagg),(np.dot(ccc,prdagg))) 
        
        defra_foot[str(yr)+'_reg'] = df(reg, index = meta['reg']['idx'], columns = Yregion[yr].columns[0:np.size(Yregion[yr],1)-1])
        defra_foot[str(yr)+'_sic'] = df(ccc, index = meta['sectors']['ind'], columns = meta['sectors']['prd'])
       
        if indicator == 'ghg':
            drct[0,i] = direct.loc['Consumer expenditure gas',yr]*np.sum(bigYregion[:,31],0)/np.sum(newY[yr].iloc[:,31],0)
            drct[1,i] = direct.loc['Consumer expenditure liquid fuel',yr]*np.sum(bigYregion[:,32],0)/np.sum(newY[yr].iloc[:,32],0)
            drct[2,i] = direct.loc['Consumer expenditure solid fuel',yr]*np.sum(bigYregion[:,33],0)/np.sum(newY[yr].iloc[:,33],0)
            drct[3,i] = direct.loc['Consumer expenditure travel',yr]*np.sum(bigYregion[:,59],0)/np.sum(newY[yr].iloc[:,59],0)
            
        if indicator == 'co2':
            drct[0,i] = direct.loc['Consumer expenditure gas',yr]*np.sum(bigYregion[:,31],0)/np.sum(newY[yr].iloc[:,31],0)
            drct[1,i] = direct.loc['Consumer expenditure liquid fuel',yr]*np.sum(bigYregion[:,32],0)/np.sum(newY[yr].iloc[:,32],0)
            drct[2,i] = direct.loc['Consumer expenditure solid fuel',yr]*np.sum(bigYregion[:,33],0)/np.sum(newY[yr].iloc[:,33],0)
            drct[3,i] = direct.loc['Consumer expenditure travel',yr]*np.sum(bigYregion[:,59],0)/np.sum(newY[yr].iloc[:,59],0)
            
        if indicator == 'nrg':
           drct[0,i] = direct.loc['Consumer expenditure electricity',yr]*np.sum(bigYregion[:,30],0)/np.sum(newY[yr].iloc[:,30],0)
           drct[1,i] = direct.loc['Consumer expenditure gas',yr]*np.sum(bigYregion[:,31],0)/np.sum(newY[yr].iloc[:,31],0)
           drct[2,i] = direct.loc['Consumer expenditure liquid fuel',yr]*np.sum(bigYregion[:,32],0)/np.sum(newY[yr].iloc[:,32],0)
           drct[3,i] = direct.loc['Consumer expenditure solid fuel',yr]*np.sum(bigYregion[:,33],0)/np.sum(newY[yr].iloc[:,33],0)
           drct[4,i] = direct.loc['Consumer expenditure heat fuel',yr]*np.sum(bigYregion[:,34],0)/np.sum(newY[yr].iloc[:,34],0)
           drct[5,i] = direct.loc['Consumer expenditure travel',yr]*np.sum(bigYregion[:,59],0)/np.sum(newY[yr].iloc[:,59],0)
            
        if indicator == 'blc':
            drct[:,i] = direct.loc[:,yr]*np.sum(bigYregion[:,26],0)/np.sum(newY[yr].iloc[:,26],0)
            
        if indicator == 'blw':
            drct[:,i] = direct.loc[:,yr]*np.sum(bigYregion[:,26],0)/np.sum(newY[yr].iloc[:,26],0)
        
        if Region_Name == 'England':
            pop[1,i] = np.sum(regpophholdsyr[yr]['North East'].loc[:,'pop'])+\
                np.sum(regpophholdsyr[yr]['North West'].loc[:,'pop'])+\
                    np.sum(regpophholdsyr[yr]['Yorkshire and The Humber'].loc[:,'pop'])+\
                        np.sum(regpophholdsyr[yr]['East Midlands'].loc[:,'pop'])+\
                            np.sum(regpophholdsyr[yr]['West Midlands'].loc[:,'pop'])+\
                                np.sum(regpophholdsyr[yr]['East'].loc[:,'pop'])+\
                                    np.sum(regpophholdsyr[yr]['London'].loc[:,'pop'])+\
                                        np.sum(regpophholdsyr[yr]['South East'].loc[:,'pop'])+\
                                            np.sum(regpophholdsyr[yr]['South West'].loc[:,'pop'])
        elif Region_Name == 'Northern Ireland':
            pop[1,i] = nipop[i]
        else:
            pop[1,i] = np.sum(regpophholdsyr[yr][Region_Name].loc[:,'pop'])
            
        pop[0,i] = np.sum(regpophholdsyr[yr]['North East'].loc[:,'pop'])+\
            np.sum(regpophholdsyr[yr]['North West'].loc[:,'pop'])+\
                np.sum(regpophholdsyr[yr]['Yorkshire and The Humber'].loc[:,'pop'])+\
                    np.sum(regpophholdsyr[yr]['East Midlands'].loc[:,'pop'])+\
                        np.sum(regpophholdsyr[yr]['West Midlands'].loc[:,'pop'])+\
                            np.sum(regpophholdsyr[yr]['East'].loc[:,'pop'])+\
                                np.sum(regpophholdsyr[yr]['London'].loc[:,'pop'])+\
                                    np.sum(regpophholdsyr[yr]['South East'].loc[:,'pop'])+\
                                        np.sum(regpophholdsyr[yr]['South West'].loc[:,'pop'])+\
                                            np.sum(regpophholdsyr[yr]['Wales'].loc[:,'pop'])+\
                                                np.sum(regpophholdsyr[yr]['Scotland'].loc[:,'pop'])+nipop[i]

                
    defra_foot['direct'] = df(drct, index = direct.index, columns = years)
    defra_foot['coicop'] = df(coicop, index = Yregion[yr].columns[0:np.size(Yregion[yr],1)-1], columns = years)
    defra_foot['coicop_mult'] = coicop_mult
    defra_foot['sic_mult'] = sic_mult
    defra_foot['population'] = df(pop, index = ['United Kingdom',Region_Name],columns = years)        
       
    return defra_foot

#################################
# used in generations_2023_main #
#################################

# Also used in defra_main_2023 (see defra_main_2023):
# makeukresults2023

#############################################
# used in sub_national_footprints_2023_main #
#############################################

# Also used in defra_main_2023 (see defra_main_2023):
# makeukresults2023
# printdefradata

def printdefradata2(region,results_filepath,years,results,indicator):
    
    defrawriter = os.path.join(results_filepath, 'CBA_results_'+region+'.xlsx')
    writer = pd.ExcelWriter(defrawriter)
   
    for item in results:
        print(item)
        sheetlabel = [indicator,item]
        results[item].to_excel(writer,('_'.join(sheetlabel)))
           
    writer.save()
    
    return

def lacheck2023(defra_region,region_las,years):
    check={}
    temppop =  np.zeros((len(region_las),len(years)))
    for j, yr in enumerate(years):
        tempcheck =  np.zeros((len(region_las),111))
        label = []
        for i, la in enumerate(region_las):
            tempcheck[i,:] = np.sum(region_las[la][str(yr)+'_reg']).values
            tempcheck[i,31] = tempcheck[i,31]+region_las[la]['direct'].loc['Consumer expenditure gas', yr]
            tempcheck[i,32] = tempcheck[i,32]+region_las[la]['direct'].loc['Consumer expenditure liquid fuel', yr]
            tempcheck[i,33] = tempcheck[i,33]+region_las[la]['direct'].loc['Consumer expenditure solid fuel', yr]
            tempcheck[i,59] = tempcheck[i,59]+region_las[la]['direct'].loc['Consumer expenditure travel', yr]      
            label.append(la)
            temppop[i,j] = region_las[la]['population'].loc[la,yr]
        tempcheck = df(tempcheck, columns = region_las[la][str(yr)+'_reg'].columns, index = label)        
        check[yr] = tempcheck       
    check['pop'] = df(temppop, columns = years, index = label)
    return check

def regioncheck2023(defra_ghg_reg,years):
    
    check={}
    temppop =  np.zeros((len(defra_ghg_reg),len(years)))
    for y, yr in enumerate(years):
        tempcheck =  np.zeros((len(defra_ghg_reg),111))
        for r, reg in enumerate(defra_ghg_reg):
            tempcheck[r,:] = np.sum(defra_ghg_reg[reg][str(yr)+'_reg']).values
        tempcheck = df(tempcheck, columns = defra_ghg_reg[reg][str(yr)+'_reg'].columns, index = defra_ghg_reg.keys())
        
        for r, reg in enumerate(defra_ghg_reg):
            tempcheck.iloc[r,31] = tempcheck.iloc[r,31]+defra_ghg_reg[reg]['direct'].loc['Consumer expenditure gas', yr]
            tempcheck.iloc[r,32] = tempcheck.iloc[r,32]+defra_ghg_reg[reg]['direct'].loc['Consumer expenditure liquid fuel', yr]
            tempcheck.iloc[r,33] = tempcheck.iloc[r,33]+defra_ghg_reg[reg]['direct'].loc['Consumer expenditure solid fuel', yr]
            tempcheck.iloc[r,59] = tempcheck.iloc[r,59]+defra_ghg_reg[reg]['direct'].loc['Consumer expenditure travel', yr]
            
            temppop[r,y] = defra_ghg_reg[reg]['population'].iloc[1,y]
        
          
        check[yr] = tempcheck.iloc[1:,:]
    check['pop'] = df(temppop[1:,:], columns = years, index = check[yr].index)
    return check

##################
# Not sorted yet #
##################

def regioncheck2022(defra_ghg_north_east,defra_ghg_north_west,defra_ghg_yorkshire_and_the_humber,defra_ghg_east_midlands,defra_ghg_west_midlands,defra_ghg_east,defra_ghg_london,defra_ghg_south_east,defra_ghg_south_west,defra_ghg_scotland,defra_ghg_wales,defra_ghg_northern_ireland,years):
    check={}
    temppop =  np.zeros((12,len(years)))
    for i, yr in enumerate(years):
        tempcheck =  np.zeros((12,115))
        tempcheck[0,:] = np.sum(defra_ghg_north_east[str(yr)+'_reg']).values
        tempcheck[1,:] = np.sum(defra_ghg_north_west[str(yr)+'_reg']).values
        tempcheck[2,:] = np.sum(defra_ghg_yorkshire_and_the_humber[str(yr)+'_reg']).values
        tempcheck[3,:] = np.sum(defra_ghg_east_midlands[str(yr)+'_reg']).values
        tempcheck[4,:] = np.sum(defra_ghg_west_midlands[str(yr)+'_reg']).values
        tempcheck[5,:] = np.sum(defra_ghg_east[str(yr)+'_reg']).values
        tempcheck[6,:] = np.sum(defra_ghg_london[str(yr)+'_reg']).values
        tempcheck[7,:] = np.sum(defra_ghg_south_east[str(yr)+'_reg']).values
        tempcheck[8,:] = np.sum(defra_ghg_south_west[str(yr)+'_reg']).values
        tempcheck[9,:] = np.sum(defra_ghg_scotland[str(yr)+'_reg']).values
        tempcheck[10,:] = np.sum(defra_ghg_wales[str(yr)+'_reg']).values
        tempcheck[11,:] = np.sum(defra_ghg_northern_ireland[str(yr)+'_reg']).values
        tempcheck = df(tempcheck, columns = defra_ghg_north_east[str(yr)+'_reg'].columns, index = ['North East','North West','Yorkshire and The Humber','East Midlands','West Midlands','East','London','South East','South West','Scotland','Wales','Northern Ireland'])
        
        tempcheck.iloc[0,31] = tempcheck.iloc[0,31]+defra_ghg_north_east['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[1,31] = tempcheck.iloc[1,31]+defra_ghg_north_west['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[2,31] = tempcheck.iloc[2,31]+defra_ghg_yorkshire_and_the_humber['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[3,31] = tempcheck.iloc[3,31]+defra_ghg_east_midlands['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[4,31] = tempcheck.iloc[4,31]+defra_ghg_west_midlands['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[5,31] = tempcheck.iloc[5,31]+defra_ghg_east['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[6,31] = tempcheck.iloc[6,31]+defra_ghg_london['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[7,31] = tempcheck.iloc[7,31]+defra_ghg_south_east['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[8,31] = tempcheck.iloc[8,31]+defra_ghg_south_west['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[9,31] = tempcheck.iloc[9,31]+defra_ghg_scotland['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[10,31] = tempcheck.iloc[10,31]+defra_ghg_wales['direct'].loc['Consumer expenditure gas', yr]
        tempcheck.iloc[11,31] = tempcheck.iloc[11,31]+defra_ghg_northern_ireland['direct'].loc['Consumer expenditure gas', yr]
        
        tempcheck.iloc[0,32] = tempcheck.iloc[0,32]+defra_ghg_north_east['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[1,32] = tempcheck.iloc[1,32]+defra_ghg_north_west['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[2,32] = tempcheck.iloc[2,32]+defra_ghg_yorkshire_and_the_humber['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[3,32] = tempcheck.iloc[3,32]+defra_ghg_east_midlands['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[4,32] = tempcheck.iloc[4,32]+defra_ghg_west_midlands['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[5,32] = tempcheck.iloc[5,32]+defra_ghg_east['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[6,32] = tempcheck.iloc[6,32]+defra_ghg_london['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[7,32] = tempcheck.iloc[7,32]+defra_ghg_south_east['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[8,32] = tempcheck.iloc[8,32]+defra_ghg_south_west['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[9,32] = tempcheck.iloc[9,32]+defra_ghg_scotland['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[10,32] = tempcheck.iloc[10,32]+defra_ghg_wales['direct'].loc['Consumer expenditure liquid fuel', yr]
        tempcheck.iloc[11,32] = tempcheck.iloc[11,32]+defra_ghg_northern_ireland['direct'].loc['Consumer expenditure liquid fuel', yr]
    
        tempcheck.iloc[0,33] = tempcheck.iloc[0,33]+defra_ghg_north_east['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[1,33] = tempcheck.iloc[1,33]+defra_ghg_north_west['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[2,33] = tempcheck.iloc[2,33]+defra_ghg_yorkshire_and_the_humber['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[3,33] = tempcheck.iloc[3,33]+defra_ghg_east_midlands['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[4,33] = tempcheck.iloc[4,33]+defra_ghg_west_midlands['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[5,33] = tempcheck.iloc[5,33]+defra_ghg_east['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[6,33] = tempcheck.iloc[6,33]+defra_ghg_london['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[7,33] = tempcheck.iloc[7,33]+defra_ghg_south_east['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[8,33] = tempcheck.iloc[8,33]+defra_ghg_south_west['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[9,33] = tempcheck.iloc[9,33]+defra_ghg_scotland['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[10,33] = tempcheck.iloc[10,33]+defra_ghg_wales['direct'].loc['Consumer expenditure solid fuel', yr]
        tempcheck.iloc[11,33] = tempcheck.iloc[11,33]+defra_ghg_northern_ireland['direct'].loc['Consumer expenditure solid fuel', yr]

        tempcheck.iloc[0,59] = tempcheck.iloc[0,59]+defra_ghg_north_east['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[1,59] = tempcheck.iloc[1,59]+defra_ghg_north_west['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[2,59] = tempcheck.iloc[2,59]+defra_ghg_yorkshire_and_the_humber['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[3,59] = tempcheck.iloc[3,59]+defra_ghg_east_midlands['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[4,59] = tempcheck.iloc[4,59]+defra_ghg_west_midlands['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[5,59] = tempcheck.iloc[5,59]+defra_ghg_east['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[6,59] = tempcheck.iloc[6,59]+defra_ghg_london['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[7,59] = tempcheck.iloc[7,59]+defra_ghg_south_east['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[8,59] = tempcheck.iloc[8,59]+defra_ghg_south_west['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[9,59] = tempcheck.iloc[9,59]+defra_ghg_scotland['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[10,59] = tempcheck.iloc[10,59]+defra_ghg_wales['direct'].loc['Consumer expenditure travel', yr]
        tempcheck.iloc[11,59] = tempcheck.iloc[11,59]+defra_ghg_northern_ireland['direct'].loc['Consumer expenditure travel', yr]
        
        temppop[0,i] = defra_ghg_north_east['population'].iloc[1,i]
        temppop[1,i] = defra_ghg_north_west['population'].iloc[1,i]
        temppop[2,i] = defra_ghg_yorkshire_and_the_humber['population'].iloc[1,i]
        temppop[3,i] = defra_ghg_east_midlands['population'].iloc[1,i]
        temppop[4,i] = defra_ghg_west_midlands['population'].iloc[1,i]
        temppop[5,i] = defra_ghg_east['population'].iloc[1,i]
        temppop[6,i] = defra_ghg_london['population'].iloc[1,i]
        temppop[7,i] = defra_ghg_south_east['population'].iloc[1,i]
        temppop[8,i] = defra_ghg_south_west['population'].iloc[1,i]
        temppop[9,i] = defra_ghg_scotland['population'].iloc[1,i]
        temppop[10,i] = defra_ghg_wales['population'].iloc[1,i]
        temppop[11,i] = defra_ghg_northern_ireland['population'].iloc[1,i]
        
        label = ['North East','North West','Yorkshire and the Humber','East Midlands','West Midlands','East','London','South East','South West','Scotland','Wales','Northern Ireland']
        check[yr] = tempcheck
        check['pop'] = df(temppop, columns = years, index = label)
    return check



def makelaresults2022(defra_region_result,reglaspropyr,regpophholdsyr,LAcodesnames,regoacsyr,oacyrmeta,oacyrspends,oac01_names,oac11_names,Region_Name,years,concs_dict):
    la_foot = {}
    la_names = reglaspropyr[2019][Region_Name].index
    for la in la_names:
        coicop = np.zeros((111,len(years)))
        reg = np.zeros((15,111))
        drct = np.zeros((4,len(years)))
        pop = np.zeros((3,len(years)))
        la_foot_temp = {}
        pop[0] = defra_region_result['population'].iloc[0,:]
        for n, yr in enumerate(years):
            coicopprop = reglaspropyr[yr][Region_Name].loc[la].values  
            reg = defra_region_result[str(yr)+'_reg']*coicopprop
            la_foot_temp[str(yr)+'_reg'] = df(reg, columns = defra_region_result[str(yr)+'_reg'].columns, index = defra_region_result[str(yr)+'_reg'].index)
            drct[0,n] = defra_region_result['direct'].loc['Consumer expenditure gas',yr]*reglaspropyr[yr][Region_Name].loc[la,[('4.5.2 Gas')]]
            drct[1,n] = defra_region_result['direct'].loc['Consumer expenditure liquid fuel',yr]*reglaspropyr[yr][Region_Name].loc[la,[('4.5.3 Liquid fuels')]]
            drct[2,n] = defra_region_result['direct'].loc['Consumer expenditure solid fuel',yr]*reglaspropyr[yr][Region_Name].loc[la,[('4.5.4 Solid fuels')]]
            drct[3,n] = defra_region_result['direct'].loc['Consumer expenditure travel',yr]*reglaspropyr[yr][Region_Name].loc[la,[('7.2.2 Fuels and lubricants for personal transport equipment')]]
            
            pop[1,n] = np.sum(regpophholdsyr[yr][Region_Name].loc[:,'pop'])
            pop[2,n] = regpophholdsyr[yr][Region_Name].loc[la,'pop']
            oac=df(regoacsyr[yr][Region_Name].loc[la,'pop'])
            if yr < 2014:
                for i, c in enumerate(regoacsyr[yr][Region_Name].loc[la].index):
                    oac.loc[c,'OAC name'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                    if oacyrmeta[yr][Region_Name].iloc[i,2] > 19:
                        oac.loc[c,'OAC code used'] = c
                        oac.loc[c,'OAC name used'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                    elif oacyrmeta[yr][Region_Name].iloc[i,1] > 19:
                        oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].index[i][1]
                        oac.loc[c,'OAC name used'] = oac01_names['Group'].loc[oacyrmeta[yr][Region_Name].index[i][1], 'Group name']
                    elif oacyrmeta[yr][Region_Name].iloc[i,0] > 19:
                        oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].index[i][0]
                        oac.loc[c,'OAC name used'] = oac01_names['Supergroup'].loc[int( oacyrmeta[yr][Region_Name].index[i][0]), 'Supergroup name']
                    else:
                        oac.loc[c,'OAC code used'] = Region_Name
                        oac.loc[c,'OAC name used'] = [Region_Name +' average']                            
            else:
                  for i, c in enumerate(regoacsyr[yr][Region_Name].loc[la].index):
                         oac.loc[c,'OAC name'] = oac11_names['Subgroup'].loc[c, 'Subgroup name']
                         if oacyrmeta[yr][Region_Name].iloc[i,2] > 19:
                             oac.loc[c,'OAC code used'] = c
                             oac.loc[c,'OAC name used'] = oac11_names['Subgroup'].loc[c, 'Subgroup name']
                         elif oacyrmeta[yr][Region_Name].iloc[i,1] > 19:
                             oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].index[i][1]
                             oac.loc[c,'OAC name used'] = oac11_names['Group'].loc[oacyrmeta[yr][Region_Name].index[i][1], 'Group name']
                         elif oacyrmeta[yr][Region_Name].iloc[i,0] > 19:
                             oac.loc[c,'OAC code used'] =oacyrmeta[yr][Region_Name].index[i][0]
                             oac.loc[c,'OAC name used'] = oac11_names['Supergroup'].loc[int( oacyrmeta[yr][Region_Name].index[i][0]), 'Supergroup name']
                         else:
                             oac.loc[c,'OAC code used'] = Region_Name 
                             oac.loc[c,'OAC name used'] = [Region_Name +' average']   
              
        la_foot_temp['coicop'] = df(coicop, index = defra_region_result['coicop'].iloc[0:115].index, columns = years)
        la_foot_temp['direct'] = df(drct, index = defra_region_result['direct'].index, columns = years)
        la_foot_temp['population'] = df(pop, index = ['United Kingdom',Region_Name,LAcodesnames.loc[la,'la_name']], columns = years)        
        
        la_foot[LAcodesnames.loc[la,'la_name']] = la_foot_temp
        
    return la_foot





def makelaresults2023(defra_region_result,reglaspropyr,regpophholdsyr,LAcodesnames,regoacsyr,oacyrmeta,oacyrspends,oac01_names,oac11_names,Region_Name,years,concs_dict):
    la_foot = {}
    oac_lookup = {}
    la_names = reglaspropyr[2020][Region_Name].index
    for la in la_names:
        coicop = np.zeros((111,len(years)))
        reg = np.zeros((15,111))
        drct = np.zeros((4,len(years)))
        pop = np.zeros((3,len(years)))
        la_foot_temp = {}
        oac_lookup_temp = {}
        pop[0] = defra_region_result['population'].iloc[0,:]
        for n, yr in enumerate(years):
            coicopprop = reglaspropyr[yr][Region_Name].loc[la].values  
            reg = defra_region_result[str(yr)+'_reg']*coicopprop
            coicop[:,n] = defra_region_result['coicop'].loc[:,yr]*coicopprop
            la_foot_temp[str(yr)+'_reg'] = df(reg, columns = defra_region_result[str(yr)+'_reg'].columns, index = defra_region_result[str(yr)+'_reg'].index)
            drct[0,n] = defra_region_result['direct'].loc['Consumer expenditure gas',yr]*reglaspropyr[yr][Region_Name].loc[la,[('4.5.2 Gas')]]
            drct[1,n] = defra_region_result['direct'].loc['Consumer expenditure liquid fuel',yr]*reglaspropyr[yr][Region_Name].loc[la,[('4.5.3 Liquid fuels')]]
            drct[2,n] = defra_region_result['direct'].loc['Consumer expenditure solid fuel',yr]*reglaspropyr[yr][Region_Name].loc[la,[('4.5.4 Solid fuels')]]
            drct[3,n] = defra_region_result['direct'].loc['Consumer expenditure travel',yr]*reglaspropyr[yr][Region_Name].loc[la,[('7.2.2 Fuels and lubricants for personal transport equipment')]]
            
            pop[1,n] = np.sum(regpophholdsyr[yr][Region_Name].loc[:,'pop'])
            pop[2,n] = regpophholdsyr[yr][Region_Name].loc[la,'pop']
            oac=df(regoacsyr[yr][Region_Name].loc[la,'pop'])
            if yr < 2014:
                for i, c in enumerate(regoacsyr[yr][Region_Name].loc[la].index):
                    oac.loc[c,'OAC name'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                    if oacyrmeta[yr][Region_Name].iloc[i,2] > 19:
                        oac.loc[c,'OAC code used'] = c
                        oac.loc[c,'OAC name used'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                    elif oacyrmeta[yr][Region_Name].iloc[i,1] > 19:
                        oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].index[i][1]
                        oac.loc[c,'OAC name used'] = oac01_names['Group'].loc[oacyrmeta[yr][Region_Name].index[i][1], 'Group name']
                    elif oacyrmeta[yr][Region_Name].iloc[i,0] > 19:
                        oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].index[i][0]
                        oac.loc[c,'OAC name used'] = oac01_names['Supergroup'].loc[int( oacyrmeta[yr][Region_Name].index[i][0]), 'Supergroup name']
                    else:
                        oac.loc[c,'OAC code used'] = Region_Name
                        oac.loc[c,'OAC name used'] = [Region_Name +' average']                                
            else:
                 for i, c in enumerate(regoacsyr[yr][Region_Name].loc[la].index):
                        oac.loc[c,'OAC name'] = oac11_names['Subgroup'].loc[c, 'Subgroup name']
                        if oacyrmeta[yr][Region_Name].iloc[i,2] > 19:
                            oac.loc[c,'OAC code used'] = c
                            oac.loc[c,'OAC name used'] = oac11_names['Subgroup'].loc[c, 'Subgroup name']
                        elif oacyrmeta[yr][Region_Name].iloc[i,1] > 19:
                            oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].index[i][1]
                            oac.loc[c,'OAC name used'] = oac11_names['Group'].loc[oacyrmeta[yr][Region_Name].index[i][1], 'Group name']
                        elif oacyrmeta[yr][Region_Name].iloc[i,0] > 19:
                            oac.loc[c,'OAC code used'] =oacyrmeta[yr][Region_Name].index[i][0]
                            oac.loc[c,'OAC name used'] = oac11_names['Supergroup'].loc[int( oacyrmeta[yr][Region_Name].index[i][0]), 'Supergroup name']
                        else:
                            oac.loc[c,'OAC code used'] = Region_Name 
                            oac.loc[c,'OAC name used'] = [Region_Name +' average']  
            oac_lookup_temp[str(yr)+'oac'] = df(oac)
            
        la_foot_temp['coicop'] = df(coicop, index = defra_region_result['coicop'].iloc[0:111].index, columns = years)
        la_foot_temp['direct'] = df(drct, index = defra_region_result['direct'].index, columns = years)
        la_foot_temp['population'] = df(pop, index = ['United Kingdom',Region_Name,LAcodesnames.loc[la,'la_name']], columns = years)        
        
        
        
        la_foot[LAcodesnames.loc[la,'la_name']] = la_foot_temp
        oac_lookup[LAcodesnames.loc[la,'la_name']] = oac_lookup_temp
        
    return (la_foot,oac_lookup)






         

def makeukresultscc(S,U,Y,newY,coicop_exp_tot,meta,stressor,direct,indicator,allyears,years):
    defra_foot = {}
    tempregagg = np.zeros((meta['fd']['len_idx'],meta['reg']['len'])) 
    tempsecagg = np.zeros((meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    for r in range(0,meta['reg']['len']):
        tempregagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],r] = np.transpose(np.ones((meta['sup_dd']['len_idx'],1)))
        tempsecagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],:] = np.identity(meta['sup_dd']['len_idx'])
    regagg = np.zeros((2*meta['fd']['len_idx'],meta['reg']['len'])) 
    secagg = np.zeros((2*meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    regagg[0:meta['fd']['len_idx'],:] = tempregagg
    secagg[0:meta['fd']['len_idx'],:] = tempsecagg
    coicop = np.zeros((np.size(newY[2012],1)-1,len(years)))
    coicop2 = np.zeros((np.size(newY[2012],1)-1,len(years)))
    drct = np.zeros((np.size(direct,0),len(allyears)))
    
    for i,yr in enumerate(allyears):
        Z = io.make_Z_from_S_U(S[yr],U[yr])
        bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
        bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
        x = io.make_x(Z,bigY)
        L = io.make_L(Z,x)
        bigstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
        bigstressor[0:np.size(Y[yr],0),:] = stressor[yr]
        e = np.sum(bigstressor,1)/x 
        eL = np.dot(np.diag(e),L)
        reg = np.zeros((meta['reg']['len'],39))
        ind = np.zeros((meta['sup_dd']['len_idx'],39))
        
        for ysec in range (0,39):
            reg[:,ysec] = np.dot(np.dot(eL,bigY[:,ysec]),regagg)
            ind[:,ysec] = np.dot(np.dot(eL,bigY[:,ysec]),secagg)
            
        if indicator == 'ghg':
            drct[:,i] = direct.loc[:,yr]
        if indicator == 'co2':
            drct[:,i] = direct.loc[:,yr]
        if indicator == 'nrg':
            drct[:,i] = direct.loc[:,yr]
        if indicator == 'wat':
            drct[:,i] = direct.loc[:,yr]
        
        defra_foot[str(yr)+'_reg'] = df(reg, index = meta['reg']['idx'], columns = Y[yr].columns[0:np.size(Y[yr],1)-1])
        defra_foot[str(yr)+'_ind'] = df(ind, index = meta['sectors']['ind'], columns = Y[yr].columns[0:np.size(Y[yr],1)-1])
    
    defra_foot['direct'] = df(drct, index = direct.index, columns = allyears)
        
    for i,yr in enumerate(years): 
        print(yr)
        Z = io.make_Z_from_S_U(S[yr],U[yr])
        bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
        bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
        newbigY = np.zeros(shape = [np.size(newY[yr],0)*2,np.size(newY[yr],1)])
        newbigY[np.size(Y[yr],0):np.size(newY[yr],0)*2,0:np.size(newY[yr],1)] = newY[yr]
        x = io.make_x(Z,bigY)
        L = io.make_L(Z,x)
        bigstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
        bigstressor[0:np.size(Y[yr],0),:] = stressor[yr]
        e = np.sum(bigstressor,1)/x 
        eL_vector = np.dot(e,L)
      
        for ysec2 in range (0,np.size(newY[yr],1)-1):
            coicop[ysec2,i] = np.dot(eL_vector,newbigY[:,ysec2])
        if indicator == 'ghg':
            coicop[62,i] = coicop[62,i] + direct.loc['Consumer expenditure - not travel',yr]
            coicop[80,i] = coicop[80,i] + direct.loc['Consumer expenditure - travel',yr]
        if indicator == 'co2':
            coicop[62,i] = coicop[62,i] + direct.loc['Consumer expenditure - not travel',yr]
            coicop[80,i] = coicop[80,i] + direct.loc['Consumer expenditure - travel',yr]
        if indicator == 'nrg':
            coicop[61,i] = coicop[61,i] + direct.loc['domestic',yr]*0.223
            coicop[62,i] = coicop[62,i] + direct.loc['domestic',yr]*0.643
            coicop[63,i] = coicop[63,i] + direct.loc['domestic',yr]*0.134
            
            coicop[80,i] = coicop[80,i] + direct.loc['private transport',yr]
        if indicator == 'wat':
            coicop[60,i] = coicop[60,i] + direct[yr]
                  
        coicop2[0:134,i] = (coicop[0:134,i]*1000000)/(coicop_exp_tot[yr].values*52)
        coicop2[134,i] =  (coicop[134,i]*1000000)/(np.sum(Y[yr].loc[:,'Non-profit\r\ninstns serving\r\nhouseholds'],0)*1000000)
        coicop2[135,i] =  (coicop[135,i]*1000000)/(np.sum(Y[yr].loc[:,'Central\r\ngovernment'],0)*1000000)
        coicop2[136,i] =  (coicop[136,i]*1000000)/(np.sum(Y[yr].loc[:,'Local\r\nAuthorities'],0)*1000000)
        coicop2[137,i] =  (coicop[137,i]*1000000)/(np.sum(Y[yr].loc[:,'Gross fixed\r\ncapital\r\nformation'],0)*1000000)
        coicop2[138,i] =  (coicop[138,i]*1000000)/(np.sum(Y[yr].loc[:,'Valuables'],0)*1000000)
        coicop2[139,i] =  (coicop[139,i]*1000000)/(np.sum(Y[yr].loc[:,'Changes in inventories'],0)*1000000)
    
    defra_foot['coicop'] = df(coicop, index = newY[yr].columns[0:np.size(newY[yr],1)-1], columns = years)
    defra_foot['coicop_conv'] = df(coicop2, index = newY[yr].columns[0:np.size(newY[yr],1)-1], columns = years)
        
    return defra_foot
         
def makeregionresultscc(S,U,Y,newY,Yregion,pop_region,oa_pop,name,meta,stressor,direct,indicator,years,concs_dict3,conv):
    defra_foot = {}
    tempregagg = np.zeros((meta['fd']['len_idx'],meta['reg']['len'])) 
    tempsecagg = np.zeros((meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    for r in range(0,meta['reg']['len']):
        tempregagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],r] = np.transpose(np.ones((meta['sup_dd']['len_idx'],1)))
        tempsecagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],:] = np.identity(meta['sup_dd']['len_idx'])
    regagg = np.zeros((2*meta['fd']['len_idx'],meta['reg']['len'])) 
    secagg = np.zeros((2*meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
    regagg[0:meta['fd']['len_idx'],:] = tempregagg
    secagg[0:meta['fd']['len_idx'],:] = tempsecagg
    coicop = np.zeros((140,len(years)))
    drct = np.zeros((np.size(direct,0),len(years)))
    pop = np.zeros((2,len(years)))
  
    for i, yr in enumerate(years):
        Z = io.make_Z_from_S_U(S[yr],U[yr])
        bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
        bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
        bigYregion = np.zeros(shape = [np.size(Yregion[yr],0)*2,np.size(Yregion[yr],1)])
        bigYregion[np.size(Yregion[yr],0):np.size(Yregion[yr],0)*2,0:np.size(Yregion[yr],1)] = Yregion[yr]       
        x = io.make_x(Z,bigY)
        L = io.make_L(Z,x)
        bigstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
        bigstressor[0:np.size(Y[yr],0),:] = stressor[yr]
        e = np.sum(bigstressor,1)/x
        eL_vector = np.dot(e,L)
        eL = np.dot(np.diag(e),L)
        reg = np.zeros((meta['reg']['len'],140))
        ind = np.zeros((meta['sup_dd']['len_idx'],140))
      
        for ysec in range (0,np.size(Yregion[yr],1)):
            reg[:,ysec] = np.dot(np.dot(eL,bigYregion[:,ysec]),regagg)
            ind[:,ysec] = np.dot(np.dot(eL,bigYregion[:,ysec]),secagg)
            coicop[ysec,i] = np.dot(eL_vector,bigYregion[:,ysec])
        
        defra_foot[str(yr)+'_reg'] = df(np.dot(reg,concs_dict3['C140_to_C39']), index = meta['reg']['idx'], columns = Y[yr].columns[0:np.size(Y[yr],1)-1])
        defra_foot[str(yr)+'_ind'] = df(np.dot(ind,concs_dict3['C140_to_C39']), index = meta['sectors']['ind'], columns = Y[yr].columns[0:np.size(Y[yr],1)-1])
    
        if indicator == 'ghg':
            drct[0,i] = direct.loc['Consumer expenditure - not travel',yr]*np.sum(bigYregion[:,62],0)/np.sum(newY[yr].iloc[:,62],0)
            drct[1,i] = direct.loc['Consumer expenditure - travel',yr]*np.sum(bigYregion[:,80],0)/np.sum(newY[yr].iloc[:,80],0)
            coicop[62,i] = coicop[62,i] + drct[0,i]
            coicop[80,i] = coicop[80,i] + drct[1,i]
        if indicator == 'co2':
            drct[0,i] = direct.loc['Consumer expenditure - not travel',yr]*np.sum(bigYregion[:,62],0)/np.sum(newY[yr].iloc[:,62],0)
            drct[1,i] = direct.loc['Consumer expenditure - travel',yr]*np.sum(bigYregion[:,80],0)/np.sum(newY[yr].iloc[:,80],0)
            coicop[62,i] = coicop[62,i] + drct[0,i]
            coicop[80,i] = coicop[80,i] + drct[1,i]
        if indicator == 'nrg':
            drct[0,i] = direct.loc['domestic',yr]*np.sum(np.sum(bigYregion[:,61:64]))/np.sum(np.sum(newY[yr].iloc[:,61:64]))
            drct[1,i] = direct.loc['private transport',yr]*np.sum(bigYregion[:,80],0)/np.sum(newY[yr].iloc[:,80],0)
            coicop[61,i] = coicop[61,i] + drct[0,i]*0.223
            coicop[62,i] = coicop[62,i] + drct[0,i]*0.643
            coicop[63,i] = coicop[63,i] + drct[0,i]*0.134
            
            coicop[80,i] = coicop[80,i] + drct[1,i]
        if indicator == 'wat':
            drct[:,i] = direct.loc[:,yr]*np.sum(bigYregion[:,60],0)/np.sum(newY[yr].iloc[:,60],0)
            coicop[60,i] = coicop[60,i] + drct[:,i] 
        pop[0,i] = pop_region[yr].loc['Pop']
        pop[1,i] = np.sum(oa_pop[yr])
        
    defra_foot['direct'] = df(drct, index = direct.index, columns = years)
    defra_foot['coicop'] = df(coicop, index = Yregion[yr].columns[0:np.size(Yregion[yr],1)], columns = years)
    defra_foot['coicop_conv'] = conv
    defra_foot['pop'] = df(pop, index = [name,'UK'],columns = years)
    

    return defra_foot

def makelaresults(defra_region_result,reglaspropyr,regpophholdsyr,LAcodesnames,regoacsyr,oacyrmeta,oac01_names,oac11_names,Region_Name,years,concs_dict2):
    la_foot = {}
    la_names = reglaspropyr[2018][Region_Name].index
    for la in la_names:
        coicop = np.zeros((307,len(years)))
        ind = np.zeros((112,33))
        prd = np.zeros((112,33))
        reg = np.zeros((15,33))
        drct = np.zeros((2,len(years)))
        pop = np.zeros((3,len(years)))
        la_foot_temp = {}
        for i, yr in enumerate(years):
            coicop[:,i] =  defra_region_result['coicop'].iloc[0:307,i].values*reglaspropyr[yr][Region_Name].loc[la].values  
            coicop33prop = np.dot(coicop[:,i],concs_dict2['C307_to_C33'])/np.sum(defra_region_result[str(yr)+'_reg'].iloc[:,0:33])
            ind = defra_region_result[str(yr)+'_ind'].iloc[:,0:33].values*coicop33prop.values
            prd = defra_region_result[str(yr)+'_prd'].iloc[:,0:33].values*coicop33prop.values
            reg = defra_region_result[str(yr)+'_reg'].iloc[:,0:33].values*coicop33prop.values
            la_foot_temp[str(yr)+'_ind'] = df(ind, columns = coicop33prop.index, index = defra_region_result[str(yr)+'_ind'].index)
            la_foot_temp[str(yr)+'_prd'] = df(prd, columns = coicop33prop.index, index = defra_region_result[str(yr)+'_prd'].index)
            la_foot_temp[str(yr)+'_reg'] = df(reg, columns = coicop33prop.index, index = defra_region_result[str(yr)+'_reg'].index)
            drct[0,i] = defra_region_result['direct'].loc['Consumer expenditure - not travel',yr]*reglaspropyr[yr][Region_Name].loc[la,'4.4.2']
            drct[1,i] = defra_region_result['direct'].loc['Consumer expenditure - travel',yr]*reglaspropyr[yr][Region_Name].loc[la,'7.2.2.1']
            pop[1,i] = np.sum(regpophholdsyr[yr][Region_Name].loc[:,'pop'])
            pop[2,i] = regpophholdsyr[yr][Region_Name].loc[la,'pop']
            oac=df(regoacsyr[yr][Region_Name].loc[la,'pop'])
            if yr < 2014:
                if yr == 2008:
                     for i, c in enumerate(regoacsyr[yr][Region_Name].loc[la].index):
                         oac.loc[c,'OAC name'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                         if oacyrmeta[yr][Region_Name].loc[c,'oac3_count'] > 19:
                            oac.loc[c,'OAC code used'] = c
                            oac.loc[c,'OAC name used'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                         elif oacyrmeta[yr][Region_Name].loc[c,'oac2_count'] > 19:
                            oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].loc[c,'oac2'].lower()
                            oac.loc[c,'OAC name used'] = oac01_names['Group'].loc[oacyrmeta[yr][Region_Name].loc[c,'oac2'].lower(), 'Group name']
                         elif oacyrmeta[yr][Region_Name].loc[c,'oac1_count'] > 19:
                            oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].loc[c,'oac1']
                            oac.loc[c,'OAC name used'] = oac01_names['Supergroup'].loc[int(oacyrmeta[yr][Region_Name].loc[c,'oac1']), 'Supergroup name']
                         else:
                            oac.loc[c,'OAC code used'] = Region_Name
                            oac.loc[c,'OAC name used'] = [Region_Name +' average']
                else:
                    for i, c in enumerate(regoacsyr[yr][Region_Name].loc[la].index):
                        oac.loc[c,'OAC name'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                        if oacyrmeta[yr][Region_Name].loc[c.upper(),'oac3_count'] > 19:
                            oac.loc[c,'OAC code used'] = c
                            oac.loc[c,'OAC name used'] = oac01_names['Subgroup'].loc[c, 'Subgroup name']
                        elif oacyrmeta[yr][Region_Name].loc[c.upper(),'oac2_count'] > 19:
                            oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].loc[c.upper(),'oac2'].lower()
                            oac.loc[c,'OAC name used'] = oac01_names['Group'].loc[oacyrmeta[yr][Region_Name].loc[c.upper(),'oac2'].lower(), 'Group name']
                        elif oacyrmeta[yr][Region_Name].loc[c.upper(),'oac1_count'] > 19:
                            oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].loc[c.upper(),'oac1']
                            oac.loc[c,'OAC name used'] = oac01_names['Supergroup'].loc[int(oacyrmeta[yr][Region_Name].loc[c.upper(),'oac1']), 'Supergroup name']
                        else:
                            oac.loc[c,'OAC code used'] = Region_Name 
                            oac.loc[c,'OAC name used'] = [Region_Name +' average']                            
            else:
                 for i, c in enumerate(regoacsyr[yr][Region_Name].loc[la].index):
                        oac.loc[c,'OAC name'] = oac11_names['Subgroup'].loc[c, 'Subgroup name']
                        if oacyrmeta[yr][Region_Name].loc[c.upper(),'oac3_count'] > 19:
                            oac.loc[c,'OAC code used'] = c
                            oac.loc[c,'OAC name used'] = oac11_names['Subgroup'].loc[c, 'Subgroup name']
                        elif oacyrmeta[yr][Region_Name].loc[c.upper(),'oac2_count'] > 19:
                            oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].loc[c.upper(),'oac2'].lower()
                            oac.loc[c,'OAC name used'] = oac11_names['Group'].loc[oacyrmeta[yr][Region_Name].loc[c.upper(),'oac2'].lower(), 'Group name']
                        elif oacyrmeta[yr][Region_Name].loc[c.upper(),'oac1_count'] > 19:
                            oac.loc[c,'OAC code used'] = oacyrmeta[yr][Region_Name].loc[c.upper(),'oac1']
                            oac.loc[c,'OAC name used'] = oac11_names['Supergroup'].loc[int(oacyrmeta[yr][Region_Name].loc[c.upper(),'oac1']), 'Supergroup name']
                        else:
                            oac.loc[c,'OAC code used'] = Region_Name 
                            oac.loc[c,'OAC name used'] = [Region_Name +' average']   
        
            la_foot_temp[str(yr)+'_oac'] = oac
        pop[0] = [59113000, 	 59365700, 	 59636700, 	 59950400, 	 60413300, 	 60827100, 	 61319100, 	 61823800, 	 62260500, 	 62759500, 	 63285100, 	 63705000, 	 64105700, 	 64596800, 	 65110000, 	 65648100, 	 66040200, 	 66435600]
        
        la_foot_temp['coicop'] = df(coicop, index = defra_region_result['coicop'].iloc[0:307].index, columns = years)
        la_foot_temp['direct'] = df(drct, index = defra_region_result['direct'].index, columns = years)
        la_foot_temp['population'] = df(pop, index = ['United Kingdom',Region_Name,LAcodesnames.loc[la,'la_name']], columns = years)        
        
        la_foot[LAcodesnames.loc[la,'la_name']] = la_foot_temp
        
    return la_foot

def regioncheck(defra_ghg_north_east,defra_ghg_north_west,defra_ghg_yorkshire_and_the_humber,defra_ghg_east_midlands,defra_ghg_west_midlands,defra_ghg_east,defra_ghg_london,defra_ghg_south_east,defra_ghg_south_west,defra_ghg_scotland,defra_ghg_wales,defra_ghg_northern_ireland,years):
    check={}
    for yr in years:
        tempcheck =  np.zeros((12,39))
        tempcheck[0,:] = np.sum(defra_ghg_north_east[str(yr)+'_prd']).values
        tempcheck[1,:] = np.sum(defra_ghg_north_west[str(yr)+'_prd']).values
        tempcheck[2,:] = np.sum(defra_ghg_yorkshire_and_the_humber[str(yr)+'_prd']).values
        tempcheck[3,:] = np.sum(defra_ghg_east_midlands[str(yr)+'_prd']).values
        tempcheck[4,:] = np.sum(defra_ghg_west_midlands[str(yr)+'_prd']).values
        tempcheck[5,:] = np.sum(defra_ghg_east[str(yr)+'_prd']).values
        tempcheck[6,:] = np.sum(defra_ghg_london[str(yr)+'_prd']).values
        tempcheck[7,:] = np.sum(defra_ghg_south_east[str(yr)+'_prd']).values
        tempcheck[8,:] = np.sum(defra_ghg_south_west[str(yr)+'_prd']).values
        tempcheck[9,:] = np.sum(defra_ghg_scotland[str(yr)+'_prd']).values
        tempcheck[10,:] = np.sum(defra_ghg_wales[str(yr)+'_prd']).values
        tempcheck[11,:] = np.sum(defra_ghg_northern_ireland[str(yr)+'_prd']).values
        tempcheck = df(tempcheck, columns = defra_ghg_north_east[str(yr)+'_prd'].columns, index = ['North East','North West','Yorkshire and The Humber','East Midlands','West Midlands','East','London','South East','South West','Scotland','Wales','Northern Ireland'])
        
        tempcheck.loc['North East']['Electricity, gas and other fuels'] = tempcheck.loc['North East']['Electricity, gas and other fuels']+defra_ghg_north_east['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['North West']['Electricity, gas and other fuels'] = tempcheck.loc['North West']['Electricity, gas and other fuels']+defra_ghg_north_west['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['Yorkshire and The Humber']['Electricity, gas and other fuels'] = tempcheck.loc['Yorkshire and The Humber']['Electricity, gas and other fuels']+defra_ghg_yorkshire_and_the_humber['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['East Midlands']['Electricity, gas and other fuels'] = tempcheck.loc['East Midlands']['Electricity, gas and other fuels']+defra_ghg_east_midlands['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['West Midlands']['Electricity, gas and other fuels'] = tempcheck.loc['West Midlands']['Electricity, gas and other fuels']+defra_ghg_west_midlands['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['East']['Electricity, gas and other fuels'] = tempcheck.loc['East']['Electricity, gas and other fuels']+defra_ghg_east['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['London']['Electricity, gas and other fuels'] = tempcheck.loc['London']['Electricity, gas and other fuels']+defra_ghg_london['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['South East']['Electricity, gas and other fuels'] = tempcheck.loc['South East']['Electricity, gas and other fuels']+defra_ghg_south_east['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['South West']['Electricity, gas and other fuels'] = tempcheck.loc['South West']['Electricity, gas and other fuels']+defra_ghg_south_west['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['Scotland']['Electricity, gas and other fuels'] = tempcheck.loc['Scotland']['Electricity, gas and other fuels']+defra_ghg_scotland['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['Wales']['Electricity, gas and other fuels'] = tempcheck.loc['Wales']['Electricity, gas and other fuels']+defra_ghg_wales['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['Northern Ireland']['Electricity, gas and other fuels'] = tempcheck.loc['Northern Ireland']['Electricity, gas and other fuels']+defra_ghg_northern_ireland['direct'].loc['Consumer expenditure - not travel', yr]
        
        tempcheck.loc['North East']['Operation of personal transport equipment'] = tempcheck.loc['North East']['Operation of personal transport equipment']+defra_ghg_north_east['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['North West']['Operation of personal transport equipment'] = tempcheck.loc['North West']['Operation of personal transport equipment']+defra_ghg_north_west['direct'].loc['Consumer expenditure - not travel', yr]
        tempcheck.loc['Yorkshire and The Humber']['Operation of personal transport equipment'] = tempcheck.loc['Yorkshire and The Humber']['Operation of personal transport equipment']+defra_ghg_yorkshire_and_the_humber['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['East Midlands']['Operation of personal transport equipment'] = tempcheck.loc['East Midlands']['Operation of personal transport equipment']+defra_ghg_east_midlands['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['West Midlands']['Operation of personal transport equipment'] = tempcheck.loc['West Midlands']['Operation of personal transport equipment']+defra_ghg_west_midlands['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['East']['Operation of personal transport equipment'] = tempcheck.loc['East']['Operation of personal transport equipment']+defra_ghg_east['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['London']['Operation of personal transport equipment'] = tempcheck.loc['London']['Operation of personal transport equipment']+defra_ghg_london['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['South East']['Operation of personal transport equipment'] = tempcheck.loc['South East']['Operation of personal transport equipment']+defra_ghg_south_east['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['South West']['Operation of personal transport equipment'] = tempcheck.loc['South West']['Operation of personal transport equipment']+defra_ghg_south_west['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['Scotland']['Operation of personal transport equipment'] = tempcheck.loc['Scotland']['Operation of personal transport equipment']+defra_ghg_scotland['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['Wales']['Operation of personal transport equipment'] = tempcheck.loc['Wales']['Operation of personal transport equipment']+defra_ghg_wales['direct'].loc['Consumer expenditure - travel', yr]
        tempcheck.loc['Northern Ireland']['Operation of personal transport equipment'] = tempcheck.loc['Northern Ireland']['Operation of personal transport equipment']+defra_ghg_northern_ireland['direct'].loc['Consumer expenditure - travel', yr]
    
        check[yr] = tempcheck
    return check

def lacheck(defra_region,region_las,years):
    check={}
    temppop =  np.zeros((len(region_las),len(years)))
    for j, yr in enumerate(years):
        tempcheck =  np.zeros((len(region_las),33))
        label = []
        for i, la in enumerate(region_las):
            tempcheck[i,:] = np.sum(region_las[la][str(yr)+'_prd']).values
            tempcheck[i,10] = tempcheck[i,10]+region_las[la]['direct'].loc['Consumer expenditure - not travel', yr]        
            tempcheck[i,20] = tempcheck[i,20]+region_las[la]['direct'].loc['Consumer expenditure - travel', yr]      
            label.append(la)
            temppop[i,j] = region_las[la]['population'].loc[la,yr]
        tempcheck = df(tempcheck, columns = region_las[la][str(yr)+'_prd'].columns, index = label)        
        check[yr] = tempcheck       
    check['pop'] = df(temppop, columns = years, index = label)
    return check