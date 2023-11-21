#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 15:57:19 2018

These functions are used with the UKMRIO_main.py code to help build the UKMRIO database

@author: earao
"""
import numpy as np
import pandas as pd
from collections import OrderedDict
import os
df = pd.DataFrame

############################
# used in ukmrio_main_2024 #
############################

def get_exiobase382(use,supply,exioyrs,meta,c_conc,i_conc,p_conc, exiobase_filepath): # used in ukmrio_main_2024

    S = {}
    U = {}
    Y = {}
    v = {}
     
    exio_i_index = []
    for c in range (0,len(c_conc.index)):
        for i in range (0,len(i_conc.index)):            
            exio_i_index.append(" ".join([c_conc.index[c], i_conc.index[i]]))
                
    exio_p_index = []
    for c in range (0,len(c_conc.index)):
        for p in range (0,len(p_conc.index)):            
            exio_p_index.append(" ".join([c_conc.index[c], p_conc.index[p]]))
    
    eprd_slc = {}
    eind_slc = {}
    efd_slc = {}
    
    for c in range(0,49):
        eprd_slc[c] = slice(c*200,(c+1)*200)
        eind_slc[c] = slice(c*163,(c+1)*163)
        efd_slc[c] = slice(c*7,(c+1)*7)
    
    uprd_slc = {}
    uind_slc = {}
    
    for c in range(0,meta['reg']['len']):
        uprd_slc[c] = slice(c*meta['sup_dd']['len_col'],(c+1)*meta['sup_dd']['len_col'])
        uind_slc[c] = slice(c*meta['use_dd']['len_col'],(c+1)*meta['use_dd']['len_col'])
  
        
    idata = np.zeros(shape = (len(exio_i_index),meta['use']['len_col']))
    pdata = np.zeros(shape = (len(exio_p_index),meta['sup']['len_col']))
   
    for i in range (0,len(c_conc.index)):
        for c in range (0,len(c_conc.columns)):
            if c_conc.loc[c_conc.index[i]][c_conc.columns[c]] == 1:
                idata[i*len(i_conc):(i+1)*len(i_conc),uind_slc[c]] = i_conc.values
                pdata[i*len(p_conc):(i+1)*len(p_conc),uprd_slc[c]] = p_conc.values
     
    exioUKconci = df(idata, index = exio_i_index, columns = meta['use']['col'])
    exioUKconcp = df(pdata, index = exio_p_index, columns = meta['sup']['col'])
    
                      
    for a in range(0,len(exioyrs)):
        
        print(a)
                      
        filepath = exiobase_filepath + "3.8.2/MRSUT_{}/".format(str(exioyrs[a]))
                
        exio_s = pd.read_csv(os.path.join(filepath, 'supply.csv'), sep='\t', header = [0,1], index_col = [0,1])
        exio_u = pd.read_csv(os.path.join(filepath, 'use.csv'), sep='\t', header = [0,1], index_col = [0,1])
        exio_y = pd.read_csv(os.path.join(filepath, 'final_demand.csv'), sep='\t', header = [0,1], index_col = [0,1])
        exio_v = pd.read_csv(os.path.join(filepath, 'value_added.csv'), sep='\t', header = [0,1], index_col = 0)
        exio_v = df.sum(exio_v.iloc[0:12,:], 0)
        
#        if exioyrs[a]> 2013:
#            exio_u.loc['US','Other Bituminous Coal','EUR']['GB','Production of electricity by coal','EUR']=0
#           
        weightedEXIOUKconci = np.zeros(shape = (len(exio_i_index),meta['use']['len_col']))
        weightedEXIOUKconcp = np.zeros(shape = (len(exio_p_index),meta['sup']['len_col']))
        
        uk_output_i = df.sum(use[(exioyrs[a])].iloc[meta['v_d']['rng'],:], 1)
        uk_output_p = df.sum(supply[(exioyrs[a])], 1)
    
        for m in range (0,len(c_conc.index)):
            for n in range (0,meta['reg']['len']):       
                if c_conc.iloc[m,n] == 1:
                    num = np.transpose(np.dot(np.diag(uk_output_i),np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])))
                    den = np.transpose(np.tile(np.dot(uk_output_i,np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])),(meta['use_dd']['len_col'],1)))                                      
                    weightedEXIOUKconci[eind_slc[m],uind_slc[n]] = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
                    
                    num = np.transpose(np.dot(np.diag(uk_output_p),np.transpose(exioUKconcp.iloc[eprd_slc[m],uprd_slc[n]])))
                    den = np.transpose(np.tile(np.dot(uk_output_p,np.transpose(exioUKconcp.iloc[eprd_slc[m],uprd_slc[n]])),(meta['sup_dd']['len_col'],1)))                                      
                    weightedEXIOUKconcp[eprd_slc[m],uprd_slc[n]] = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
                    
        #exioUKregconc = np.zeros(shape=(len(np.transpose(exio_y)),2*meta['fd_dd']['len_col']))
        exioUKregconc = np.zeros(shape=(len(np.transpose(exio_y)),2*7))
        
        for m in range (0,len(c_conc.index)):
            for n in range (0,meta['reg']['len']):       
                if c_conc.iloc[m,0] == 1:
                    exioUKregconc[efd_slc[m],0:7] = np.identity(7)
                    #exioUKregconc[efd_slc[m],0:meta['fd_dd']['len_col']] = np.identity(meta['fd_dd']['len_col'])
                else:
                    exioUKregconc[efd_slc[m],7:7*2] = np.identity(7)
                    #exioUKregconc[efd_slc[m],meta['fd_dd']['len_col']:meta['fd_dd']['len_col']*2] = np.identity(meta['fd_dd']['len_col'])
     
        EXIOUKsupply = np.dot(np.transpose(weightedEXIOUKconci),np.dot(np.transpose(exio_s),weightedEXIOUKconcp))
        EXIOUKuse = np.dot(np.transpose(weightedEXIOUKconcp),np.dot(exio_u,weightedEXIOUKconci))
        EXIOUKV = np.dot(np.transpose(exio_v),weightedEXIOUKconci)
        EXIOUKY = np.dot(np.transpose(weightedEXIOUKconcp),np.dot(exio_y,exioUKregconc))  
            
        S[exioyrs[a]] = df(EXIOUKsupply/1000000, index = meta['sup']['idx'], columns = meta['sup']['col'])
        U[exioyrs[a]] = df(EXIOUKuse/1000000, index = meta['use']['idx'], columns = meta['use']['col'])       
        Y[exioyrs[a]] = df(EXIOUKY/1000000, index = meta['sup']['idx'])
        v[exioyrs[a]] = df(EXIOUKV/1000000, index = meta['v']['col'])     
                    
    return (S,U,Y,v)

def make_old_exio(nS,nU,nY,nv,filepath): # used in ukmrio_main_2023
    
    oS = {}
    oU = {}
    oY = {}
    ov = {}
    
    file = os.path.join(filepath, 'currency.xlsx')
    gdp = pd.read_excel(file, sheet_name = 'GDP', header = 4, index_col=0, usecols="A:E")

    oS[1994] = nS[1995]*gdp.loc[1994,'multiplier']
    oS[1993] = oS[1994]*gdp.loc[1993,'multiplier']
    oS[1992] = oS[1993]*gdp.loc[1992,'multiplier']
    oS[1991] = oS[1992]*gdp.loc[1991,'multiplier']
    oS[1990] = oS[1991]*gdp.loc[1990,'multiplier']
    
    oU[1994] = nU[1995]*gdp.loc[1994,'multiplier']
    oU[1993] = oU[1994]*gdp.loc[1993,'multiplier']
    oU[1992] = oU[1993]*gdp.loc[1992,'multiplier']
    oU[1991] = oU[1992]*gdp.loc[1991,'multiplier']
    oU[1990] = oU[1991]*gdp.loc[1990,'multiplier']
    
    oY[1994] = nY[1995]*gdp.loc[1994,'multiplier']
    oY[1993] = oY[1994]*gdp.loc[1993,'multiplier']
    oY[1992] = oY[1993]*gdp.loc[1992,'multiplier']
    oY[1991] = oY[1992]*gdp.loc[1991,'multiplier']
    oY[1990] = oY[1991]*gdp.loc[1990,'multiplier']
    
    ov[1994] = nv[1995]*gdp.loc[1994,'multiplier']
    ov[1993] = ov[1994]*gdp.loc[1993,'multiplier']
    ov[1992] = ov[1993]*gdp.loc[1992,'multiplier']
    ov[1991] = ov[1992]*gdp.loc[1991,'multiplier']
    ov[1990] = ov[1991]*gdp.loc[1990,'multiplier']
    
    return (oS,oU,oY,ov)

def combine_exio(nS,nU,nY,nv,oS,oU,oY,ov,exioyrs): # used in ukmrio_main_2023
    
    S = {}
    U = {}
    Y = {}
    v = {}

    for yr in range(1990,1995):
        S[yr] = oS[yr]
        U[yr] = oU[yr]
        Y[yr] = oY[yr]
        v[yr] = ov[yr]
        
    for yr in exioyrs:
        S[yr] = nS[yr]
        U[yr] = nU[yr]
        Y[yr] = nY[yr]
        v[yr] = nv[yr]
    
    return (S,U,Y,v)

def split_y(Y,final_demand,yrs,meta): # used in ukmrio_main_2023
    
    for a in range(0,len(yrs)):
        temp = np.zeros(shape=(meta['fd']['len_idx'],meta['fd']['len_col']));
    
        temp[:,0:2] = Y[yrs[a]].iloc[:,0:2];                        
        temp[:,2] = np.multiply(np.divide(df.sum(final_demand[yrs[a]].iloc[:,2]),(df.sum(df.sum(final_demand[yrs[a]].iloc[:,2:4])))),Y[yrs[a]].iloc[:,2]);
        temp[:,3] = np.multiply(np.divide(df.sum(final_demand[yrs[a]].iloc[:,3]),(df.sum(df.sum(final_demand[yrs[a]].iloc[:,2:4])))),Y[yrs[a]].iloc[:,2]);
        temp[:,4:7] = Y[yrs[a]].iloc[:,3:6];
        temp[:,7] = df.sum(Y[yrs[a]].iloc[:,7:14],1)
                
        Y[yrs[a]] = df(temp, index = meta['fd']['idx'],columns = meta['fd']['col']) 
        
    return(Y)

def convert_to_gbp(S,U,Y,v,yrs,meta,currency): # used in ukmrio_main_2023
    
    vtemp = np.zeros(shape=(meta['v']['len']))
    Ytemp = np.zeros(shape=(meta['fd']['len_idx'],meta['fd']['len_col']))
    Stemp = np.zeros(shape=(meta['sup']['len_idx'],meta['sup']['len_col']))
    Utemp = np.zeros(shape=(meta['use']['len_idx'],meta['use']['len_col']))
    
    for a in range(0,len(yrs)):
        vtemp = v[yrs[a]]*np.tile(currency.loc[yrs[a]],(meta['v']['len'],1))
        Ytemp = Y[yrs[a]]*np.tile(currency.loc[yrs[a]],(meta['fd']['len_idx'],meta['fd']['len_col']))
        Stemp = S[yrs[a]]*np.tile(currency.loc[yrs[a]],(meta['sup']['len_idx'],meta['sup']['len_col']))
        Utemp = U[yrs[a]]*np.tile(currency.loc[yrs[a]],(meta['use']['len_idx'],meta['use']['len_col']))
             
        Y[yrs[a]] = df(Ytemp, index = meta['fd']['idx'], columns = meta['fd']['col']) 
        S[yrs[a]] = df(Stemp, index = meta['sup']['idx'], columns = meta['sup']['col'])
        U[yrs[a]] = df(Utemp, index = meta['use']['idx'], columns = meta['use']['col'])
        v[yrs[a]] = df(np.sum(vtemp,1), index = meta['v']['col'])
              
    return (S,U,Y,v)
        
def make_exio_props(U,Y,yrs,meta): # used in ukmrio_main_2023
    
    imp_to_dom_prop = np.zeros(shape=(meta['use_id']['len_idx'],meta['use_id']['len_col'],len(yrs)))
    exp_fm_dom_prop = np.zeros(shape=(meta['use_id']['len_col'],meta['use_id']['len_idx']+1,len(yrs)))
    imp_to_dfd_prop = np.zeros(shape=(meta['use_id']['len_idx'],7,len(yrs)))
       
    for a in range(0,len(yrs)):  
        # the col sum of all imports to UK domestic from EXIO
        imp_to_dom_sum = df.sum(U[yrs[a]].iloc[meta['use_id']['rng']], 0)
        np.nan_to_num(imp_to_dom_sum, copy = False)
        
        # the row sum of all exports from UK domestic from EXIO
        exp_fm_dom_sum = np.sum(U[yrs[a]].iloc[meta['use_di']['rng']],1)+np.sum(Y[yrs[a]].iloc[meta['fd_di']['rng']],1)
        
        # the sum of FD not from UK from EXIO
        imp_to_dfd_sum = df.sum(df.sum(Y[yrs[a]].iloc[112:,0:7]))
        np.nan_to_num(imp_to_dfd_sum, copy = False)
               
        num = U[yrs[a]].iloc[meta['use_id']['rng']]
        den = np.tile(imp_to_dom_sum,(meta['use_id']['len_idx'],1))
        imp_to_dom_prop[:,:,a] = np.divide(num,den,out = np.zeros_like(num), where = den!=0)
        
        num = np.zeros((112,1569))
        num[:,0:1568] = U[yrs[a]].iloc[meta['use_di']['rng']]
        num[:,1568:1569] = Y[yrs[a]].iloc[meta['fd_di']['rng']]
        den = np.transpose(np.tile(exp_fm_dom_sum,(meta['use_id']['len_idx']+1,1)))
        exp_fm_dom_prop[:,:,a] = np.divide(num,den,out = np.zeros_like(num), where = den!=0)
        
        num = Y[yrs[a]].iloc[112:,0:7]
        den = np.tile(imp_to_dfd_sum,(1568,7))
        imp_to_dfd_prop[:,:,a] = np.divide(num,den,out = np.zeros_like(num), where = den!=0)
                               
    return(imp_to_dom_prop, exp_fm_dom_prop, imp_to_dfd_prop)

def split_tables(yrs,domprop,use,supply,imp_to_dom_prop,meta,final_demand,exp_fm_dom_prop,exports,imp_to_dfd_prop): # used in ukmrio_main_2023
 
   domuse = {}
   rowuse = {}
   dom_dom_fd = {}
   ex_from_dom = {}
   imp_dom_fd = {}

   for a in range(0,len(yrs)):
       temp_dd = np.multiply(domprop[str(yrs[a])].iloc[0:112,0:112],use[yrs[a]].iloc[0:112,0:112])
       imports = np.sum(supply[yrs[a]].values,1)-np.sum(temp_dd.values,0)-use[yrs[a]].loc['GVA (production measure)'].values
       imports[imports < 0] = 0.
       temp_di = np.multiply(imp_to_dom_prop[:,:,a],np.tile(imports,(meta['use_id']['len_idx'],1)))
       temp_ddfd = np.multiply(domprop[str(yrs[a])].iloc[0:112,112:119],final_demand[yrs[a]].iloc[0:112,0:7])
       temp_d_exp = np.multiply(exp_fm_dom_prop[:,:,a],np.transpose(np.tile(domprop[str(yrs[a])].iloc[0:112,119].values*exports[yrs[a]].iloc[0:112].values,(1569,1))))
       temp_idfd = np.multiply(imp_to_dfd_prop[:,:,a],np.tile(np.multiply(domprop[str(yrs[a])].iloc[112,112:119],final_demand[yrs[a]].iloc[112,0:7].values),(1568,1)))
       
       temp_dd[temp_dd < 0] = 0.000000001
       temp_di[temp_di < 0] = 0.000000001
       temp_ddfd[temp_ddfd < 0] = 0.000000001
       temp_d_exp[temp_d_exp < 0] = 0.000000001
       temp_idfd[temp_idfd < 0] = 0.000000001
       
       domuse[yrs[a]] = df(temp_dd)
       rowuse[yrs[a]] = df(temp_di)
       dom_dom_fd[yrs[a]] = df(temp_ddfd)
       ex_from_dom[yrs[a]] = df(temp_d_exp)
       imp_dom_fd[yrs[a]] = df(temp_idfd)
           
   return(domuse,rowuse,dom_dom_fd,ex_from_dom,imp_dom_fd)
    
def balancer_prep(use,supply,v,domuse,rowuse,dom_dom_fd,ex_from_dom,imp_dom_fd,Y,U,S,yrs,meta): # used in ukmrio_main_2023
   
    all_va = np.zeros(shape=(meta['v']['len']))
    all_va1 = np.zeros(shape=(meta['v']['len']))
    all_fd = np.zeros(shape=(meta['fd']['len_idx'],meta['fd']['len_col']))
    all_fd1 = np.zeros(shape=(meta['fd']['len_idx'],meta['fd']['len_col']))
    balancer = {}
    true_row_sum = {}  
    true_col_sum = {}

    for yr in yrs:
        
        tempbalancer = np.zeros((1681,1688))
        temptrue_row_sum = np.zeros(1681)   
        temptrue_col_sum = np.zeros(1688)
       
        all_va[meta['v_d']['rng']] = use[yr].loc['GVA (production measure)']
        all_va[meta['v_i']['rng']] = np.sum(v[yr].iloc[meta['v_i']['rng']].values,1)
        
        all_fd[meta['fd_dd']['rng']] = dom_dom_fd[yr]
        all_fd[meta['fd_di']['rng']] = Y[yr].iloc[0:112,7:8]
        all_fd[meta['fd_id']['rng']] = imp_dom_fd[yr]
        all_fd[meta['fd_ii']['rng']] = Y[yr].iloc[112:,7:8]
        
        prop_fd = np.divide(np.sum(all_fd[:,:]),np.sum(all_va[:]))
        all_va1[:] = np.multiply(all_va[:],np.tile(prop_fd,(meta['v']['len'])))
        
        prop_va = np.divide(np.sum(all_va[:]),np.sum(all_fd[:,:]))
        all_fd1[:] = np.multiply(all_fd[:,:],np.tile(prop_va,(meta['fd']['len_idx'],meta['fd']['len_col'])))
                 
        tempbalancer[0:112,0:112] = domuse[yr]
        tempbalancer[112:1680,0:112] = rowuse[yr]
        tempbalancer[0:112,112:1680] = ex_from_dom[yr].iloc[:,0:1568]
        tempbalancer[112:1680,112:1680] = U[yr].iloc[meta['use_ii']['rng']]
         
        tempbalancer[1680,0:1680] = np.transpose(all_va1)
        tempbalancer[0:1680,1680:1688] = all_fd[:,:]
        
        temptrue_col_sum[0:112] = np.sum(supply[yr],1)
        temptrue_col_sum[112:1680] = np.sum(S[yr].iloc[112:1680,112:1680],1)
        
        temptrue_row_sum[0:112] = np.sum(supply[yr],0)
        temptrue_row_sum[112:1680] = np.sum(S[yr].iloc[112:1680,112:1680],0)
        
        temptrue_col_sum[1680:1688] = np.sum(all_fd1,0)
        temptrue_row_sum[1680] = np.sum(all_va)
        
        np.nan_to_num(true_col_sum, copy = False) 
        np.nan_to_num(true_row_sum, copy = False) 
        
        balancer[yr] = df(tempbalancer)
        true_col_sum[yr] = df(temptrue_col_sum)
        true_row_sum[yr] = df(np.transpose(temptrue_row_sum))
  
    return (balancer,true_row_sum,true_col_sum)
   
    all_va = np.zeros(shape=(meta['v']['len']))
    all_va1 = np.zeros(shape=(meta['v']['len']))
    all_fd = np.zeros(shape=(meta['fd']['len_idx'],meta['fd']['len_col']))
    balancer = np.zeros(shape=(meta['bal']['len_idx'],meta['bal']['len_col'],len(yrs)))
    true_row_sum = np.zeros(shape =(meta['bal']['len_idx'],np.size(yrs)))   
    true_col_sum = np.zeros(shape =(np.size(yrs),meta['bal']['len_col']))  
    
    for a in range(0,np.size(yrs)):   
         
        all_va[meta['v_d']['rng']] = np.add(use[yrs[a]].loc['GVA (production measure)'],prd_tax[yrs[a]])
        all_va[meta['v_i']['rng']] = np.sum(v[yrs[a]].iloc[meta['v_i']['rng']].values,1)
        all_fd[meta['fd_dd']['rng']] = dom_dom_fd[yrs[a]]
        all_fd[meta['fd_id']['rng']] = dom_row_fd[yrs[a]]
        all_fd[meta['fd_di']['rng']] = Y[yrs[a]].iloc[meta['fd_di']['rng']]
        all_fd[meta['fd_ii']['rng']] = Y[yrs[a]].iloc[meta['fd_ii']['rng']]
        
        prop_fd = np.divide(np.sum(all_fd[:,:]),np.sum(all_va[:]))
        all_va1[:] = np.multiply(all_va[:],np.tile(prop_fd,(meta['v']['len'])))
                
        balancer[meta['bal_use_dd']['rng_idx'],meta['bal_use_dd']['rng_col'],a] = domuse[yrs[a]];
        balancer[meta['bal_use_id']['rng_idx'],meta['bal_use_id']['rng_col'],a] = rowuse[yrs[a]];
        balancer[meta['bal_use_di']['rng_idx'],meta['bal_use_di']['rng_col'],a] = ex_from_dom[yrs[a]];
        balancer[meta['bal_use_ii']['rng_idx'],meta['bal_use_ii']['rng_col'],a] = U[yrs[a]].iloc[meta['use_ii']['rng']];
        balancer[meta['bal_sup_dd']['rng_idx'],meta['bal_sup_dd']['rng_col'],a] = supply[yrs[a]];
        balancer[meta['bal_sup_ii']['rng_idx'],meta['bal_sup_ii']['rng_col'],a] = S[yrs[a]].iloc[meta['sup_ii']['rng']];
        
        balancer[meta['bal']['len_idx']-1,meta['bal_use']['rng_col'],a] = np.transpose(all_va1);
        balancer[meta['bal_fd']['rng_idx'],meta['bal_fd']['rng_col'],a] = all_fd[:,:]
        
        balancer[balancer < 0] = 0.000000001
        
        true_row_sum[meta['bal_sup_dd']['rng_idx'],a] = np.sum(supply[yrs[a]],1); 
        true_col_sum[a,meta['bal_use_dd']['rng_col']] = np.transpose(true_row_sum[meta['bal_sup_dd']['rng_idx'],a])
        
        true_row_sum[meta['bal_sup_ii']['rng_idx'],a] = np.sum(S[yrs[a]].iloc[meta['sup_ii']['rng']],1); 
        true_col_sum[a,meta['bal_use_ii']['rng_col']] = np.transpose(true_row_sum[meta['bal_sup_ii']['rng_idx'],a]);
        
        true_col_sum[a,meta['bal_sup_dd']['rng_col']] = np.sum(supply[yrs[a]],0);
        true_row_sum[meta['bal_use_dd']['rng_idx'],a] = np.transpose(true_col_sum[a,meta['bal_sup_dd']['rng_col']]);
        
        true_col_sum[a,meta['bal_sup_ii']['rng_col']] = np.sum(S[yrs[a]].iloc[meta['sup_ii']['rng']],0);
        true_row_sum[meta['bal_use_id']['rng_idx'],a] = np.transpose(true_col_sum[a,meta['bal_sup_ii']['rng_col']]);
       
        true_col_sum[a,meta['bal_fd']['rng_col']] = np.transpose(np.sum(all_fd,0))
        true_row_sum[meta['bal']['len_idx']-1,a] = np.sum(all_va1); 
        
        np.nan_to_num(true_col_sum, copy = False) 
        np.nan_to_num(true_row_sum, copy = False) 
  
    return (balancer,true_row_sum,true_col_sum)

def apply_ras50(balancer,true_row_sum,true_col_sum,yrs,meta,supply,S): # used in ukmrio_main_2023
    
    tU50={}
    tY50={}
    tv50={}
    tS50={}
        
    for yr in yrs:
        print(yr)
        true_row_sum[yr][true_row_sum[yr]<=0] = 1
        true_col_sum[yr][true_col_sum[yr]<=0] = 1
        to_balance = balancer[yr]
        trs = true_row_sum[yr]
        tcs = np.transpose(true_col_sum[yr])
        for _ in range(0, 50):
            col_sum_now = np.zeros(shape = (1,len(tcs)))
            col_sum_now = np.sum(to_balance,0)
            col_sum_now[col_sum_now<=0] = 1
                
            for row in range(0, len(trs)):
                to_balance.iloc[row,:] = np.multiply(np.divide(to_balance.iloc[row,:], col_sum_now), np.sum(tcs,0).values);
                    
            row_sum_now = np.sum(to_balance,1)
            row_sum_now[row_sum_now<=0] = 1
                
            for col in range(0, len(tcs)):
                to_balance.iloc[:,col] = np.multiply(np.divide(to_balance.iloc[:,col], row_sum_now), np.sum(trs,1).values)
     
        balanced = to_balance
        temp_tS = np.zeros((1680,1680))
        temp_tS[0:112,0:112] = supply[yr]
        temp_tS[112:,112:] = S[yr].iloc[112:,112:]
        
        tU50[yr] = df(balanced.iloc[0:1680,0:1680].values,index = meta['use']['idx'],columns = meta['use']['col'])
        tv50[yr] = df(balanced.iloc[1680,0:1680].values,index = meta['use']['col'])
        tY50[yr] = df(balanced.iloc[0:1680,1680:1688].values,index = meta['fd']['idx'],columns = meta['fd']['col'])
        tS50[yr] = df(temp_tS,index = meta['use']['col'],columns = meta['use']['idx'])
                
    return(tU50,tv50,tY50,tS50) 

def footprint(stressor,U,S,Y,yrs,meta): # used in ukmrio_main_2023
    
    f = np.zeros(shape=(meta['reg']['len'],np.size(yrs)))
    x = {}
    uind_slc = {}    
    for c in range(0,meta['reg']['len']):
         uind_slc[c] = slice(c*meta['use_dd']['len_col'],(c+1)*meta['use_dd']['len_col'])
    
    for a, yr in enumerate(yrs):
        Z = np.zeros(shape = (meta['Z']['len_idx'],meta['Z']['len_col']))
        long_Y = np.zeros(shape = (meta['Z']['len_idx'],meta['fd']['len_col']))
        long_stressor = np.zeros(shape = (meta['Z']['len_idx']))
        Z[meta['bal_use']['rng_idx'],meta['bal_use']['rng_col']] = U[yr]
        Z[meta['bal_sup']['rng_idx'],meta['bal_sup']['rng_col']] = S[yr]
        long_Y[meta['bal_fd']['rng_idx'],meta['fd']['rng_col']] = Y[yr]
        long_stressor[0:np.size(U[yr],0)] = stressor[yr].iloc[:,0]
        x = np.sum(Z,1)+np.sum(long_Y,1)
        x[x==0] = 0.000000001
        recix = 1/x
        diag_recix = np.diag(recix) 
        L = np.linalg.inv(np.eye(Z.shape[0])-np.dot(Z,diag_recix))
        intensity = long_stressor*recix
        bigfoot = np.dot(np.dot(np.diag(intensity),L),np.diag(np.sum(long_Y[:,meta['fd_dd']['rng_col']],1)))
        for b in range(0,meta['reg']['len']):
            f[b,a] = np.sum(bigfoot[uind_slc[b],:])

    footprint = df(f, columns = yrs)
    
    return(footprint)

def correctY(Y,final_demand,yrs):
    
    nY = {}
    
    for yr in yrs:
        
        tempY = Y[yr].values
        
        nonhhshares = np.divide(final_demand[yr].iloc[0:112,1:6].values,np.transpose(np.tile(np.sum(final_demand[yr].iloc[0:112,1:6],1),(5,1))))
        nonhhshares[np.isnan(nonhhshares)] = 0
        nonhhshares = np.tile(nonhhshares,(15,1))
        fdzerolist = np.transpose(np.tile(final_demand[yr].iloc[0:112,0],(1,15)))
        
        for n in range(0,1680):
            if fdzerolist[n] == 0:
                tempY[n,1] = tempY[n,1]+tempY[n,0]*nonhhshares[n,0]
                tempY[n,2] = tempY[n,2]+tempY[n,0]*nonhhshares[n,1]
                tempY[n,3] = tempY[n,3]+tempY[n,0]*nonhhshares[n,2]
                tempY[n,4] = tempY[n,4]+tempY[n,0]*nonhhshares[n,3]
                tempY[n,5] = tempY[n,5]+tempY[n,0]*nonhhshares[n,4]
                tempY[n,0] = 0
                
        nY[yr] = df(tempY, index = Y[yr].index, columns = Y[yr].columns)   
    
    return(nY)

def make_domprop(com_use,dom_use,use,final_demand,ayears,yrs,conc):
   
    tempdp_ayears = {}
    tempdp = {}
    temp = np.zeros((113,120))
    domprop = OrderedDict()
    
    
    for a in range(0,3):
        
        temp = np.zeros((113,120))
        
        conc['annxb_v'].columns = com_use[str(ayears[a])].index
        conc['annxb_h'].columns = com_use[str(ayears[a])].columns
        tempdom = df.dot(conc['annxb_v'],df.dot(dom_use[str(ayears[a])],df.transpose(conc['annxb_h'])))
        tempcom = df.dot(conc['annxb_v'],df.dot(com_use[str(ayears[a])],df.transpose(conc['annxb_h'])))
        
        temp[0:112,0:112] = tempdom.iloc[0:112,0:112]/tempcom.iloc[0:112,0:112]
        temp[0:112,112:119] = tempdom.iloc[0:112,113:120].values/tempcom.iloc[0:112,113:120].values
        temp[112:113,0:112] = np.transpose(tempdom.iloc[113,0:112].values/tempcom.iloc[112,0:112].values)
        temp[112:113,112:119] = np.transpose(tempdom.iloc[113,113:120].values/tempcom.iloc[112,113:120].values)
        temp[0:113,119] = tempdom.iloc[0:113,120]/tempcom.iloc[0:113,120].values
        
        temp = df(temp)
        
        temp.fillna(1, inplace=True)
        temp[temp<0]=0
        temp[temp>1]=1
        
        tempdp_ayears[ayears[a]]=temp
        
    for a in range(3,7):
        
        temp = np.zeros((113,120))
    
        temp[0:112,0:112] = dom_use[str(ayears[a])].iloc[0:112,0:112]/com_use[str(ayears[a])].iloc[0:112,0:112]
        temp[0:112,112:119] = dom_use[str(ayears[a])].iloc[0:112,113:120].values/com_use[str(ayears[a])].iloc[0:112,113:120].values
        temp[112:113,0:112] = np.transpose(dom_use[str(ayears[a])].iloc[113,0:112].values/com_use[str(ayears[a])].iloc[112,0:112].values)
        temp[112:113,112:119] = np.transpose(dom_use[str(ayears[a])].iloc[113,113:120].values/com_use[str(ayears[a])].iloc[112,113:120].values)
        temp[0:113,119] = dom_use[str(ayears[a])].iloc[0:113,120]/com_use[str(ayears[a])].iloc[0:113,120].values
    
        temp = df(temp)
        
        temp.fillna(1, inplace=True)
        temp[temp<0]=0
        temp[temp>1]=1
        
        tempdp_ayears[ayears[a]]=temp
        
           
    for a in range(0,len(ayears)-1):
        diff_in_prop = np.divide((tempdp_ayears[ayears[a+1]].values-tempdp_ayears[ayears[a]].values),(ayears[a+1]-ayears[a]))
        
        for b in range((ayears[a]-ayears[0]),ayears[a+1]-ayears[0]+1):
            year_mult = b-(ayears[a]-ayears[0])
            temp =  np.add(np.multiply(diff_in_prop,year_mult),tempdp_ayears[ayears[a]])
            tempdp[str(1990+b)] = temp
            del(temp)
    
    for a in range(0,(ayears[len(ayears)-1]-yrs[0])):       
        domprop[str(yrs[a])] = tempdp[str(yrs[a])]
    for a in range(ayears[len(ayears)-1]-yrs[0],len(yrs)):         
        domprop[str(yrs[a])] = tempdp[str(ayears[len(ayears)-1])]  
    
    domprop['1990']=domprop['1995']
    domprop['1991']=domprop['1995']
    domprop['1992']=domprop['1995']
    domprop['1993']=domprop['1995']
    domprop['1994']=domprop['1995']
        
    return domprop

def make_domprop2(com_use,dom_use,use,final_demand,ayears,yrs,conc):
   
    tempdp_ayears = {}
    tempdp = {}
    temp = np.zeros((113,120))
    domprop = OrderedDict()
    
    
    for a in range(0,3):
        
        temp = np.zeros((113,120))
        
        conc['annxb_v'].columns = com_use[str(ayears[a])].index
        conc['annxb_h'].columns = com_use[str(ayears[a])].columns
        tempdom = df.dot(conc['annxb_v'],df.dot(dom_use[str(ayears[a])],df.transpose(conc['annxb_h'])))
        tempcom = df.dot(conc['annxb_v'],df.dot(com_use[str(ayears[a])],df.transpose(conc['annxb_h'])))
        
        temp[0:112,0:112] = tempdom.iloc[0:112,0:112]/tempcom.iloc[0:112,0:112]
        temp[0:112,112:119] = tempdom.iloc[0:112,113:120].values/tempcom.iloc[0:112,113:120].values
        temp[112:113,0:112] = np.transpose(tempdom.iloc[113,0:112].values/tempcom.iloc[112,0:112].values)
        temp[112:113,112:119] = np.transpose(tempdom.iloc[113,113:120].values/tempcom.iloc[112,113:120].values)
        temp[0:113,119] = tempdom.iloc[0:113,120]/tempcom.iloc[0:113,120].values
        
        temp = df(temp)
        
        temp.fillna(1, inplace=True)
        temp[temp<0]=0
        temp[temp>1]=1
        
        tempdp_ayears[ayears[a]]=temp
        
    for a in range(3,7):
        
        temp = np.zeros((113,120))
    
        temp[0:112,0:112] = dom_use[str(ayears[a])].iloc[0:112,0:112]/com_use[str(ayears[a])].iloc[0:112,0:112]
        temp[0:112,112:119] = dom_use[str(ayears[a])].iloc[0:112,113:120].values/com_use[str(ayears[a])].iloc[0:112,113:120].values
        temp[112:113,0:112] = np.transpose(dom_use[str(ayears[a])].iloc[113,0:112].values/com_use[str(ayears[a])].iloc[112,0:112].values)
        temp[112:113,112:119] = np.transpose(dom_use[str(ayears[a])].iloc[113,113:120].values/com_use[str(ayears[a])].iloc[112,113:120].values)
        temp[0:113,119] = dom_use[str(ayears[a])].iloc[0:113,120]/com_use[str(ayears[a])].iloc[0:113,120].values
    
        temp = df(temp)
        
        temp.fillna(1, inplace=True)
        temp[temp<0]=0
        temp[temp>1]=1
        
        tempdp_ayears[ayears[a]]=temp
        
           
    for a in range(0,len(ayears)-1):
        diff_in_prop = np.divide((tempdp_ayears[ayears[a+1]].values-tempdp_ayears[ayears[a]].values),(ayears[a+1]-ayears[a]))
        
        for b in range((ayears[a]-ayears[0]),ayears[a+1]-ayears[0]+1):
            year_mult = b-(ayears[a]-ayears[0])
            temp =  np.add(np.multiply(diff_in_prop,year_mult),tempdp_ayears[ayears[a]])
            tempdp[str(1990+b)] = temp
            del(temp)
    
    for a in range(0,(ayears[len(ayears)-1]-yrs[0])):       
        domprop[str(yrs[a])] = tempdp[str(yrs[a])]
    for a in range(ayears[len(ayears)-1]-yrs[0],len(yrs)):         
        domprop[str(yrs[a])] = tempdp[str(ayears[len(ayears)-1])]  
    
    domprop['1990']=domprop['1995']
    domprop['1991']=domprop['1995']
    domprop['1992']=domprop['1995']
    domprop['1993']=domprop['1995']
    domprop['1994']=domprop['1995']
        
    return domprop

###########################
# used in defra_main_2023 #
###########################

def deflate_io_regions(S,U,Y,newY,y_regions,allyears,years,io_deflators,meta): # used in defra_main_2023
    
    S_d = {}
    U_d = {}
    Y_d = {}
    newY_d = {}
    y_regions_d = {}
          
    for yr in allyears:
        
        temp =np.diag(pd.Series.repeat(io_deflators.loc[int(yr)],meta['reg']['len']))
        io_deflators_s = df(temp, index = S[yr].columns, columns = S[yr].index)
        io_deflators_u = df(temp, index = U[yr].columns, columns = U[yr].index)
            
        S_d[yr] = df.dot(S[yr],io_deflators_s)
        U_d[yr] = df.dot(io_deflators_u,U[yr])
        Y_d[yr] = df.dot(io_deflators_u,Y[yr])
    
    for yr in years:
        
        temp =np.diag(pd.Series.repeat(io_deflators.loc[int(yr)],meta['reg']['len']))
        io_deflators_s = df(temp, index = S[yr].columns, columns = S[yr].index)
        io_deflators_u = df(temp, index = U[yr].columns, columns = U[yr].index)
        
        newY_d[yr] = df.dot(io_deflators_u,newY[yr])
        for reg in range(1,13):
            y_regions_d[str(yr)+'_'+str(reg)] = df.dot(io_deflators_u,y_regions[str(yr)+'_'+str(reg)])

    return (S_d,U_d,Y_d,newY_d,y_regions_d)
    
def make_v(U,Y,yrs,meta):  # used in defra_main_2023
    v = {}
    for a in range(0,np.size(yrs)):
        tempV = np.sum(U[yrs[a]],1)+np.sum(Y[yrs[a]],1)-np.sum(U[yrs[a]],0)
        v[yrs[a]] = df(tempV, index = meta['v']['col'])
        
    return v
        
def get_deflator_data_2023(deflator_filepath): # used in defra_main_2023
    
    file = deflator_filepath + 'deflators_AO2023.xlsx'
    io = pd.read_excel(file, sheet_name='UK_domestic', index_col=0)
    cc = pd.read_excel(file, sheet_name='COICOP')
    
    io_deflators = 1/io*100
    cc_deflators = 1/cc*100
       
    return (io_deflators,cc_deflators) 

##########################################
# used in defra_uk_devolved_regions_2023 #
##########################################

# Also used in defra_main_2023 (see defra_main_2023):
# get_deflator_data_2023
# deflate_io_regions
# make_v(U,Y,yrs,meta)
