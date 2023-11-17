# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 17:01:48 2019

@author: earao
"""

import numpy as np
import pandas as pd
import os
df = pd.DataFrame

############################
# used in ukmrio_main_2024 #
############################

def make_exio382_wat(use,exioyrs,meta,c_conc,i_conc,exiobase_filepath): # used in ukmrio_main_2024

    exioWATgrn_cons = {}
    exioWATblu_cons = {}
    exioWATblu_wdrl = {}
    uk_wat_blu_cons_direct = []
    uk_wat_blu_wdrl_direct = []
       
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
        
        filepath = exiobase_filepath + "/3.8.2/MRSUT_{}/".format(str(exioyrs[a]))
        stressor = pd.read_csv(os.path.join(filepath, 'F.txt'), sep='\t', header = [0,1], index_col = 0)
        
        tempwat_g_c = stressor.loc['Water Consumption Green - Agriculture']
        tempwat_b_c = stressor.loc['Water Consumption Blue - Total']
        tempwat_b_w = stressor.loc['Water Withdrawal Blue - Total']
        
        weightedEXIOUKconci = np.zeros(shape=(7987,meta['use']['len_col']));
         
        uk_output_i = df.sum(use[exioyrs[a]].iloc[meta['v_d']['rng'],:], 1)
                             
        for m in range (0,49):
            for n in range (0,meta['reg']['len']):       
                if c_conc.iloc[m,n] == 1:
                    num = np.transpose(np.dot(np.diag(uk_output_i),np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])))
                    den = np.transpose(np.tile(np.dot(uk_output_i,np.transpose(exioUKconci.iloc[eind_slc[m],uind_slc[n]])),(meta['use_dd']['len_col'],1)))                                      
                    weightedEXIOUKconci[eind_slc[m],uind_slc[n]] = np.divide(num,den, out=np.zeros_like(num), where=den!=0)
        
        exioWATgrn_cons[exioyrs[a]] = df(np.dot(tempwat_g_c,weightedEXIOUKconci), index = meta['v']['col'])
        exioWATblu_cons[exioyrs[a]] = df(np.dot(tempwat_b_c,weightedEXIOUKconci), index = meta['v']['col'])
        exioWATblu_wdrl[exioyrs[a]] = df(np.dot(tempwat_b_w,weightedEXIOUKconci), index = meta['v']['col'])
        
        hh_stressor = pd.read_csv(os.path.join(filepath, 'F_Y.txt'), sep='\t', header = [0,1], index_col = 0)
        tempwat_b_c_hh = hh_stressor.loc['Water Consumption Blue - Domestic', ('GB', 'Final consumption expenditure by households')]
        tempwat_b_w_hh = hh_stressor.loc['Water Withdrawal Blue - Domestic', ('GB', 'Final consumption expenditure by households')]
        
        uk_wat_blu_cons_direct.append(tempwat_b_c_hh)
        uk_wat_blu_wdrl_direct.append(tempwat_b_w_hh)
    
    uk_wat_blu_cons_direct = df(uk_wat_blu_cons_direct,index = exioyrs)
    uk_wat_blu_wdrl_direct = df(uk_wat_blu_wdrl_direct,index = exioyrs)
        
    return(exioWATgrn_cons,exioWATblu_cons,exioWATblu_wdrl,uk_wat_blu_cons_direct,uk_wat_blu_wdrl_direct)