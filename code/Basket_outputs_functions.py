#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 1 2023
"""
import numpy as np
from itertools import permutations
import math
import pandas as pd
df = pd.DataFrame

def make_total(data, years, idx_dict):
    total_data = pd.DataFrame()
    for year in years:
        temp = pd.DataFrame(data[year].sum(axis=0)).T
        temp['year'] = year
        total_data = total_data.append(temp)
    total_data = total_data.set_index('year').T.fillna(0)
    total_data.index = [x.split(' ')[0] for x in total_data.index]
    total_data = total_data.rename(index=idx_dict)
    
    return total_data

def import_cpi(data_filepath, idx_dict):  
    cpi = {}
    cpi['2008-2022'] = pd.read_excel(data_filepath + 'raw/CPI/202311_consumerpriceinflationdetailedreferencetables.xlsx', sheet_name='Table 9', index_col=2, header=5)\
        .drop(['Unnamed: 0', 'Unnamed: 1'], axis=1).dropna(how='all')
        
    cpi['2001-2007'] = pd.read_excel(data_filepath + 'raw/CPI/202311_consumerpriceinflationdetailedreferencetables_tcm77-4192423.xls', sheet_name='Table 9', index_col=2, header=5)\
        .drop(['Unnamed: 0', 'Unnamed: 1'], axis=1).dropna(how='all')
    
    lookup = pd.read_csv(data_filepath + 'lookups/cpi_to_lcfs.csv').set_index('CPI')
    
    for item in list(cpi.keys()):
        # macth lcfs columns
        cpi[item] = cpi[item].join(lookup, how='right').set_index('LCFS').mean(axis=0, level=0).fillna(0)
        # change to 2010 as base year
        cpi_comp = cpi[item][2010]
        cpi[item] = cpi[item].apply(lambda x: x / cpi_comp * 100)
    cpi = cpi['2001-2007'].loc[:, :2007].join(cpi['2008-2022']).fillna(0)
    cpi.index = [x.split(' ')[0] for x in cpi.index]
    cpi = cpi.rename(index=idx_dict)
    
    return cpi


def sda(sda_0,sda_1):
    
    terms = len(sda_0)
    template = np.zeros(shape = [math.factorial(terms),terms,terms])
    solutionmatrix = np.zeros(shape = [terms,terms])
    result = np.zeros(shape=[len(template)+5,terms+1])
    temp = {}
    combs_items = np.zeros(shape=terms)
    
    for a in range(0,terms):
        solutionmatrix[a,a] = 2
    
    for a in range(0,terms):
        for b in range(0,terms):
            if solutionmatrix[a,b]==2:
                solutionmatrix[a,b+1:terms] = 1

    for a in range(0,terms):
        combs_items[a] = a    
    
    combs = np.array(list(set(permutations(combs_items))))          
        
    for a in range(0,terms):        
        for b in range(0,math.factorial(terms)):
            for c in range(0,terms):
                
                if a==c:                    
                    template[b,c,a] = 2
                    
                else:                   
                    for d in range(0,terms):                        
                        for e in range(0,terms):                            
                            if combs[b,d] == a:                                
                                if combs[b,e] == c:                                    
                                    template[b,c,a]=solutionmatrix[d,e]
        
    for a in range(0,terms):     
        for b  in range(0,len(template)):            
            for c in range(0,terms):
                
                if template[b,c,a] == 2:                    
                    temp[c] = sda_0[c]-sda_1[c]
                    
                elif template[b,c,a] == 1:                    
                    temp[c] = sda_0[c]
                    
                elif template[b,c,a] == 0:                    
                    temp[c] = sda_1[c]            
            
            tempresults = 1
            
            for c in range(terms):                
               tempresults = np.dot(tempresults, temp[c])
               print(tempresults)
            
            result[b,a+1] = tempresults
                         
    for a in range(0,len(template)):        
        result[a,0] = np.sum(result[a,1:terms+1])
    
    for a in range(0,terms):
        
        result[len(template)+0,a+1]=np.mean(result[0:len(template),a+1])
        result[len(template)+1,a+1]=np.max(result[0:len(template),a+1])
        result[len(template)+2,a+1]=np.min(result[0:len(template),a+1])
        result[len(template)+3,a+1]=np.std(result[0:len(template),a+1])
        result[len(template)+4,a+1]=result[len(template),a+1]/result[0,0]*100
        
    result = pd.DataFrame(result)
    result.index = ['SDA_' + str(x) for x in range(len(template))] + ['mean', 'max', 'min', 'std', 'pct_mean']
    
    return result
