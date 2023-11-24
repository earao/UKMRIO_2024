#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:09:25 2018

This function is used with household_demand_main.py

@author: earao
"""

import pandas as pd
import numpy as np
df = pd.DataFrame

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


##########################################
# used in defra_uk_devolved_regions_2023 #
##########################################

# Also used in defra_main_2023 (see defra_main_2023):
# make_Yhh_109_34

#################################
# used in generations_2023_main #
#################################

# Also used in defra_main_2023 (see defra_main_2023):
# make_Yhh_109_34

    