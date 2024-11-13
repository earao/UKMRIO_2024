
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 14:33:17 2022
 
@author: earao
"""
import pandas as pd
import numpy as np
df = pd.DataFrame
import pickle
import io_functions as io
from sys import platform
import os
import LCF_functions as lcf


# set working directory
# make different path depending on operating system
if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
inputs_filepath = wd + 'data/model_inputs/'
# results_filepath = wd + 'outputs/results_2023/'
results_filepath = wd + 'outputs/results_2024/'
lcf_filepath =  wd +'data/raw/LCFS/'
deflator_filepath =  wd +'data/raw/ONS/ONS deflators/'
census_filepath =  wd +'data/raw/census data/'
energy_filepath = wd + 'data/processed/uk energy/'

years = np.array([int(x) for x in range(1990,2022)])

S = pickle.load( open(results_filepath + "S.p", "rb" ) )
U = pickle.load( open(results_filepath + "U.p", "rb" ) )
Y = pickle.load( open(results_filepath + "Y.p", "rb" ) )
meta = pickle.load( open(results_filepath + "meta.p", "rb" ) )
ghg = pickle.load( open(results_filepath + "ghg.p", "rb" ) )
uk_ghg_direct = pickle.load( open(results_filepath + "uk_ghg_direct.p", "rb" ) )

        
concs_dict = pd.read_excel(os.path.join(inputs_filepath, 'ONS_to_COICOP_LCF_concs_2024.xlsx'), sheet_name=None, header = 0, index_col=0)

Y2 = lcf.convert43to41(Y,concs_dict,years)

Y=Y2

hhold_foot_cc = {}
hhold_foot_sic = {}
npish_foot = {}
cengv_foot = {}
locgv_foot = {}
gfixc_foot = {}
valua_foot = {}
chinv_foot = {}
total_foot = {}

for yr in years:
    print(yr)
    Z = io.make_Z_from_S_U(S[yr],U[yr])
    bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
    bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
    x = io.make_x(Z,bigY)
    L = io.make_L(Z,x)
    bigghgstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
    bigghgstressor[0:np.size(Y[yr],0),:] = ghg[yr]
    eghg = np.sum(bigghgstressor,1)/x 
    eghgL = np.dot(np.diag(eghg),L)
    temp_hhold_cc = np.zeros((3360,34))
    for i in range(0,34):
        temp_hhold_cc[:,i] = np.sum(np.dot(eghgL,np.diag(bigY[:,i])),1)
    hhold_foot_cc[yr] = df(temp_hhold_cc[0:1680,:],index = S[yr].index,columns = Y[yr].columns[0:34])
    
    temp_hhold_sic = np.dot(eghgL,np.diag(np.sum(bigY[:,0:34],1)))
    hhold_foot_sic[yr] = df(temp_hhold_sic[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    temp_npish = np.dot(eghgL,np.diag(bigY[:,34]))
    npish_foot[yr] = df(temp_npish[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    temp_cengv = np.dot(eghgL,np.diag(bigY[:,35]))
    cengv_foot[yr] = df(temp_cengv[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    temp_locgv = np.dot(eghgL,np.diag(bigY[:,36]))
    locgv_foot[yr] = df(temp_locgv[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    temp_gfixc = np.dot(eghgL,np.diag(bigY[:,37]))
    gfixc_foot[yr] = df(temp_gfixc[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    temp_valua = np.dot(eghgL,np.diag(bigY[:,38]))
    valua_foot[yr] = df(temp_valua[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    temp_chinv = np.dot(eghgL,np.diag(bigY[:,39]))
    chinv_foot[yr] = df(temp_chinv[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    temp_total = np.dot(eghgL,np.diag(np.sum(bigY[:,0:40],1)))
    total_foot[yr] = df(temp_total[0:1680,1680:],index = S[yr].index,columns = S[yr].columns)
    

for yr in years:
    hhold_foot_cc[yr].to_csv(results_filepath + str(yr) + '_hhold_cc.csv')
    hhold_foot_sic[yr].to_csv(results_filepath + str(yr) + '_hhold_sic.csv')
    npish_foot[yr].to_csv(results_filepath + str(yr) + '_npish.csv')
    cengv_foot[yr].to_csv(results_filepath + str(yr) + '_cengv.csv')
    locgv_foot[yr].to_csv(results_filepath + str(yr) + '_locgv.csv')
    gfixc_foot[yr].to_csv(results_filepath + str(yr) + '_gfixc.csv')
    valua_foot[yr].to_csv(results_filepath + str(yr) + '_valua.csv')
    chinv_foot[yr].to_csv(results_filepath + str(yr) + '_chinv.csv')
    total_foot[yr].to_csv(results_filepath + str(yr) + '_total.csv')
    