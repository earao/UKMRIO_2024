# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:04:43 2024

@author: earao
"""

import pandas as pd 
import numpy as np
df = pd.DataFrame
import pickle
import io_functions as io
from sys import platform

# set working directory
# make different path depending on operating system
if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
inputs_filepath = wd + 'data/model_inputs/'
results_filepath = wd + 'outputs/results_2024/'
lcf_filepath =  wd +'data/raw/LCFS/'
deflator_filepath =  wd +'data/raw/ONS/ONS deflators/'
census_filepath =  wd +'data/raw/census data/'
energy_filepath = wd + 'data/processed/uk energy/'

years = np.array([int(x) for x in range(2001,2022)])

# load the data
S = pickle.load( open(results_filepath + "S.p", "rb" ) )
U = pickle.load( open(results_filepath + "U.p", "rb" ) )
Y = pickle.load( open(results_filepath + "Y.p", "rb" ) )
meta = pickle.load( open(results_filepath + "meta.p", "rb" ) )
ghg = pickle.load( open(results_filepath + "ghg.p", "rb" ) )

# make an aggregator that compresses the data from all regions into a single set of 112 sectors
tempsecagg = np.zeros((meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
for r in range(0,meta['reg']['len']):
    tempsecagg[r*meta['sup_dd']['len_idx']:(r+1)*meta['sup_dd']['len_idx'],:] = np.identity(meta['sup_dd']['len_idx'])
prdagg = np.zeros((2*meta['fd']['len_idx'],meta['sup_dd']['len_idx']))
prdagg[meta['fd']['len_idx']:,:] = tempsecagg

# for loop set to run for 2021 only
for y,yr in enumerate(range(2021,2022)):
    Z = io.make_Z_from_S_U(S[yr],U[yr])
    bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
    bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
    x = io.make_x(Z,bigY)
    L = io.make_L(Z,x)
    bigghgstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
    bigghgstressor[0:np.size(Y[yr],0),:] = ghg[yr]
    eghg = np.sum(bigghgstressor,1)/x 
    
    # this would be the multipliers for all sectors in all regions - we want to combine all regions together
    eghgL = np.dot(np.diag(eghg),L)
    
    # make uk footprint
    ghgfoot = np.dot(eghgL,np.diag(np.sum(bigY[:,0:42],1)))
    
    # aggregate to 112 sectors
    ghgfoot_112sic = np.sum(np.dot(ghgfoot,prdagg),0)
    
    # aggregate Y to  112 sectors
    y_112sic = np.sum(np.dot(np.transpose(prdagg),bigY)[:,0:42],1)
    
    # divide emissions by final damand to make average multiplier by 112 sectors
    ghg_sic_mult = ghgfoot_112sic/y_112sic
    
    
    

   