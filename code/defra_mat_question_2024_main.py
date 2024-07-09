
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
allyears =  np.array([int(x) for x in range(1990,2022)])
oac_2001_years = np.array([int(x) for x in range(2007,2014)])
oac_2011_years = np.array([int(x) for x in range(2014,2022)])

S = pickle.load( open(results_filepath + "S.p", "rb" ) )
U = pickle.load( open(results_filepath + "U.p", "rb" ) )
Y = pickle.load( open(results_filepath + "Y.p", "rb" ) )
meta = pickle.load( open(results_filepath + "meta.p", "rb" ) )
ghg = pickle.load( open(results_filepath + "ghg.p", "rb" ) )
mat = pickle.load( open(results_filepath + "mat.p", "rb" ) )
uk_ghg_direct = pickle.load( open(results_filepath + "uk_ghg_direct.p", "rb" ) )

prod_conc = np.zeros((1680,2))
ind_conc = np.zeros((12,1680))

t_prod_conc = np.zeros((112,2))
t_ind_conc = np.zeros((6,112))

t_prod_conc[0:112,1] = np.ones((112))
t_prod_conc[58:61,0] = np.ones((3))
t_prod_conc[58:61,1] = np.zeros((3))

t_ind_conc[5,0:112] = np.ones((112))
t_ind_conc[0,3:6] = np.ones((3))
t_ind_conc[5,3:6] = np.zeros((3))
t_ind_conc[3,22] = np.ones((1))
t_ind_conc[5,22] = np.zeros((1))
t_ind_conc[4,25] = np.ones((1))
t_ind_conc[5,25] = np.zeros((1))
t_ind_conc[3,33] = np.ones((1))
t_ind_conc[5,33] = np.zeros((1))
t_ind_conc[2,34] = np.ones((1))
t_ind_conc[5,34] = np.zeros((1))
t_ind_conc[3,35:38] = np.ones((3))
t_ind_conc[5,35:38] = np.zeros((3))
t_ind_conc[4,52:54] = np.ones((2))
t_ind_conc[5,52:54] = np.zeros((2))
t_ind_conc[1,58:61] = np.ones((3))
t_ind_conc[5,58:61] = np.zeros((3))


for r in range(0,15):
    prod_conc[r*112:(r+1)*112,0:2] = t_prod_conc
    
ind_conc[0:6,0:112] = t_ind_conc

for r in range(1,15):
    ind_conc[6:12,r*112:(r+1)*112] = t_ind_conc

temp_ghg_cons_prod_by_src = np.zeros((1680,96))
temp_mat_cons_prod_by_src = np.zeros((1680,96))
temp_ghg_cons_ind_by_dest = np.zeros((96,1681))

int_ghg = np.zeros((2,32))
int_mat = np.zeros((2,32))

for y,yr in enumerate(range(1990,2022)):
    Z = io.make_Z_from_S_U(S[yr],U[yr])
    bigY = np.zeros(shape = [np.size(Y[yr],0)*2,np.size(Y[yr],1)])
    bigY[np.size(Y[yr],0):np.size(Y[yr],0)*2,0:np.size(Y[yr],1)] = Y[yr]
    x = io.make_x(Z,bigY)
    L = io.make_L(Z,x)
    A = io.make_A(Z,x)
    A_uk_const_zeroed = A
    A_uk_const_zeroed[58:61,:] = np.zeros((3,3360))
    A_uk_const_zeroed[58+1680:61+1680,:] = np.zeros((3,3360))
    A_uk_const_zeroed[:,58:61] = np.zeros((3360,3))
    A_uk_const_zeroed[:,58+1680:61+1680] = np.zeros((3360,3))
    L_zeros = np.linalg.inv(np.identity(len(L))-A_uk_const_zeroed)
    bigghgstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
    bigmatstressor = np.zeros(shape = [np.size(Y[yr],0)*2,1])
    bigghgstressor[0:np.size(Y[yr],0),:] = ghg[yr]
    bigmatstressor[0:np.size(Y[yr],0),:] = mat[yr]
    eghg = np.sum(bigghgstressor,1)/x 
    eghgL = np.dot(np.diag(eghg),L)
    eghgL_zeros= np.dot(np.diag(eghg),L_zeros)
    emat = np.sum(bigmatstressor,1)/x 
    ematL = np.dot(np.diag(emat),L)

    ghgfoot = df(np.dot(eghgL,np.diag(np.sum(bigY[:,0:42],1)))[0:1680,1680:3360],index=S[2020].index,columns=S[2020].columns)
    ghgfoot_zeros = df(np.dot(eghgL_zeros,np.diag(np.sum(bigY[:,0:42],1)))[0:1680,1680:3360],index=S[2020].index,columns=S[2020].columns)
    ghgfoot_no_cons = ghgfoot-ghgfoot_zeros
    ghgfoot_no_cons.iloc[58:61,:] = np.zeros((3,1680))
    ghgfoot_no_cons.iloc[:,58:61] = np.zeros((1680,3))
    
    int_ghg[0,y] = np.sum(np.sum(ghgfoot_no_cons))
    
    matfoot = df(np.dot(ematL,np.diag(np.sum(bigY[:,0:42],1)))[0:1680,1680:3360],index=S[2020].index,columns=S[2020].columns)
    
    ghgexport = df(np.dot(eghgL,np.diag(bigY[:,42]))[0:1680,1680:3360],index=S[2020].index,columns=S[2020].columns)
    ghgexport_zeros = df(np.dot(eghgL_zeros,np.diag(bigY[:,42]))[0:1680,1680:3360],index=S[2020].index,columns=S[2020].columns)
    ghgexport_no_cons = ghgexport-ghgexport_zeros
    ghgexport_no_cons.iloc[58:61,:] = np.zeros((3,1680))
    ghgexport_no_cons.iloc[:,58:61] = np.zeros((1680,3))
    
    int_ghg[1,y] = np.sum(np.sum(ghgexport_no_cons))
    
    matexport = df(np.dot(ematL,np.diag(bigY[:,42]))[0:1680,1680:3360],index=S[2020].index,columns=S[2020].columns)
    
    temp_ghg_cons_prod_by_src[:,y*3:(y+1)*3] = ghgfoot.iloc[:,58:61]
    temp_mat_cons_prod_by_src[:,y*3:(y+1)*3] = matfoot.iloc[:,58:61]
    temp_ghg_cons_ind_by_dest[y*3:(y+1)*3,0:1680] = ghgfoot.iloc[58:61,:]
    temp_ghg_cons_ind_by_dest[y*3:(y+1)*3,1680] = np.sum(ghgexport.iloc[58:61,:],1)

    temp_ghg_result = np.dot(np.dot(ind_conc,ghgfoot),prod_conc)
    temp_mat_result = np.dot(np.dot(ind_conc,matfoot),prod_conc)
    
    ghg_result[yr] = df(temp_ghg_result)
    mat_result[yr] = df(temp_mat_result)

UKscUKc = np.sum(ghgfoot['UK Construction Of Buildings'])

UKscUKc = np.sum(ghgfoot.iloc[58:61,58:61])
UK 
 for ysec in range (0,40):
     reg[:,ysec] = np.dot(np.dot(eL,bigY[:,ysec]),regagg)
     prd[:,ysec] = np.dot(np.dot(np.sum(eL,0),np.diag(bigY[:,ysec])),prdagg)
 ccc = np.dot(np.transpose(secagg),np.dot(np.dot(eL,np.diag(np.sum(bigY[:,0:-1],1))),prdagg))
 sic[:,i] = np.sum(prd,1)/ygg
