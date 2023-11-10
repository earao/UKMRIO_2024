#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  21 14:01:34 2019

This code builds the UKMRIO database for the 2023 release (1990-2020) and the extension vectors. 

@authors: Anne Owen & Lena Kilian  
"""

import numpy as np
import pandas as pd
import ons as ons
import ukmrio_functions as uk
import nrg_stressor_functions as energy
import co2_stressor_functions as gas
import mat_stressor_functions as mat
import wat_stressor_functions as water
import metadata
import os
import pickle
from sys import platform
import copy as cp


# set working directory
# make different path depending on operating system
if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
mat_filepath = wd + 'data/raw/WU_materials/'
ons_filepath = wd + 'data/raw/ONS/'
iea_filepath = wd + 'data/raw/IEA/'
uk_energy_filepath = wd + 'data/processed/uk energy/'
inputs_filepath = wd + 'data/model_inputs/'
results_filepath = wd + 'outputs/results_2023/'
ons_name = 'SU_114_BB22_1997-2020.xlsx'
ons_year = '2023'
exiobase_filepath = wd + 'EXIOBASE/'
edgar_filepath = wd + 'data/raw/EDGAR/'

ayears = np.array([1990,1995,2005,2010,2013,2014,2015])
oldyrs = np.array(range(1990,1997))
exioyrs = np.array(range(1995,2021))
newyrs = np.array(range(1997,2021))
yrs = np.array(range(1990,2021))

(n_supply,n_use,n_final_demand,n_exports,dom_use,com_use,conc) = ons.load_io_data(wd, ons_filepath,newyrs,ons_year,ons_name)
(o_supply,o_use,o_final_demand,o_exports) = ons.load_old_io_data(ons_filepath)
(n_use,n_final_demand) = ons.remove_fd_negatives(n_use,n_final_demand,newyrs,112)

(dom_use,com_use) = ons.align_analytic_data(dom_use,com_use,conc)
 
for yr in range(1992,1997):
    o_final_demand[yr].loc['Total Intermediate consumption'] = np.sum(o_final_demand[yr],0)
    o_exports[yr].loc['Total Intermediate consumption'] = np.sum(o_exports[yr],0)
    
(o_supply,o_use,o_final_demand,o_exports) = ons.align_old_SUT_data(o_supply,o_use,o_final_demand,o_exports,n_supply[list(n_supply.keys())[0]],conc)

(o_use,o_final_demand) = ons.remove_fd_negatives(o_use,o_final_demand,oldyrs,112)

(supply,use,final_demand,exports) = ons.combine_data(o_supply,o_use,o_final_demand,o_exports,n_supply,n_use,n_final_demand,n_exports,oldyrs,newyrs)     

prod_sectors = supply[1997].columns
prod_use_sectors = use[1997]. columns
ind_sectors = supply[1997].index

file = os.path.join(inputs_filepath, 'EXIOBASE_data.xlsx')

c_conc = pd.read_excel(file, sheet_name = 'c_conc_15', index_col = 0)
i_conc = pd.read_excel(file, sheet_name = 'i_conc', header=[0,1], index_col=0)
p_conc = pd.read_excel(file, sheet_name = 'p_conc', header=[0,1], index_col=0)
    
regions = c_conc.columns

meta = metadata.make_meta(ind_sectors,prod_sectors,final_demand,regions)

(nS,nU,nY,nv) = uk.get_exiobase382(use,supply,exioyrs,meta,c_conc,i_conc,p_conc, exiobase_filepath)
(oS,oU,oY,ov) = uk.make_old_exio(nS,nU,nY,nv,inputs_filepath)
(S,U,Y,v) = uk.combine_exio(nS,nU,nY,nv,oS,oU,oY,ov,exioyrs)

for y in yrs:
     S[y].values[S[y].values<0] = 0
     U[y].values[U[y].values<0] = 0
     Y[y].values[Y[y].values<0] = 0 
     v[y].values[v[y].values<0] = 0          

Y = uk.split_y(Y,final_demand,yrs,meta)

file = os.path.join(inputs_filepath, 'currency.xlsx')
currency = pd.read_excel(file, sheet_name = 'currency', header = 4, index_col=0, usecols="A:B")

(S,U,Y,v) = uk.convert_to_gbp(S,U,Y,v,yrs,meta,currency)

domprop = uk.make_domprop(com_use,dom_use,use,final_demand,ayears,yrs,conc)

(imp_to_dom_prop, exp_fm_dom_prop, imp_to_dfd_prop) = uk.make_exio_props(U,Y,yrs,meta)

(domuse,rowuse,dom_dom_fd,ex_from_dom,imp_dom_fd) = uk.split_tables(yrs,domprop,use,supply,imp_to_dom_prop,meta,final_demand,exp_fm_dom_prop,exports,imp_to_dfd_prop)
  
(balancer,true_row_sum,true_col_sum) = uk.balancer_prep(use,supply,v,domuse,rowuse,dom_dom_fd,ex_from_dom,imp_dom_fd,Y,U,S,yrs,meta)

(U,v,Y,S) = uk.apply_ras50(balancer,true_row_sum,true_col_sum,yrs,meta,supply,S)

Y = uk.correctY(Y,final_demand,yrs)

#make emissions stressors
(exioCO2,exioGHG) = gas.make_exio382_stressor(use,exioyrs,exiobase_filepath,meta,c_conc,i_conc)
(exioCO2,exioGHG) = gas.make_old_exio_stressor_382(exioCO2,exioGHG,edgar_filepath,regions)
(uk_ghg_sectors,uk_co2_sectors,uk_ghg_direct,uk_co2_direct) = gas.make_UK_emissions(ons_filepath,yrs)
ghg = gas.make_ghg_382(exioGHG,uk_ghg_sectors,ons_filepath,yrs,meta)
co2 = gas.make_co2_382(exioCO2,uk_co2_sectors,ons_filepath,yrs,meta)

#make water stressors
(WATgrn_cons,WATblu_cons,WATblu_wdrl,uk_wat_blu_cons_direct,uk_wat_blu_wdrl_direct) = water.make_exio382_wat(use,exioyrs,meta,c_conc,i_conc,exiobase_filepath)
for yr in range(1990, 1996):
    WATgrn_cons[yr] = WATgrn_cons[1996]
    WATblu_cons[yr] = WATblu_cons[1996]
    WATblu_wdrl[yr] = WATblu_wdrl[1996]

    uk_wat_blu_cons_direct.loc[yr] = uk_wat_blu_cons_direct.loc[1996]
    uk_wat_blu_wdrl_direct.loc[yr] = uk_wat_blu_wdrl_direct.loc[1996]

#make energy stressor
(iea_fe_data,aviation,shipping) = energy.make_IEA_data(iea_filepath)
fullexioNRG = energy.iea_to_full_exio382(inputs_filepath,exiobase_filepath,iea_fe_data,aviation,shipping)
exioNRG = energy.uk_exio_nrg(fullexioNRG,use,yrs,meta,c_conc,i_conc) 
uk_nrg_direct = pd.read_excel(os.path.join(uk_energy_filepath,'UKenergy2023.xlsx'), sheet_name = 'direct')
nrg = energy.make_nrg_2023(uk_energy_filepath,exioNRG,S,yrs,meta)
       
#make material stressors
exioMAT = mat.make_exio_stressor_382(mat_filepath,yrs)
exioMAT = mat.make_UK_exioMAT(exioMAT,use,yrs,meta,c_conc,i_conc)
uk_MAT_sectors = mat.make_uk_stressor(ons_filepath,yrs)
(mat,bio,ore,nmm,ffl) = mat.make_mat(uk_MAT_sectors,exioMAT,S,yrs,meta)

# calculate footprints
uk_foot = {}
for item in ['ghg', 'co2', 'nrg', 'mat', 'bio', 'ore', 'nmm', 'ffl', 'WATblu_cons', 'WATblu_wdrl', 'WATgrn_cons']:
    uk_foot[item] = uk.footprint(eval(item), U, S, Y, yrs, meta)

#stressor data done

n_final_demand_hh = ons.load_hh_data(inputs_filepath,ons_filepath,final_demand,newyrs)
o_final_demand_hh = ons.make_old_fd_coicop(n_final_demand_hh,Y,oldyrs,meta)
final_demand_hh = ons.combine_fd_hh(o_final_demand_hh,n_final_demand_hh,oldyrs,newyrs)

final_demand_2 = ons.combine_fd(final_demand,final_demand_hh,yrs)

meta = metadata.make_meta(ind_sectors,prod_sectors,final_demand_2,regions)

hh_prop = ons.make_hh_prop(final_demand_hh,yrs)

Y2 = ons.make_wide_Y(Y,hh_prop,meta,yrs)
Y = Y2

#economic data done

uk_foot['ghg_Y2'] = uk.footprint(ghg,U,S,Y2,yrs,meta)

# save all files

# save as excel file
# save direct emissions
for data in ['uk_ghg_direct', 'uk_wat_blu_wdrl_direct', 'uk_wat_blu_cons_direct']:
    temp = cp.copy(eval(data))
    temp.to_excel(inputs_filepath + data + '.xlsx')
    
# save others
item_list = ['S', 'U', 'Y', 'co2', 'uk_co2_direct', 'ghg', 'mat', 'bio', 'ore', 'nmm', 'ffl', 
             'WATgrn_cons', 'WATblu_cons', 'WATblu_wdrl', 'nrg', 'uk_nrg_direct']

for data in item_list:
    stressor_writer = pd.ExcelWriter(inputs_filepath + data + '.xlsx')
    temp = cp.copy(eval(data))
    for yr in yrs:
        temp[yr].to_excel(stressor_writer, str(yr))
    stressor_writer.save() 

# save as pickle file
for data in item_list[:3]:
    pickle.dump(eval(data), open(inputs_filepath + data + ".p", "wb" ) )
    
for data in item_list[3:] + ['uk_ghg_direct', 'uk_wat_blu_wdrl_direct', 'uk_wat_blu_cons_direct']:
    pickle.dump(eval(data), open(inputs_filepath + data.lower() + ".p", "wb" ) )


# save meta
meta['fd']['rng_col']=slice(0, 43, None)  
meta['fd']['len_col']=43  
meta['fd_dd']['rng_col']=slice(0, 42, None)  
meta['fd_dd']['len_col']=42   

pickle.dump(meta, open(inputs_filepath + "meta.p", "wb" ) )