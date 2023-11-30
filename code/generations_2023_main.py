
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 14:33:17 2022
 
@author: earao
"""
import pandas as pd
import numpy as np
import os
df = pd.DataFrame
import LCF_functions as lcf
import pickle
import demand_functions as dm
import defra_functions as defra
import io_functions as io
import copy
from sys import platform

# set working directory
# make different path depending on operating system
if platform[:3] == 'win':
    wd = 'O://UKMRIO_Data/'
else:
    wd = r'/Volumes/a72/UKMRIO_Data/'

# define filepaths
inputs_filepath = wd + 'data/model_inputs/'
# results_filepath = wd + 'outputs/results_2023/'
results_filepath = wd + 'temp_anne_outputs/results_2023/'
lcf_filepath =  wd +'data/raw/LCFS/'
deflator_filepath =  wd +'data/raw/ONS/ONS deflators/'
census_filepath =  wd +'data/raw/census data/'
energy_filepath = wd + 'data/processed/uk energy/'

years = np.array([int(x) for x in range(2001,2021)])
allyears =  np.array([int(x) for x in range(1990,2021)])
oac_2001_years = np.array([int(x) for x in range(2007,2014)])
oac_2011_years = np.array([int(x) for x in range(2014,2021)])

S = pickle.load( open(inputs_filepath + "S.p", "rb" ) )
U = pickle.load( open(inputs_filepath + "U.p", "rb" ) )
Y = pickle.load( open(inputs_filepath + "Y.p", "rb" ) )
hhspenddata = pickle.load( open(inputs_filepath + "hhspenddata.p", "rb" ) )
meta = pickle.load( open(inputs_filepath + "meta.p", "rb" ) )
ghg = pickle.load( open(inputs_filepath + "ghg.p", "rb" ) )
uk_ghg_direct = pickle.load( open(inputs_filepath + "uk_ghg_direct.p", "rb" ) )

hhspenddata2 = copy.deepcopy(hhspenddata)
hhspenddata2 = lcf.removeoutliers(hhspenddata2,years)
        
concs_dict = pd.read_excel(os.path.join(inputs_filepath, 'ONS_to_COICOP_LCF_concs.xlsx'), sheet_name=None, header = 0, index_col=0)

Y2 = lcf.convert43to41(Y,concs_dict,allyears)

Y=Y2

total_Yhh_112 = dm.make_Yhh_112_34(Y,years,meta)

coicop_exp_tot = lcf.make_totals_2023(hhspenddata2,years)
coicop_exp_tot2 = lcf.convert_exp_tot_sizes(coicop_exp_tot, concs_dict, years, '456_to_105')
coicop_exp_tot3 = lcf.make_balanced_totals_2023(coicop_exp_tot2,total_Yhh_112,concs_dict,years)

yhh_wide = lcf.make_y_hh_105(Y,coicop_exp_tot3,years,concs_dict,meta)
    
newY = lcf.make_new_Y_105(Y,yhh_wide,years)

hhspenddata3 = lcf.convert_hhspend_sizes(hhspenddata2,concs_dict,years, '456_to_105')
for yr in list(hhspenddata3.keys()):
    hhspenddata3[yr]['OA class 1'] = hhspenddata3[yr]['OA class 1'].astype(str)

defra_ghg_uk = defra.makeukresults2023(wd,S,U,Y,newY,coicop_exp_tot2,meta,ghg,uk_ghg_direct,'ghg',allyears,years)
# ghg_mults = defra_ghg_uk['coicop_mult'].iloc[0:105,:]

coicop_deflators = pd.read_excel(deflator_filepath + 'coicop.xlsx', sheet_name = 'COICOP', index_col = 0)
coicop_deflators = coicop_deflators.loc[years,:]/100
 

hhspenddata3[2008]=hhspenddata3[2008].drop(index = 2788)
hhspenddata3[2008]=hhspenddata3[2008].drop(index = 89)
hhspenddata3[2020].loc[:,'4.2.1 Imputed rentals of owner occupiers']=hhspenddata3[2020].loc[:,'4.2.1 Imputed rentals of owner occupiers']/150

temp = np.zeros((105,20))
temp2 = np.zeros((105,20))
for i,yr in enumerate(years):
    temp[:,i] = (np.sum(hhspenddata3[yr].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.']))*52
    temp2[:,i] = defra_ghg_uk['coicop'].iloc[0:105,i]
ghg_mults = df(temp2/temp*1000,columns = years, index = hhspenddata3[2020].columns[20:125])

ghg_mults = ghg_mults.fillna(0)

ghg_mults_deflated = np.multiply(ghg_mults,np.transpose(coicop_deflators))

hhspenddata_deflated =  copy.deepcopy(hhspenddata3)
for year in years:
    temp = np.divide(hhspenddata3[year].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'],np.tile(coicop_deflators.loc[year,:],(np.size(hhspenddata3[year],0),1)))
    hhspenddata_deflated[year].loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'] = temp

footdata_def = lcf.processdataforfoots(ghg_mults_deflated,hhspenddata_deflated,years)
footdata = lcf.processdataforfoots(ghg_mults,hhspenddata3,years)

totalbygen = lcf.totalbygen(footdata,years)
totalbygen_d = lcf.totalbygen(footdata_def,years)  

pop_agecohort_foots = lcf.popagecohortfoots(footdata,years)
pop_agecohort_foots_d = lcf.popagecohortfoots(footdata_def,years)

# make SDA variables

temp = np.zeros((4,20))
for y, year in enumerate(years):
    temp2 = footdata[year].iloc[:,0:44].groupby(by='gen code').sum()     
    temp[:,y] = temp2.loc[:,'pop x weight']
    

pop = df(temp,columns = years, index = ['Silent and older','Babyboomer','Gen X','Millenial and younger'])

total_spend = np.zeros((4,20))
prop_spend_gen= {}
total_spend_all= {}
for y, year in enumerate(years):
    temp = hhspenddata_deflated[year]
    temp = temp.reset_index()
    for i in range(0,len(temp)):
        temp.loc[i,'birth year'] = year - temp.loc[i,'age hrp']
    for i in range(0,len(temp)):
        if(temp.loc[i,'birth year']>1981):
            temp.loc[i,'gen code'] = 3
        elif(temp.loc[i,'birth year']>1965):
            temp.loc[i,'gen code'] = 2
        elif(temp.loc[i,'birth year']>1946):
            temp.loc[i,'gen code'] = 1
        else:
            temp.loc[i,'gen code'] = 0
        
    spends = temp.iloc[:,21:128].groupby(by='gen code').sum()
    propspends = spends/np.transpose(np.tile(np.sum(spends.loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'],1),(106,1)))
    total_spend[:,y] = np.sum(spends.loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'],1)
    total_spend_all[year] = df(spends.loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'])*52   
    prop_spend_gen[year] = df(propspends.loc[:,'1.1.1 Bread and cereals':'12.7.1 Other services n.e.c.'])
total_spend = df(total_spend,columns = years, index = ['Silent and older','Babyboomer','Gen X','Millenial and younger'])
 
per_cap_spend = total_spend/pop*52

prop_spend = np.zeros((105,20)) 
new_prop_spend_gen= {} 

g=0
for y, year in enumerate(years):
    prop_spend[:,y] = np.transpose(prop_spend_gen[year].iloc[0,:])
new_prop_spend_gen['Silent and older'] = df(prop_spend, columns = years)

prop_spend = np.zeros((105,20))
g=1
for y, year in enumerate(years):
    prop_spend[:,y] = np.transpose(prop_spend_gen[year].iloc[1,:])
new_prop_spend_gen['Babyboomer'] = df(prop_spend, columns = years)

prop_spend = np.zeros((105,20))
g=2
for y, year in enumerate(years):
    prop_spend[:,y] = np.transpose(prop_spend_gen[year].iloc[2,:])
new_prop_spend_gen['Gen X'] = df(prop_spend, columns = years)

prop_spend = np.zeros((105,20))
g=3
for y, year in enumerate(years):     
    prop_spend[:,y] = np.transpose(prop_spend_gen[year].iloc[3,:])
new_prop_spend_gen['Millenial and younger'] = df(prop_spend, columns = years)


# SDA
# silent
sda_0 = {}
sda_1 = {}
silent = {}

for year in range(2002,2021):
    pop_0 = df(np.tile(pop.loc['Silent and older',2001],(1,105)))
    pop_1 = df(np.tile(pop.loc['Silent and older',year],(1,105)))

    exp_0 = df(np.diagflat(np.tile(per_cap_spend.loc['Silent and older',2001],(1,105))))
    exp_1 = df(np.diagflat(np.tile(per_cap_spend.loc['Silent and older',year],(1,105))))

    prp_0 = df(np.diag(new_prop_spend_gen['Silent and older'].loc[:,2001]))
    prp_1 = df(np.diag(new_prop_spend_gen['Silent and older'].loc[:,year]))

    ghg_0 = ghg_mults_deflated.loc[:,2001].values
    ghg_1 = ghg_mults_deflated.loc[:,year].values

    foot_0 = df.dot(pop_0,exp_0)
    foot_0 = df.dot(foot_0,prp_0)
    foot_0 = df.dot(foot_0,ghg_0)
    
    foot_1 = df.dot(pop_1,exp_1)
    foot_1 = df.dot(foot_1,prp_1)
    foot_1 = df.dot(foot_1,ghg_1)

    sda_0 = {}
    sda_0[0] = pop_0
    sda_0[1] = exp_0
    sda_0[2] = prp_0
    sda_0[3] = ghg_0
    
    sda_1 = {}
    sda_1[0] = pop_1
    sda_1[1] = exp_1
    sda_1[2] = prp_1
    sda_1[3] = ghg_1

    silent[year] = io.sda(sda_1,sda_0)
    
# babyboomer
sda_0 = {}
sda_1 = {}
babyboomer = {}

for year in range(2002,2021):
    pop_0 = df(np.tile(pop.loc['Babyboomer',2001],(1,105)))
    pop_1 = df(np.tile(pop.loc['Babyboomer',year],(1,105)))

    exp_0 = df(np.diagflat(np.tile(per_cap_spend.loc['Babyboomer',2001],(1,105))))
    exp_1 = df(np.diagflat(np.tile(per_cap_spend.loc['Babyboomer',year],(1,105))))

    prp_0 = df(np.diag(new_prop_spend_gen['Babyboomer'].loc[:,2001]))
    prp_1 = df(np.diag(new_prop_spend_gen['Babyboomer'].loc[:,year]))

    ghg_0 = ghg_mults_deflated.loc[:,2001].values
    ghg_1 = ghg_mults_deflated.loc[:,year].values

    foot_0 = df.dot(pop_0,exp_0)
    foot_0 = df.dot(foot_0,prp_0)
    foot_0 = df.dot(foot_0,ghg_0)
    
    foot_1 = df.dot(pop_1,exp_1)
    foot_1 = df.dot(foot_1,prp_1)
    foot_1 = df.dot(foot_1,ghg_1)

    sda_0 = {}
    sda_0[0] = pop_0
    sda_0[1] = exp_0
    sda_0[2] = prp_0
    sda_0[3] = ghg_0
    
    sda_1 = {}
    sda_1[0] = pop_1
    sda_1[1] = exp_1
    sda_1[2] = prp_1
    sda_1[3] = ghg_1

    babyboomer[year] = io.sda(sda_1,sda_0)

# Gen_x
sda_0 = {}
sda_1 = {}
Gen_x = {}

for year in range(2002,2021):
    pop_0 = df(np.tile(pop.loc['Gen X',2001],(1,105)))
    pop_1 = df(np.tile(pop.loc['Gen X',year],(1,105)))

    exp_0 = df(np.diagflat(np.tile(per_cap_spend.loc['Gen X',2001],(1,105))))
    exp_1 = df(np.diagflat(np.tile(per_cap_spend.loc['Gen X',year],(1,105))))

    prp_0 = df(np.diag(new_prop_spend_gen['Gen X'].loc[:,2001]))
    prp_1 = df(np.diag(new_prop_spend_gen['Gen X'].loc[:,year]))

    ghg_0 = ghg_mults_deflated.loc[:,2001].values
    ghg_1 = ghg_mults_deflated.loc[:,year].values

    foot_0 = df.dot(pop_0,exp_0)
    foot_0 = df.dot(foot_0,prp_0)
    foot_0 = df.dot(foot_0,ghg_0)
    
    foot_1 = df.dot(pop_1,exp_1)
    foot_1 = df.dot(foot_1,prp_1)
    foot_1 = df.dot(foot_1,ghg_1)

    sda_0 = {}
    sda_0[0] = pop_0
    sda_0[1] = exp_0
    sda_0[2] = prp_0
    sda_0[3] = ghg_0
    
    sda_1 = {}
    sda_1[0] = pop_1
    sda_1[1] = exp_1
    sda_1[2] = prp_1
    sda_1[3] = ghg_1

    Gen_x[year] = io.sda(sda_1,sda_0)
    
# Millenial
sda_0 = {}
sda_1 = {}
Millenial = {}

for year in range(2002,2021):
    pop_0 = df(np.tile(pop.loc['Millenial and younger',2001],(1,105)))
    pop_1 = df(np.tile(pop.loc['Millenial and younger',year],(1,105)))

    exp_0 = df(np.diagflat(np.tile(per_cap_spend.loc['Millenial and younger',2001],(1,105))))
    exp_1 = df(np.diagflat(np.tile(per_cap_spend.loc['Millenial and younger',year],(1,105))))

    prp_0 = df(np.diag(new_prop_spend_gen['Millenial and younger'].loc[:,2001]))
    prp_1 = df(np.diag(new_prop_spend_gen['Millenial and younger'].loc[:,year]))

    ghg_0 = ghg_mults_deflated.loc[:,2001].values
    ghg_1 = ghg_mults_deflated.loc[:,year].values

    foot_0 = df.dot(pop_0,exp_0)
    foot_0 = df.dot(foot_0,prp_0)
    foot_0 = df.dot(foot_0,ghg_0)
    
    foot_1 = df.dot(pop_1,exp_1)
    foot_1 = df.dot(foot_1,prp_1)
    foot_1 = df.dot(foot_1,ghg_1)

    sda_0 = {}
    sda_0[0] = pop_0
    sda_0[1] = exp_0
    sda_0[2] = prp_0
    sda_0[3] = ghg_0
    
    sda_1 = {}
    sda_1[0] = pop_1
    sda_1[1] = exp_1
    sda_1[2] = prp_1
    sda_1[3] = ghg_1

    Millenial[year] = io.sda(sda_1,sda_0)
    
silent_sda = np.zeros((19,4))
boomer_sda = np.zeros((19,4))
gen_x__sda = np.zeros((19,4))
millen_sda = np.zeros((19,4))

for y, year in enumerate(range(2002,2021)):
    silent_sda[y,:] = silent[year][24,1:5]/1000
    boomer_sda[y,:] = babyboomer[year][24,1:5]/1000
    gen_x__sda[y,:] = Gen_x[year][24,1:5]/1000
    millen_sda[y,:] = Millenial[year][24,1:5]/1000  
 

# make SDA variables

temp = np.zeros((1,20))
for y, year in enumerate(years):
    temp[:,y] = np.sum(footdata[year].loc[:,'pop x weight'],0)
    
total_pop = df(temp,columns = years)

pop_prop = pop/np.tile(total_pop,(4,1))

sda_0 = {}
sda_1 = {}
sda2 = {}

for year in range(2002,2021):
    
    ghg_0 = np.transpose(df(ghg_mults_deflated.loc[:,2001].values))
    ghg_1 = np.transpose(df(ghg_mults_deflated.loc[:,year].values))
    
    prp_0 = df(np.transpose(prop_spend_gen[2001]).values)
    prp_1 = df(np.transpose(prop_spend_gen[year]).values)
    
    exp_0 = df(np.diag(per_cap_spend.loc[:,2001]))
    exp_1 = df(np.diag(per_cap_spend.loc[:,year]))
       
    popp_0 = df(np.diag(pop_prop.loc[:,2001]))
    popp_1 = df(np.diag(pop_prop.loc[:,year]))
    
    tpop_0 = df(np.tile(total_pop.loc[:,2001],(4,1)))
    tpop_1 = df(np.tile(total_pop.loc[:,year],(4,1)))


    foot_0 = df.dot(ghg_0,prp_0)
    foot_0 = df.dot(foot_0,exp_0)
    foot_0 = df.dot(foot_0,popp_0)
    foot_0 = df.dot(foot_0,tpop_0)
    
    foot_1 = df.dot(ghg_1,prp_1)
    foot_1 = df.dot(foot_1,exp_1)
    foot_1 = df.dot(foot_1,popp_1)
    foot_1 = df.dot(foot_1,tpop_1)

    sda_0 = {}
    sda_0[0] = ghg_0
    sda_0[1] = prp_0
    sda_0[2] = exp_0
    sda_0[3] = popp_0
    sda_0[4] = tpop_0
    
    sda_1 = {}
    sda_1[0] = ghg_1
    sda_1[1] = prp_1
    sda_1[2] = exp_1
    sda_1[3] = popp_1
    sda_1[4] = tpop_1
    
    sda2[year] = io.sda(sda_1,sda_0)

sda_2 = np.zeros((19,5))
for y, year in enumerate(range(2002,2021)):
   sda_2[y,:] = sda2[year][120,1:6]/1000
