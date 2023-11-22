#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 17:03:39 2018

@author: earao
"""
import os
import pandas as pd
import numpy as np
df = pd.DataFrame

############################
# used in ukmrio_main_2024 #
############################

def load_io_data(wd, ons_filepath, newyrs, ons_year, ons_name): # used in ukmrio_main_2023

    # load supply, use, final demand, exports data
    supply = {}
    use = {}
    final_demand = {}
    exports = {}
    
    file = os.path.join(ons_filepath, ons_year, ons_name)
    for yr in newyrs:
        temp = pd.read_excel(file, sheet_name=(str(yr)+' Supply'),header=8, index_col=0, usecols="B:DJ",nrows = 112)
        supply[yr] = df(np.transpose(temp.values),index = temp.columns, columns=temp.index)
        use[yr] = pd.read_excel(file, sheet_name=(str(yr)+ ' Use'),header=8, index_col=0, usecols="B:DJ",nrows = 130)
        final_demand[yr] = pd.read_excel(file, sheet_name=(str(yr)+ ' Use'),header=8, index_col=None, usecols="DL:DR",nrows = 113)
        final_demand[yr].index = use[yr].index[0:113]
        exports[yr] = np.sum(pd.read_excel(file, sheet_name=(str(yr)+ ' Use'),header=8, index_col=None, usecols="DS:DV",nrows = 113),1)
        exports[yr].index = use[yr].index[0:113]
        
    # load dom_use, com_use, conc data
    # import lookup needed
    lookup = pd.read_csv(wd + 'data/lookups/io_data_import.csv').set_index(['year', 'data', 'variable', 'conc_type'])\
        .unstack(level=['year', 'data', 'conc_type']).droplevel(axis=1, level=0)
    
    # import dom_use and  com_use data
    dom_use = {}
    com_use = {}
    
    all_data = {}
    for yr in lookup.columns.levels[0].tolist():
        all_data[yr] = {}
        temp = lookup[yr].drop('conc', axis=1).droplevel(axis=1, level=1)
        for dataset in temp.columns.tolist():
            file = ons_filepath + temp.loc['file', dataset]
            header = int(temp.loc['header', dataset])
            index_col = int(temp.loc['index_col', dataset])
            nrows = int(temp.loc['nrows', dataset])
            usecols = temp.loc['usecols', dataset]
            sheet_name = temp.loc['sheet_name', dataset]
            
            all_data[yr][dataset] = pd.read_excel(file, sheet_name=sheet_name, header=header, index_col=index_col, usecols=usecols, nrows=nrows)
            all_data[yr][dataset] = all_data[yr][dataset].apply(lambda x: pd.to_numeric(x, errors='coerce'))
            all_data[yr][dataset].fillna(0, inplace=True)

    # make 1990 com_use from dom_use and imp_use
    all_data[1990]['com_use'] = all_data[1990]['dom_use'].iloc[0:124,:] + all_data[1990]['imp_use'].iloc[0:124,:]
    cols = ['Imports of goods and services', 'Sales by final demand', 'Taxes on expenditure less subsidies', 'Income from employment',
            'Gross profits etc', 'Total inputs']
    
    for i in range(len(cols)):
        item = cols[i]
        all_data[1990]['com_use'].loc[item] = all_data[1990]['dom_use'].iloc[124 + i,:] 
    all_data[1990]['com_use'].fillna(0, inplace=True)
    
    # save in different format
    for yr in lookup.columns.levels[0].tolist():
        com_use[str(yr)] = all_data[yr]['com_use']
        dom_use[str(yr)] = all_data[yr]['dom_use'] 
        
    # import conc data
    conc = {}
    temp = lookup.swaplevel(axis=1, i=0, j=1)['conc']
    file = ons_filepath + 'analytical tables/concordances112.xlsx'
    
    for yr in lookup.columns.levels[0].tolist():
        for item in ['dv', 'dh', 'cv', 'ch']: 
            # import data
            usecols = temp.loc['usecols', (yr, item)]
            conc[str(yr) + '_' + item] = pd.read_excel(file, sheet_name= str(yr) + item[0] + '_' + item[1], header = 1, index_col=0, usecols=usecols)
            
        # fix index
        if yr == 1990:
            conc[str(yr) + '_cv'].index = dom_use[str(yr)].index
        else:
            conc[str(yr) + '_cv'].index = com_use[str(yr)].index
        conc[str(yr) + '_ch'].index = com_use[str(yr)].columns
        conc[str(yr) + '_dv'].index = dom_use[str(yr)].index
        conc[str(yr) + '_dh'].index = dom_use[str(yr)].columns
    
    conc['annxb_v'] = pd.read_excel(file, sheet_name='AnnexB_v',header = 1, index_col=0, usecols="B:EU")
    conc['annxb_h'] = pd.read_excel(file, sheet_name='AnnexB_h',header = 1, index_col=0, usecols="B:EX")
    conc['annxb_v_u'] = pd.read_excel(file, sheet_name='AnnexB_v_u',header = 1, index_col=0, usecols="B:EA")
    conc['annxb_h_u'] = pd.read_excel(file, sheet_name='AnnexB_h_u',header = 1, index_col=0, usecols="B:DW")
    conc['annxb_v_y'] = pd.read_excel(file, sheet_name='AnnexB_v_y',header = 1, index_col=0, usecols="B:DV")
    conc['annxb_s'] = pd.read_excel(file, sheet_name='AnnexB_s',header = 1, index_col=0, usecols="B:DU")
    
    # fix columns
    conc['1995_dh'].columns = conc['annxb_h'].columns
    conc['1995_dv'].columns = conc['annxb_v'].columns
    conc['2005_dh'].columns = conc['annxb_h'].columns
    conc['2005_dv'].columns = conc['annxb_v'].columns
       
    conc['1995_ch'].columns = conc['annxb_h'].columns
    conc['1995_cv'].columns = conc['annxb_v'].columns
    conc['2005_ch'].columns = conc['annxb_h'].columns
    conc['2005_cv'].columns = conc['annxb_v'].columns

    return (supply,use,final_demand,exports,dom_use,com_use,conc)

def load_io_data2(wd, ons_filepath, newyrs, ons_year, ons_name): # used in ukmrio_main_2023

    # load supply, use, final demand, exports data
    supply = {}
    use = {}
    final_demand = {}
    exports = {}
    
    file = os.path.join(ons_filepath, ons_year, ons_name)
    for yr in newyrs:
        temp = pd.read_excel(file, sheet_name=(str(yr)+' Supply'),header=8, index_col=0, usecols="B:DJ",nrows = 112)
        supply[yr] = df(np.transpose(temp.values),index = temp.columns, columns=temp.index)
        use[yr] = pd.read_excel(file, sheet_name=(str(yr)+ ' Use'),header=8, index_col=0, usecols="B:DJ",nrows = 130)
        final_demand[yr] = pd.read_excel(file, sheet_name=(str(yr)+ ' Use'),header=8, index_col=None, usecols="DL:DR",nrows = 113)
        final_demand[yr].index = use[yr].index[0:113]
        exports[yr] = np.sum(pd.read_excel(file, sheet_name=(str(yr)+ ' Use'),header=8, index_col=None, usecols="DS:DV",nrows = 113),1)
        exports[yr].index = use[yr].index[0:113]
    
    # import lookup needed
    lookup = pd.read_csv(wd + 'data/lookups/io_data_import_2024.csv').set_index(['year', 'data', 'variable', 'conc_type'])\
        .unstack(level=['year', 'data', 'conc_type']).droplevel(axis=1, level=0)
    
    # import dom_use, imp_use, transition and com_use
    
    analytic_data = {}
    for yr in lookup.columns.levels[0].tolist():
        analytic_data[yr] = {}
        temp = lookup[yr].drop('conc', axis=1).droplevel(axis=1, level=1)
        for dataset in temp.columns.tolist():
            file = ons_filepath + temp.loc['file', dataset]
            header = int(temp.loc['header', dataset])
            index_col = int(temp.loc['index_col', dataset])
            nrows = int(temp.loc['nrows', dataset])
            usecols = temp.loc['usecols', dataset]
            sheet_name = temp.loc['sheet_name', dataset]
            
            analytic_data[yr][dataset] = pd.read_excel(file, sheet_name=sheet_name, header=header, index_col=index_col, usecols=usecols, nrows=nrows)
            analytic_data[yr][dataset] = analytic_data[yr][dataset].apply(lambda x: pd.to_numeric(x, errors='coerce'))
            analytic_data[yr][dataset].fillna(0, inplace=True)    
        
    # import conc data
    conc = {}
    temp = lookup.swaplevel(axis=1, i=0, j=1)['conc']
    file = ons_filepath + 'analytical tables/concordances112_2024.xlsx'
    
    for yr in lookup.columns.levels[0].tolist():
        if yr < 2015:
            for item in ['dv', 'dh', 'cv', 'ch']: 
                # import data
                usecols = temp.loc['usecols', (yr, item)]
                conc[str(yr) + '_' + item] = pd.read_excel(file, sheet_name= str(yr) + item[0] + '_' + item[1], header = 1, index_col=0, usecols=usecols)
            # fix index
            conc[str(yr) + '_cv'].index = analytic_data[yr]['com_use'].index
            conc[str(yr) + '_ch'].index = analytic_data[yr]['com_use'].columns
            conc[str(yr) + '_dv'].index = analytic_data[yr]['dom_use'].index
            conc[str(yr) + '_dh'].index = analytic_data[yr]['dom_use'].columns
        elif yr == 2015:
            for item in ['dv', 'dh', 'cv', 'ch', 'iv', 'ih', 'tv', 'th']: 
                # import data
                usecols = temp.loc['usecols', (yr, item)]
                conc[str(yr) + '_' + item] = pd.read_excel(file, sheet_name= str(yr) + item[0] + '_' + item[1], header = 1, index_col=0, usecols=usecols)
            # fix index
            conc[str(yr) + '_cv'].index = analytic_data[yr]['com_use'].index
            conc[str(yr) + '_ch'].index = analytic_data[yr]['com_use'].columns
            conc[str(yr) + '_dv'].index = analytic_data[yr]['dom_use'].index
            conc[str(yr) + '_dh'].index = analytic_data[yr]['dom_use'].columns
            conc[str(yr) + '_iv'].index = analytic_data[yr]['imp_use'].index
            conc[str(yr) + '_ih'].index = analytic_data[yr]['imp_use'].columns
            conc[str(yr) + '_tv'].index = analytic_data[yr]['transition'].index
            conc[str(yr) + '_th'].index = analytic_data[yr]['transition'].columns
        else:
            for item in ['dv', 'dh', 'iv', 'ih']: 
                # import data
                usecols = temp.loc['usecols', (yr, item)]
                conc[str(yr) + '_' + item] = pd.read_excel(file, sheet_name= str(yr) + item[0] + '_' + item[1], header = 1, index_col=0, usecols=usecols)
            # fix index
            conc[str(yr) + '_dv'].index = analytic_data[yr]['dom_use'].index
            conc[str(yr) + '_dh'].index = analytic_data[yr]['dom_use'].columns
            conc[str(yr) + '_iv'].index = analytic_data[yr]['imp_use'].index
            conc[str(yr) + '_ih'].index = analytic_data[yr]['imp_use'].columns
        
    conc['annxb_v'] = pd.read_excel(file, sheet_name='AnnexB_v',header = 1, index_col=0, usecols="B:EU")
    conc['annxb_h'] = pd.read_excel(file, sheet_name='AnnexB_h',header = 1, index_col=0, usecols="B:EX")
    conc['annxb_v_u'] = pd.read_excel(file, sheet_name='AnnexB_v_u',header = 1, index_col=0, usecols="B:EA")
    conc['annxb_h_u'] = pd.read_excel(file, sheet_name='AnnexB_h_u',header = 1, index_col=0, usecols="B:DW")
    conc['annxb_v_y'] = pd.read_excel(file, sheet_name='AnnexB_v_y',header = 1, index_col=0, usecols="B:DV")
    conc['annxb_s'] = pd.read_excel(file, sheet_name='AnnexB_s',header = 1, index_col=0, usecols="B:DU")
    
    # fix columns
    conc['1995_dh'].columns = conc['annxb_h'].columns
    conc['1995_dv'].columns = conc['annxb_v'].columns
    conc['2005_dh'].columns = conc['annxb_h'].columns
    conc['2005_dv'].columns = conc['annxb_v'].columns
       
    conc['1995_ch'].columns = conc['annxb_h'].columns
    conc['1995_cv'].columns = conc['annxb_v'].columns
    conc['2005_ch'].columns = conc['annxb_h'].columns
    conc['2005_cv'].columns = conc['annxb_v'].columns
    
    
    return (supply,use,final_demand,exports,analytic_data,conc)

def load_old_io_data(filepath): # used in ukmrio_main_2023
    
    supply = {}
    use = {}
    final_demand = {}
    exports = {}
     
    for yr in range(1992, 1997):
        (supply[yr],use[yr],final_demand[yr],exports[yr]) = manip_old_io_data(os.path.join(filepath, 'ONS supply and use tables/Supply_Use_' + str(yr)[-2:] + '_bb2002.xls'))   

    return (supply,use,final_demand,exports)

def remove_fd_negatives(use,final_demand,yrs,n):  # used in ukmrio_main_2023
    
    for yr in yrs:
        for a in range(0,n):
            for b in range (0,7):
                if final_demand[yr].iloc[a,b] < 0:
                    temp = use[yr].loc['GVA (production measure)'][a]+(final_demand[yr].iloc[a,b]*-1)
                    use[yr].loc['GVA (production measure)',use[yr].columns[a]] = temp
                    final_demand[yr].iloc[a,b] = 0
        final_demand[yr].loc['Total Intermediate consumption'] = np.sum(final_demand[yr].iloc[0:n,:],0)

    return (use,final_demand)

def align_analytic_data(dom_use,com_use,conc): # used in ukmrio_main_2023
    
    for yr in list(dom_use.keys()):
        dom_use[yr] = df.dot(df.dot(df.transpose(conc[yr + '_dv']),dom_use[yr]),conc[yr + '_dh'])
        com_use[yr] = df.dot(df.dot(df.transpose(conc[yr + '_cv']),com_use[yr]),conc[yr + '_ch'])
      
    return (dom_use,com_use)

def align_old_SUT_data(o_supply,o_use,o_final_demand,o_exports,supply,conc): # used in ukmrio_main_2023
    
    conc['annxb_v_y'].columns = o_final_demand[1992].index
    conc['annxb_h_u'].columns = o_use[1992].columns
    conc['annxb_v_u'].columns = o_use[1992].index
    conc['annxb_s'].columns = o_supply[1992].index
    
    supply_prop = df(np.divide(supply.values, np.transpose(np.tile(np.sum(supply,1),(112,1)))))
    
    supply_prop = supply_prop.replace(np.nan, 0)
        
    for yr in list(o_final_demand.keys()):
        o_final_demand[yr].index = o_final_demand[1992].index
        o_exports[yr].index = o_exports[1992].index
        o_supply[yr].index = o_supply[1992].index
    for yr in list(o_use.keys()):
        o_use[yr] =  df.dot(conc['annxb_v_u'],df.dot(o_use[yr],df.transpose(conc['annxb_h_u'])))
        o_final_demand[yr] =  df.dot(conc['annxb_v_y'],o_final_demand[yr])
        o_exports[yr] =  df.dot(conc['annxb_v_y'],o_exports[yr])
        o_supply[yr] = df.dot(conc['annxb_s'],o_supply[yr])
        temp=np.multiply(supply_prop.values,o_supply[yr].values)
        o_supply[yr] = df(temp,index=supply.index,columns=supply.columns)
        o_use[yr].columns = supply.index
    
    o_use[1991] = o_use[1992]*0.968751703
    o_use[1990] = o_use[1991]*0.956573406

    o_supply[1991] = o_supply[1992]*0.968751703
    o_supply[1990] = o_supply[1991]*0.956573406

    o_final_demand[1991] = o_final_demand[1992]*0.968751703
    o_final_demand[1990] = o_final_demand[1991]*0.956573406
    
    o_exports[1991] = o_exports[1992]*0.968751703
    o_exports[1990] = o_exports[1991]*0.956573406
           
    return (o_supply,o_use,o_final_demand,o_exports)

def combine_data(o_supply,o_use,o_final_demand,o_exports,n_supply,n_use,n_final_demand,n_exports,oldyrs,newyrs): # used in ukmrio_main_2023
    
    supply = {}
    use = {}
    final_demand = {}
    exports = {}

    for yr in oldyrs:
        supply[yr] = o_supply[yr]
        use[yr] = o_use[yr]
        final_demand[yr] = o_final_demand[yr]
        exports[yr] = o_exports[yr]
        
    for yr in newyrs:
        supply[yr] = n_supply[yr]
        use[yr] = n_use[yr]
        final_demand[yr] = n_final_demand[yr]
        exports[yr] = n_exports[yr]
        
    
    return(supply,use,final_demand,exports)

def load_hh_data(inputs_filepath,ons_filepath,n_final_demand,newyrs): # used in ukmrio_main_2024
    
    n_final_demand_hh = {}
    file = os.path.join(inputs_filepath, 'COICOP_concs.xlsx')
    tempconc1 = pd.read_excel(file, sheet_name = '103_112_2021',index_col=0)
     
    file = os.path.join(ons_filepath, '2024/supublicationtablesbb23v2.xlsx')
    for yr in newyrs:
        n_final_demand_hh[yr] = pd.read_excel(file, sheet_name=('Table 3 - HHFCe '+str(yr)),header=3, index_col=None, usecols="C:AL",nrows = 104)
        n_final_demand_hh[yr].index = tempconc1.columns
        n_final_demand_hh[yr] = np.divide(n_final_demand_hh[yr],np.transpose(np.tile(np.sum(n_final_demand_hh[yr].values,1),(36,1))))
        n_final_demand_hh[yr].fillna(0, inplace=True)
        n_final_demand_hh[yr]  = df.dot(tempconc1,n_final_demand_hh[yr])
        
        tempElecequip = n_final_demand_hh[yr].loc['Electrical equipment','Furniture, furnishings, carpets etc']+n_final_demand_hh[yr].loc['Electrical equipment','Household appliances']+n_final_demand_hh[yr].loc['Electrical equipment','Audio-visual, photo and info processing equipment']+n_final_demand_hh[yr].loc['Electrical equipment','Other recreational equipment etc']+n_final_demand_hh[yr].loc['Electrical equipment','Miscellaneous goods and services']
        tempMachequip = n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Tools and equipment for house and garden']+n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Other recreational equipment etc']+n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Miscellaneous goods and services']
       
        n_final_demand_hh[yr].loc['Electrical equipment','Furniture, furnishings, carpets etc']=n_final_demand_hh[yr].loc['Electrical equipment','Furniture, furnishings, carpets etc']/tempElecequip
        n_final_demand_hh[yr].loc['Electrical equipment','Household appliances']=n_final_demand_hh[yr].loc['Electrical equipment','Household appliances']/tempElecequip
        n_final_demand_hh[yr].loc['Electrical equipment','Tools and equipment for house and garden']=0
        n_final_demand_hh[yr].loc['Electrical equipment','Audio-visual, photo and info processing equipment']=n_final_demand_hh[yr].loc['Electrical equipment','Audio-visual, photo and info processing equipment']/tempElecequip
        n_final_demand_hh[yr].loc['Electrical equipment','Other recreational equipment etc']= n_final_demand_hh[yr].loc['Electrical equipment','Other recreational equipment etc']/tempElecequip
        n_final_demand_hh[yr].loc['Electrical equipment','Miscellaneous goods and services']=n_final_demand_hh[yr].loc['Electrical equipment','Miscellaneous goods and services']/tempElecequip
       
        n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Furniture, furnishings, carpets etc']=0
        n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Household appliances']=0
        n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Tools and equipment for house and garden']=n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Tools and equipment for house and garden']/tempMachequip
        n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Audio-visual, photo and info processing equipment']=0
        n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Other recreational equipment etc']=n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Other recreational equipment etc']/tempMachequip
        n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Miscellaneous goods and services']=n_final_demand_hh[yr].loc['Machinery and equipment n.e.c.','Miscellaneous goods and services']/tempMachequip
     
        n_final_demand_hh[yr].loc['Alcoholic beverages','Alcoholic beverages']=1
        n_final_demand_hh[yr].loc['Alcoholic beverages','Tobacco']=0
        n_final_demand_hh[yr].loc['Tobacco products','Alcoholic beverages']=0
        n_final_demand_hh[yr].loc['Tobacco products','Tobacco']=1
        
        n_final_demand_hh[yr].loc['Postal and courier services','Postal services']=1
        n_final_demand_hh[yr].loc['Postal and courier services','Restaurants and hotels']=0
        n_final_demand_hh[yr].loc['Accommodation services','Postal services']=0
        n_final_demand_hh[yr].loc['Accommodation services','Restaurants and hotels']=1
        
        temp = np.dot(np.diag(n_final_demand[yr].loc[:,'Households'].values),n_final_demand_hh[yr].values)
        n_final_demand_hh[yr] = df(temp,index = n_final_demand_hh[yr].index,columns = n_final_demand_hh[yr].columns)
  
    return n_final_demand_hh

def make_old_fd_coicop(n_final_demand_hh,Y,oldyrs,meta): # used in ukmrio_main_2023
    
    prop = n_final_demand_hh[1997]/np.transpose(np.tile(np.sum(n_final_demand_hh[1997],1),[36,1]))
    prop.fillna(0,inplace=True)
    
    o_final_demand_hh={}
        
    for yr in oldyrs:
        temp = np.zeros( meta['fd_dd']['len_idx'])
        for r in range(0,meta['reg']['len']):
            temp = temp + Y[yr].loc[:,'Households'].iloc[r* meta['fd_dd']['len_idx']:(r+1)* meta['fd_dd']['len_idx']].values
        
        o_final_demand_hh[yr] = df(np.dot(np.diag(temp),prop.iloc[0:meta['fd_dd']['len_idx'],:]),index = n_final_demand_hh[1997][0:meta['fd_dd']['len_idx']].index,columns = n_final_demand_hh[1997].columns ) 
        o_final_demand_hh[yr].loc['Total Intermediate Demand'] = np.sum(o_final_demand_hh[yr])

    return o_final_demand_hh

def combine_fd_hh(o_final_demand_hh,n_final_demand_hh,oldyrs,newyrs): # used in ukmrio_main_2023

    final_demand_hh = {}

    for yr in oldyrs:
        final_demand_hh[yr] = o_final_demand_hh[yr]
        
    for yr in newyrs:
        final_demand_hh[yr] = n_final_demand_hh[yr]

    return final_demand_hh
 
def make_hh_prop(final_demand_hh,yrs): # used in ukmrio_main_2023
    hh_prop = {}
    for yr in yrs:
        temp = final_demand_hh[yr]/np.transpose(np.tile(np.sum(final_demand_hh[yr],1),(36,1)))
        temp=temp.fillna(0)
        hh_prop[yr] = df(temp, final_demand_hh[yr].index, columns = final_demand_hh[yr].columns)

    return hh_prop
 
def combine_fd(final_demand,final_demand_hh,yrs): # used in ukmrio_main_2023

    col = []
    newfd = {}
    col[0:36] = final_demand_hh[1997].columns
    col[36:36+6] = final_demand[1997].columns[1:7]
    for yr in yrs:
        temp = np.zeros(shape = [final_demand_hh[1997].shape[0],36+6])
        prop = final_demand_hh[yr]/np.transpose(np.tile(np.sum(final_demand_hh[yr],1),(36,1)))
        prop=prop.fillna(0)
        temp[:,0:36] = np.transpose(np.tile(final_demand[yr].iloc[:,0],(36,1)))*prop.values
        temp[:,36:42] = final_demand[yr].iloc[:,1:7]
        newfd[yr] = df(temp, index = final_demand[yr].index, columns = col)
            
    return newfd

def make_wide_Y(Y,hh_prop,meta,yrs): # used in ukmrio_main_2023
    wideY = {}
    col = []
    col[0:36] = hh_prop[yrs[0]].columns
    col[36:36+7] = Y[yrs[0]].columns[1:8]
    
    for yr in yrs:
        temp = np.zeros(shape = [meta['fd']['len_idx'],36+7]) 
        bigprop  = np.tile(hh_prop[yr].iloc[meta['v_d']['rng'],:],[meta['reg']['len'],1])
        temp[:,0:36] = np.multiply(np.transpose(np.tile(Y[yr].iloc[:,0],[36,1])),bigprop)
        temp[:,36:43] = Y[yr].iloc[:,1:8]
        wideY[yr] = df(temp, index = Y[yr].index, columns = col)
    
    return wideY

##################
# Not sorted yet #
##################

def manip_old_io_data(file):
    
    u = pd.read_excel(file, sheet_name='Table 3 int',header=8, index_col=0, usecols="C:DX",nrows = 129)
    u = u.replace('-',0)
    u.rename(index={'Gross value added at basic prices':'GVA (production measure)'}, inplace = True)
    s = pd.read_excel(file, sheet_name='Table 2',header=8, index_col=0, usecols="C:D",nrows = 123)
    s = s.replace('-',0)
    tempy = pd.read_excel(file, sheet_name='Table 3 fd',header=8, index_col=0, usecols="C:Y",nrows = 123)
    y = tempy.iloc[:,0:4]
    y['Gross fixed capital formation'] = tempy.iloc[:,6]
    y['Valuables'] = tempy.iloc[:,7]
    y['Changes in inventoru'] = tempy.iloc[:,8]
    y = y.replace('-',0)
    tempe = pd.read_excel(file, sheet_name='Table 3 fd',header=8, index_col=0, usecols="C:U",nrows = 123)
    e = tempe.iloc[:,17]
    e = e.replace('-',0)
 
    return(s,u,y,e)

