#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 16:10:53 2018

This function is used with ukmrio_main to build a metadata object

@author: earao
"""

############################
# used in ukmrio_main_2024 #
############################

def make_meta(ind_sectors,prod_sectors,final_demand,regions): # used in ukmrio_main_2023
    fd_sectors = []
    
    year_1 = list(final_demand.keys())[0]
    
    fd_sectors = final_demand[year_1].columns
    
    fd_sectors = fd_sectors.insert(len(final_demand[year_1].columns),"row final demand")
    
    i_index = []
    for c in range (0,len(regions)):
        for i in range (0,len(ind_sectors)):            
            i_index.append(" ".join([regions[c], ind_sectors[i]]))
       
    p_index = []
    for c in range (0,len(regions)):
        for p in range (0,len(prod_sectors)):            
            p_index.append(" ".join([regions[c], prod_sectors[p]]))
            
    dom_ind_rng = slice(0, len(ind_sectors))   
    imp_ind_rng = slice(len(ind_sectors), len(i_index))
    dom_prd_rng = slice(0, len(prod_sectors))
    imp_prd_rng = slice(len(prod_sectors), len(p_index))
    dom_fd_rng = slice(0, len(fd_sectors)-1)
    imp_fd_rng = slice(len(fd_sectors)-1, len(fd_sectors)+1)
        
    meta = {}
    meta['Z'] = {
        'len_idx': len(p_index)+len(i_index),
        'len_col': len(p_index)+len(i_index),
        }
    meta['bal'] = {
        'len_idx': len(p_index)+len(i_index)+1,
        'len_col': len(p_index)+len(i_index)+len(fd_sectors)
        }
    meta['bal_use'] = {
        'rng_idx': slice(len(i_index),(len(i_index)+len(p_index))),
        'rng_col': slice(0,(len(i_index))),
        }
    meta['bal_sup'] = {
        'rng_idx': slice(0,(len(i_index))),
        'rng_col': slice(len(i_index),(len(i_index)+len(p_index))),
        }
    meta['bal_fd'] = {
        'rng_idx': slice(len(i_index),(len(i_index)+len(p_index))),
        'rng_col': slice(len(i_index)+len(p_index),len(i_index)+len(p_index)+len(fd_sectors)),
        }
    meta['bal_use_dd'] = {
        'rng_idx': slice(len(i_index),(len(i_index)+len(prod_sectors))),
        'rng_col': dom_ind_rng,
        }
    meta['bal_use_id'] = {
        'rng_idx': slice(len(i_index)+len(prod_sectors),len(i_index)+len(p_index)),
        'rng_col': dom_ind_rng,
        }
    meta['bal_use_di'] = {
        'rng_idx': slice(len(i_index),len(i_index)+len(prod_sectors)),
        'rng_col': imp_ind_rng,
        }
    meta['bal_use_ii'] = {
        'rng_idx': slice(len(i_index)+len(prod_sectors),len(i_index)+len(p_index)),
        'rng_col': imp_ind_rng,
        }
    meta['bal_sup_dd'] = {
        'rng_idx': dom_ind_rng,
        'rng_col': slice(len(i_index),(len(i_index)+len(prod_sectors))),
        }
    meta['bal_sup_ii'] = {
        'rng_idx': imp_ind_rng,
        'rng_col': slice(len(i_index)+len(prod_sectors),len(i_index)+len(p_index)),
        }
    meta['bal_fd_dd'] = {
        'rng_idx': slice(len(i_index),(len(i_index)+len(prod_sectors))),
        'rng_col': dom_fd_rng,
        }
    meta['bal_fd_id'] = {
        'rng_idx': slice(len(i_index)+len(prod_sectors),len(i_index)+len(p_index)),
        'rng_col': dom_fd_rng,
        }
    meta['bal_fd_di'] = {
        'rng_idx': slice(len(i_index),len(i_index)+len(prod_sectors)),
        'rng_col': 1,
        }
    meta['bal_fd_ii'] = {
        'rng_idx': slice(len(i_index)+len(prod_sectors),len(i_index)+len(prod_sectors)),
        'rng_col': 1,
        } 
    meta['use'] = {
        'len_idx': len(p_index),
        'len_col': len(i_index),
        'rng': (slice(0, len(p_index)),slice(0, len(i_index))),
        'idx': p_index,
        'col': i_index
        }
    meta['use_dd']= {
        'len_idx': len(prod_sectors),
        'len_col': len(ind_sectors),
        'rng': (dom_prd_rng,dom_ind_rng),
        'rng_col': dom_ind_rng,
        'rng_idx': dom_prd_rng,
        'idx': p_index[dom_prd_rng],
        'col': i_index[dom_ind_rng]
        }
    meta['use_id'] = {
        'len_idx': len(prod_sectors)*(len(regions)-1),
        'len_col': len(ind_sectors),
        'rng': (imp_prd_rng,dom_ind_rng),
        'rng_col': dom_ind_rng,
        'rng_idx': imp_prd_rng,
        'idx': p_index[imp_prd_rng],
        'col': i_index[dom_ind_rng]
        }
    meta['use_di'] = {
        'len_idx': len(prod_sectors),
        'len_col': len(ind_sectors)*(len(regions)-1),
        'rng': (dom_prd_rng,imp_ind_rng),
        'rng_col': imp_ind_rng,
        'rng_idx': dom_prd_rng,
        'idx': p_index[dom_prd_rng],
        'col': i_index[imp_ind_rng]
        }
    meta['use_ii'] = {
        'len_idx': len(prod_sectors)*(len(regions)-1),
        'len_col': len(ind_sectors)*(len(regions)-1),
        'rng': (imp_prd_rng,imp_ind_rng),
        'rng_col': imp_ind_rng,
        'rng_idx': imp_prd_rng,
        'idx': p_index[imp_prd_rng],
        'col': i_index[imp_ind_rng]
        }  
    meta['sup'] = {
        'len_idx': len(i_index),
        'len_col': len(p_index),
        'rng': (slice(0, len(i_index)),slice(0, len(p_index))),
        'idx': i_index,
        'col': p_index
        }
    meta['sup_dd'] =  {
        'len_idx': len(ind_sectors),
        'len_col': len(prod_sectors),
        'rng': (dom_ind_rng,dom_prd_rng),
        'idx': p_index[dom_ind_rng],
        'col': i_index[dom_prd_rng]
      }
    meta['sup_id'] = {
        'len_idx': len(ind_sectors)*(len(regions)-1),
        'len_col': len(prod_sectors),
        'rng': (imp_ind_rng,dom_prd_rng),
        'idx': p_index[imp_ind_rng],
        'col': i_index[dom_prd_rng]
      }
    meta['sup_di'] = {
        'len_idx': len(ind_sectors),
        'len_col': len(prod_sectors)*(len(regions)-1),
        'rng': (dom_ind_rng,imp_prd_rng),
        'idx': p_index[dom_ind_rng],
        'col': i_index[imp_prd_rng]
      }
    meta['sup_ii'] = {
        'len_idx': len(ind_sectors)*(len(regions)-1),
        'len_col': len(prod_sectors)*(len(regions)-1),
        'rng': (imp_ind_rng,imp_prd_rng),
        'idx': p_index[imp_ind_rng],
        'col': i_index[imp_prd_rng]
      }
    meta['fd'] = {
        'len_idx': len(p_index),
        'len_col': len(fd_sectors),
        'rng' : (slice(0, len(p_index)),slice(0, len(fd_sectors))),
        'rng_idx': slice(0, len(p_index)),
        'rng_col': slice(0, len(fd_sectors)),
        'idx': p_index,
        'col': fd_sectors
      }
    meta['fd_dd'] = {
        'len_idx': len(prod_sectors),
        'len_col': len(fd_sectors)-1,
        'rng': (dom_prd_rng,dom_fd_rng),
        'rng_col': dom_fd_rng,
        'rng_idx': dom_prd_rng,
        'idx': p_index[dom_prd_rng],
        'col': fd_sectors[dom_fd_rng]
      }
    meta['fd_id'] = {
        'len_idx': len(prod_sectors)*(len(regions)-1),
        'len_col': len(fd_sectors)-1,
        'rng': (imp_prd_rng,dom_fd_rng),
        'rng_col': dom_fd_rng,
        'rng_idx': imp_prd_rng,
        'idx': p_index[imp_prd_rng],
        'col': fd_sectors[dom_fd_rng]
      }
    meta['fd_di'] = {
        'len_idx': len(prod_sectors),
        'len_col': 1,
        'rng': (dom_prd_rng,imp_fd_rng),
        'idx': p_index[dom_prd_rng],
        'col': fd_sectors[imp_fd_rng]
      }
    meta['fd_ii'] = {
        'len_idx': len(prod_sectors)*(len(regions)-1),
        'len_col': 1,
        'rng': (imp_prd_rng,imp_fd_rng),
        'idx': p_index[imp_prd_rng],
        'col': fd_sectors[imp_fd_rng]
      } 
    meta['v'] = {
        'len': len(i_index),
        'rng': slice(0, len(i_index)),
        'col': i_index
      }
    meta['v_d'] = {
        'len': len(ind_sectors),
        'rng': dom_ind_rng,
        'col': i_index[dom_ind_rng]
      }
    meta['v_i'] = {
        'len': len(ind_sectors)*(len(regions)-1),
        'rng': imp_ind_rng,
        'col': i_index[imp_ind_rng]
      }
    meta['sectors'] = {
      'ind': ind_sectors,
      'prd': prod_sectors
      }
    meta['reg'] = {
      'len': len(regions),
      'rng': slice(0, len(regions)),
      'idx': regions
      }
    
   
    
    return meta

##################
# Not sorted yet #
##################