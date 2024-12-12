# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
df = pd.DataFrame
import pickle
from sys import platform
import time

start_time = time.time()

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
oac_none_years = np.array([int(x) for x in range(2001,2007)])
oac_2001_years = np.array([int(x) for x in range(2007,2014)])
oac_2011_years = np.array([int(x) for x in range(2014,2022)])

S = pickle.load( open(results_filepath + "S.p", "rb" ) )
U = pickle.load( open(results_filepath + "U.p", "rb" ) )
Y = pickle.load( open(results_filepath + "Y.p", "rb" ) )

U2 = np.dot(U[2020],U[2020])

end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)

