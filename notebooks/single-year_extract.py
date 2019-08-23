# Single-year extraction from historical ensemble 1

###################
##### UK-ESM ######
###################



import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import iris.plot as iplt
import warnings
import numpy as np
import pandas as pd

import sys
sys.path.append('../')

import warnings
warnings.filterwarnings('ignore')

from   libs.plot_maps    import *


# ### Setting up variables


dir = '/home/users/mbrown/UKESM/u-bc179/ap5/'
dir_poro = '/home/users/mbrown/'
outfile = '/home/users/mbrown/outputs/'

# The year ranges that you want
years = range(2000,2015)

months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
files = []
for year in years:
    for month in months:
        files.append('bc179a.p5' + str(year) + month +'.pp')
        
        
d = 12 # 12 # The number of months to skip for alphaMax


# In[5]:


stash_conFIRE = {'vegcover'           : 'm01s03i317',
#                  'alpha'              : 'm01s08i223',
#                  'lightning'          : 'm01s50i082',
#                'population_density' : 'population_density2000-2014.nc',
#                  'relative_humidity'  : 'm01s03i245'}


# In[ ]:

# Removing deciduous tree [101, 201] from treeCover
treeCover = [102, 103, 202]
cropland = [301, 401]
pasture = [302, 402]
vegcover = treeCover + cropland + pasture + [3, 4, 501, 502, 101, 201]

name_codes = [treeCover, cropland, pasture, vegcover]
name = ['treeCover', 'cropland', 'pasture', 'vegcover']
# for x in range(0, len(name)):
#     print(name_codes[x])


# ### Extracting variables from the files
l = 'vegcover' 
var_type ='treeCover'
stash_constraint = iris.AttributeConstraint(STASH = stash_conFIRE[l])
print('Now loading: ' + l)

# Load all 
aList =[]

cube_list = iris.cube.CubeList()
for f in files: 
    dat = iris.load_cube(dir + f, stash_constraint)
    aList.append(dat)

# Merge all cubes together
cube_list = iris.cube.CubeList(aList)
cube_fractional = cube_list.merge_cube() 
                 
index = [cube_fractional.coord('pseudo_level').points == x  for x in name_codes[var_type]]

# This combines all the boolean arrays together. True + False = True
index = np.any(index, axis = 0)
print('Indices for ' + name[var_type])

# Extracts just the layers we want and saves
cube = cube_fractional[:,index]

# For skipping the first x months
#xxx
cube = cube[d:,:,:,:].collapsed(['pseudo_level'], iris.analysis.SUM)

out = outfile + name[var_type] + '_evergreen' + str(years[1]) + '-' + str(years[len(years)-1]) + '.nc'
iris.save(cube, out)
print(name[var_type] + ' has been saved')
        