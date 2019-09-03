#!/usr/bin/env python
# coding: utf-8

# ## Retrieving input from UK-ESM: mulitple years
# 
# This is to be run on JASMIN and points to the revelent directories in that working space. For a more detailed walk through, go to retrieve_stash

# In[ ]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import iris.plot as iplt
import warnings
import numpy as np
import pandas as pd

# Not sure these are needed
import sys
sys.path.append('../')

import warnings
warnings.filterwarnings('ignore')

# from   libs.plot_maps    import *


# ### Setting up variables

# In[24]:


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
                 'alpha'              : 'm01s08i223',
                 'lightning'          : 'm01s50i082',
#                'population_density' : 'population_density2000-2014.nc',
                 'relative_humidity'  : 'm01s03i245'}


# In[ ]:


treeCover = [101, 102, 103, 201, 202]
cropland = [301, 401]
pasture = [302, 402]
vegcover = treeCover + cropland + pasture + [3, 4, 501, 502]

name_codes = [treeCover, cropland, pasture, vegcover]
name = ['treeCover', 'cropland', 'pasture', 'vegcover']
for x in range(0, len(name)):
    print(name_codes[x])


# ### Extracting variables from the files
# 
# I'm going to try and condense this down into one cell on Jupter notebooks (it may not work though)

# In[ ]:


stash_l = [ 'lightning', 'relative_humidity']

for l in stash_conFIRE.keys():
    
    # Extracting lightning and relative_humidity
    if l == 'lightning' or l == 'relative_humidity':
        stash_constraint = iris.AttributeConstraint(STASH = stash_conFIRE[l])
        print('Now loading: ' + l)

        # Load all cubes
        aList =[]
        cube_list = iris.cube.CubeList()
        for f in files: 
            dat = iris.load_cube(dir + f, stash_constraint)
            aList.append(dat)
        #    print(str(f) + ' file loaded')

        # Merge all cubes together
        cube_list = iris.cube.CubeList(aList)
        cubes = cube_list.merge_cube()
        
        # Convert units
        if l == 'relative_humidity':
            cubes.convert_units(1)

        # For skipping the first x months
        #xxx
        cubes = cubes[d:,:,:]

        print(l + ' has been saved')
        out = outfile + l + str(years[0]) + '-' + str(years[len(years)-1]) + '.nc'
        iris.save(cubes, out)
        
        
    # For vegcover, treecover, pasture and cropland
    elif l == 'vegcover':
        stash_constraint = iris.AttributeConstraint(STASH = stash_conFIRE[l])
        print('Now loading: ' + l)

        # Load all cubes
        aList =[]
        cube_list = iris.cube.CubeList()
        for f in files: 
            dat = iris.load_cube(dir + f, stash_constraint)
            aList.append(dat)
        #    print(str(f) + ' file loaded')

        # Merge all cubes together
        cube_list = iris.cube.CubeList(aList)
        cube_fractional = cube_list.merge_cube() 


        for var_type in range(0,len(name_codes)):
            index = [cube_fractional.coord('pseudo_level').points == x  for x in name_codes[var_type]]

            # This combines all the boolean arrays together. True + False = True
            index = np.any(index, axis = 0)
            print('Indices for ' + name[var_type])
            #print(index)

            # Extracts just the layers we want and saves
            cube = cube_fractional[:,index]

            # For skipping the first x months
            #xxx
            cube = cube[d:,:,:,:].collapsed(['pseudo_level'], iris.analysis.SUM)

            out = outfile + name[var_type] + str(years[0]) + '-' + str(years[len(years)-1]) + '.nc'
            iris.save(cube, out)
            print(name[var_type] + ' has been saved')
            
            
    # For alpha & alphaMax        
    elif l == 'alpha':
        stash_constraint = iris.AttributeConstraint(STASH = stash_conFIRE[l])
        print('Now loading: ' + l)

        # Load all cubes
        aList =[]
        cube_list = iris.cube.CubeList()
        for f in files: 
            dat = iris.load_cube(dir + f, stash_constraint)
            aList.append(dat)
        #    print(str(f) + ' has loaded')

        # Merge all cubes together
        cube_list = iris.cube.CubeList(aList)
        cube_alpha = cube_list.merge_cube() 

        # Extract just the top soil
        index_soil = [cube_alpha.coord('depth').points == 0.05]
        index_soil = np.any(index_soil, axis = 0) # Still keep this in - it makes the cube happy
        cube_soil = cube_alpha[:, index_soil]
        cube_soil = cube_soil[:,0,:,:]
        cube_soil.long_name = 'alpha'


        # Turning soil moisture into alpha: alpha = soil_moisture * soil_porosity * 1.2 (to scale it) / 50 (convert units)
        porosity = iris.load(dir_poro + 'qrparm.soil.nc')[5] # 5 = soil porosity
        time = len(cube_soil.coord("time").points)
        for t in range(time):
            cube_soil.data[t,:,:] = cube_soil.data[t,:,:] * porosity.data * 1.2 / 50

        # Adjust units to 1
        cube_soil.units = 1 
        
        #xxx
        cube_soil_skip_year = cube_soil[d:,:,:]

        # Save alpha
        out = outfile + cube_soil.long_name + str(years[0]) + '-' + str(years[len(years)-1]) + '.nc'
        iris.save(cube_soil_skip_year, out)
        print(str(l) + ' has been saved')

        # Calculating alphaMax
        #xxx
        cube2 = cube_soil[d:,:,:]
        cube3 = cube_soil[d:,:,:]
        alphaMax = cube_soil[d:,:,:]

        nmonths = len(cube2.coord("time").points)

        #xxx
        for m in range( nmonths):
            cube2.data[m,:,:] = cube_soil[m:m+d,:,:].collapsed(["time"], iris.analysis.MEAN).data
            cube3.data[m,:,:] = cube_soil[m:m+d,:,:].collapsed(["time"], iris.analysis.MAX).data
            alphaMax.data[m,:,:] = (cube3.data[m,:,:] / cube2.data[m,:,:]) - 1


        # Saving alphaMax
        alphaMax.long_name = 'alphaMax'
        out = outfile + alphaMax.long_name + str(years[0]) + '-' + str(years[len(years)-1]) + '.nc'
        iris.save(alphaMax, out)
        print(alphaMax.long_name + ' has been saved')


# In[ ]:




