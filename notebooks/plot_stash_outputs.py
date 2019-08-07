#!/usr/bin/env python
# coding: utf-8

# ## Plotting outputs from retrieved stash codes

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

import sys
sys.path.append('../')

import warnings
warnings.filterwarnings('ignore')

import os
from   io     import StringIO
import numpy  as np
import pandas as pd
import csv

import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy.ma as ma
get_ipython().run_line_magic('matplotlib', 'inline')
import cartopy.crs as ccrs
from   libs.plot_maps    import *
from netCDF4 import Dataset


# In[2]:


dir = '/home/users/mbrown/outputs/'
files = {'vegcover'           : 'vegcover2000-2014.nc',
         'alphaMax'           : 'alphaMax2000-2014.nc',
         'alpha'              : 'alpha2000-2014.nc',
         'relative_humidity'  : 'relative_humidity2000-2014.nc',
         'treeCover'          : 'treecover2000-2014.nc',
         'lightning'          : 'lightning2000-2014.nc',
         'pasture'            : 'pasture2000-2014.nc',
 #        'population_density' : 'population_density2000-2014.nc',
         'cropland'           : 'cropland2000-2014.nc'}


# In[ ]:


input_data = {}
for key, file in files.items(): 
    data = iris.load_cube(dir + file)
    input_data[key] = data


# In[ ]:


nd = 0

plt.figure(figsize = (10, 7.5))

for key, dat in input_data.items():
    nd = nd + 1
    dat = dat.collapsed('time', iris.analysis.MEAN)
    dat.long_name = key
    plot_lonely_cube(dat, 3, 4, nd, cmap = 'magma', levels = None)    


# In[ ]:


# Or use this to plot


# In[ ]:


# iplt.contourf(input_data['vegcover'][0,:,:], cmap = 'magma')
# iris.quickplot.pcolormesh(input_data['vegcover'][0,:,:], cmap = 'magma')
# plt.title('Vegcover: full resolution')

