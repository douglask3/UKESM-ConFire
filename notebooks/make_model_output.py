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
import matplotlib.pyplot as plt
import numpy.ma as ma
get_ipython().run_line_magic('matplotlib', 'inline')
import cartopy.crs as ccrs
from   libs.plot_maps    import *


# ## Input information
# 
# Just needs to say where all the data is stored and fill out a netcdf or pp file for input.


dir = '../data/retrieved_codes/2000-2014/'

files = {'vegcover'           : 'vegcover2000-2014.nc',
         'alphaMax'           : 'alphaMax2000-2014.nc',
         'alpha'              : 'alpha2000-2014.nc',
         'relative_humidity'  : 'relative_humidity_convert_scaled_2000-2014.nc',
         # Modelled:
#          'treeCover'          : 'treeCover2000-2014.nc',
         # Observed:
		 'treeCover'          : 'treecover2000-2014_masked.nc',
         'lightning'          : 'lightning_c2000-2014.nc',
         'pasture'            : 'pasture2000-2014.nc',
         'population_density' : 'pop_dens2000-2014.nc',
         'cropland'           : 'cropland2000-2014.nc'}


# Observation files
# dir = '../data/obs/'
# #
# files = {'alphaMax'           : 'alpha_12monthMax2000-2014_masked.nc',
# 		 'alpha'              : 'alpha2000-2014_masked.nc',
# 		 'cropland'           : 'cropland2000-2014_masked.nc',
# 		 'fire'               : 'fire2000-2014_masked.nc',
# 		 'lightning'          : 'lightning_ignitions2000-2014_masked.nc',
# 		 'pasture'            : 'pasture2000-2014_masked.nc',
# 		 'population_density' : 'population_density2000-2014_masked.nc',
# 		 'relative_humidity'  : 'relative_humidity2000-2014_masked.nc',
# 		 'treeCover'          : 'treecover2000-2014_masked.nc',
# 		 'vegcover'           : 'vegcover2000-2014_masked.nc'}
# 
param_file = '../outputs/params_RH4.csv'
title_output = 'no_tree_scaled_RH_light_c'

fig = True # False
dir_fig = '../figures/2000-2014/' + title_output + '/'

outfile = '../outputs/model_runs/2000-2014/'
File = title_output + '.nc'

# Open data. The model takes data in the same dict class as above.

# In[17]:


##open data
input_data = {}
for key, file in files.items():
    data = iris.load_cube(dir + file)
    input_data[key] = data
    print(key + ': ' + str(input_data[key].units))
    

params = pd.read_csv(param_file)


# Plotting annual averages just to make sure the data looks senible

# In[10]:


nd = 0

plt.figure(figsize = (10, 7.5))

for key, dat in input_data.items():
    nd = nd + 1
    dat = dat.collapsed('time', iris.analysis.MEAN)
    dat.long_name = key
    plot_lonely_cube(dat, 3, 4, nd, cmap = 'magma', levels = None)
    plt.suptitle(title_output, fontsize=16)
if fig:
    plt.savefig(dir_fig + "input_data.png")

# ## The model
# The model is now defined. See documentation paper in NCC for full model equations. This could be moved into a library at some point, but I've defined it here so you can have a proper look.
# 
# The model calculates a number of things needed to predict burnt area, and a few metric (potential limitation and sensitivity) on the fly. The things needed to calculate burnt area are ``ConFIRE.``:
# 
# * **Controls**:
#     * ``fuel``: fuel continuity
#     * ``moisture``: fuel mositure content 
#     * ``ignitions``: potential ignitions
#     * ``suppression``: human fire suppression and landscape fragmentation
#     
# * **Limitation from controls**, the maximum allowed fire considering limitation from:
#     * ``standard_fuel``:  fuel 
#     * ``standard_moisture``:  moisture
#     * ``standard_ignitions``:  ignitions
#     * ``standard_suppression``:  suppression
#     
# * and ``burnt_area``: burnt area from all limitations
# 
# Things calculated on the fly are:
# 
# * **Potential limitation**, the increase in burnt area if limitation where removed from:
#     * ``potential_fuel``:  fuel 
#     * ``potential_moisture``:  moisture
#     * ``potential_ignitions``:  ignitions
#     * ``potential_suppression``:  suppression
# 
# * **Sensitvity**, the rate of change in burnt area for a given control, relative to the maximum possible rate of change for that controls:
#     * ``sensitivity_fuel``:  fuel 
#     * ``sensitivity_moisture``:  moisture
#     * ``sensitivity_ignitions``:  ignitions
#     * ``sensitivity_suppression``:  suppression

# In[11]:


class ConFIRE(object):
    def __init__(self, data, params): 
        """
        Initalise parameters and calculates the key variables needed to calculate burnt area.
        """
        self.params = params
        
        
        ## finds controls
        self.fuel = self.control_fuel(data['vegcover'], data['alphaMax'], self.params['fuel_pw'],
                                      self.params['fuel_pg'])
        
        self.moisture = self.control_moisture(data['alpha'], data['relative_humidity'], data['treeCover'],
                                              self.params['cM'], self.params['cMT'])
        
        self.ignitions = self.control_ignitions(data['lightning'], data['pasture'], data['population_density'],
                                              self.params['cP'], self.params['cD1'])
        
        self.suppression = self.control_suppression(data['cropland'], data['population_density'],
                                              self.params['cD2'])
        
        ## calculates limiting factor of each control.
        self.standard_fuel        = self.sigmoid(self.fuel       ,
                                            self.params[       'fuel_x0'], self.params[       'fuel_k'])
        self.standard_moisture    = self.sigmoid(self.moisture   ,
                                            self.params[   'moisture_x0'], -self.params[   'moisture_k'])
        self.standard_ignitions   = self.sigmoid(self.ignitions   ,
                                            self.params[  'ignitions_x0'], self.params[  'ignitions_k'])
        self.standard_suppression = self.sigmoid(self.suppression,
                                            self.params['suppression_x0'], -self.params['suppression_k'])
        
        ## burnt area us just limitation of each control muliplied together. 
        self.burnt_area = self.params['max_f'] * self.standard_fuel * self.standard_moisture *                           self.standard_ignitions * self.standard_suppression
            
        self.standard_moisture    = self.standard_moisture    /                                     self.sigmoid(0.0, self.params['moisture_x0'], 
                                                 -self.params['moisture_k'])
        self.standard_suppression = self.standard_suppression /                                     self.sigmoid(0.0, self.params['suppression_x0'], 
                                                 -self.params['suppression_k'])
    
        ## if the inputs are iris cubes, we can add some useful metadata
        try:
            self.burnt_area.long_name = "burnt_area"
            
            self.fuel.long_name = "fuel continuity"
            self.fuel.units = '1'
            
            self.moisture.long_name = "moisture content"
            self.moisture.units = '1'
            
            self.ignitions.long_name = "ignitions"
            self.ignitions.units = 'km-2'
            
            self.suppression.long_name = "suppression"
            self.suppression.units = '1'
            
            self.standard_fuel.long_name = "standard_fuel"
            self.standard_moisture.long_name = "standard_moisture"
            self.standard_ignitions.long_name = "standard_ignitions"
            self.standard_suppression.long_name = "standard_suppression"
            
            self.standard_fuel.units = '1'
            self.standard_moisture.units = '1'
            self.standard_ignitions.units = '1'
            self.standard_suppression.units = '1'
        except:
            pass
            

    def control_fuel(self, vegcover, alphaMax, fuel_pw, fuel_pg):
        """
        Definition to describe fuel load: while return the input; capability to be modified later.
        """        
        return (vegcover**fuel_pw) * (fuel_pg * (alphaMax-1) + 1) / (1 + fuel_pg)

    
    def control_moisture(self, alpha, emc, treeCover, cM, cMT):
        """
        Definition to describe moisture
        """
        return (alpha + cM*emc + cMT * treeCover) / (1 + cM + cMT)

    
    def control_ignitions(self, lightning, pasture, population_density, cP, cD1):
        """
        Definition for the measure of ignition
        """
        ignite = lightning + cP*pasture + cD1*population_density

        return ignite

    
    def control_suppression(self, cropland, population_density, cD2):
        """
        Definition for the measure of fire supression
        """
        return cropland + cD2*population_density
    
        """
        Defines potential limitation for each control in turn
        """
    def potential_fuel(self):
        return self.potential(self.standard_fuel, "potential_fuel")
    
    
    def potential_moisture(self):
        return self.potential(self.standard_moisture, "potential_moisture")
    
    
    def potential_ignitions(self):
        return self.potential(self.standard_ignitions, "potential_ignitions")
    
    
    def potential_suppression(self):
        return self.potential(self.standard_suppression, "potential_suppression")
    
    
    def sensitivity_fuel(self):
        return self.sensitivity(self.fuel, self.params['fuel_x0'], self.params['fuel_k'],
                                self.standard_fuel, "sensitivity_fuel")
    
    """
    Defines sensitivity for each control in turn
    """
    def sensitivity_moisture(self):
        return self.sensitivity(self.moisture, self.params['moisture_x0'], -self.params['moisture_k'], 
                                self.standard_moisture, "sensitivity_moisture")
    
    
    def sensitivity_ignitions(self):
        return self.sensitivity(self.ignitions, self.params['ignitions_x0'], self.params['ignitions_k'], 
                                self.standard_ignitions, "sensitivity_ignitions")
    
    
    def sensitivity_suppression(self):
        return self.sensitivity(self.suppression, self.params['suppression_x0'], -self.params['suppression_k'] , 
                                self.standard_suppression, "sensitivity_suppression")
    
          
    def sensitivity(self, x, x0, k, fi, long_name = None):  
        
        gradient = self.gradient(x, x0, k)        
        sens = gradient * self.control_removal(fi) 
        
        try: sens.units = '1'
        except: pass 
        
        if long_name is not None: 
            try: sens.long_name = long_name
            except: pass
        return sens
    
    
    def control_removal(self, fi):
        return self.burnt_area/fi
    
    
    def potential(self, fi, long_name = None):
        out = fi.copy()
        out.data = self.burnt_area.data * ((1/out.data) - 1)
        
        try: out.units = '1'
        except: pass
        
        if long_name is not None: 
            try: out.long_name = long_name
            except: pass
        return out
    
    
    def gradient(self, x, x0, k, dx = 0.0001):

        f1 = self.sigmoid(x + dx, x0, k)
        f2 = self.sigmoid(x - dx, x0, k)
                
        n1 = self.sigmoid(x0 + dx, x0, k)
        n2 = self.sigmoid(x0 - dx, x0, k)
        
        try:
            f1 = (f1 - f2)/(n1 - n2)
        except:
            f1.data = (f1.data - f2.data) / (n1.data - n2.data)
        
        return f1
    
    
    def sigmoid(self, x, x0,k):
        """
        Sigmoid function to describe limitation using tensor
        """        
        try:
            out = x.copy()
            out.data = -k*(x.data - x0)
            out.data = 1.0/(1.0 + np.exp(out.data))
            x = out
        except:
            x = -k * (x - x0)
            x = 1.0/(1.0 + np.exp(x))
        return x


# ## Check that everything is working okay
# Here, we run the model with the median of each parameter and plot each of the outputs mentioned abpver to make sure everything is happy

# In[12]:

model = ConFIRE(input_data, params.loc[params["sigma"].idxmin()]) # We're using the row which has minimum sigma


# ### Plotting
# #### Burnt area

# In[14]:


burnt_area = model.burnt_area.collapsed('time', iris.analysis.MEAN)
burnt_area.long_name = "Annual burnt area (%)"
burnt_area.data = burnt_area.data * 1200
print(type(burnt_area))
plot_lonely_cube(burnt_area, levels = [0, 1, 2, 5, 10, 20, 50, 100], cmap = "brewer_YlOrRd_09")
plt.suptitle(title_output, fontsize=16)
if fig:
    plt.savefig(dir_fig + 'burnt_area.png')


# #### Compare with 'observed' burnt area

# In[17]:


# plt.rcParams['figure.figsize'] = [12, 6]
# obs_BA = iris.load_cube(dir + 'fire2000-2014_masked.nc')
# obs_BA.long_name = " 'Observed' burnt area (%)"

# dat = obs_BA.collapsed('time', iris.analysis.MEAN)
# dat.data = dat.data * 1200 # To make annual and a percentage
# plot_lonely_cube(dat, 1, 2, 1, cmap = 'brewer_YlOrRd_09', levels = [0, 1, 2, 5, 10, 20, 50, 100])


# #### Controls

# In[18]:


def plotModComponet(comp, n, name = None, levels = None, scale = None, cmap = "brewer_YlOrRd_09", *args, **kws):
    comp = comp.collapsed('time', iris.analysis.MEAN)
    if scale is not None: comp.data = comp.data * scale
    if name  is not None: comp.long_name = name
    plot_lonely_cube(comp, 2, 2, n, levels = levels, cmap = cmap, *args, **kws)

cmap_fuel = 'brewer_Greens_09' 
cmap_moisture = 'brewer_PuBu_09'
cmap_ignitions = 'brewer_Reds_09'
cmap_suppression = 'brewer_Greys_09'
    
plt.figure(figsize = (10, 7.5))   
plotModComponet(model.fuel, 1, levels = [0, 0.2, 0.4, 0.6, 0.8, 1.0], cmap = cmap_fuel)
plotModComponet(model.moisture, 2, cmap = cmap_moisture)
plotModComponet(model.ignitions, 3, cmap = cmap_ignitions)
plotModComponet(model.suppression, 4, cmap = cmap_suppression)
plt.suptitle(title_output, fontsize=16)
if fig:
    plt.savefig(dir_fig + 'controls.png')

# #### Standard Limitation

# In[19]:


plt.figure(figsize = (10, 7.5))   
plotModComponet(model.standard_fuel, 1, cmap = cmap_fuel)
plotModComponet(model.standard_moisture, 2, cmap = cmap_moisture)
plotModComponet(model.standard_ignitions, 3, cmap = cmap_ignitions)
plotModComponet(model.standard_suppression, 4, cmap = cmap_suppression)
plt.suptitle(title_output, fontsize=16)
if fig:
    plt.savefig(dir_fig + 'standard_limitation.png')

# #### Potential limitation

# In[20]:


plt.figure(figsize = (10, 7.5))  
levels = [0, 0.1, 1, 2, 5, 10, 20, 50, 100]
plotModComponet(model.potential_fuel(), 1, levels = levels, scale = 100,
                cmap = cmap_fuel)
plotModComponet(model.potential_moisture(), 2, levels = levels, scale = 100,
                cmap = cmap_moisture)
plotModComponet(model.potential_ignitions(), 3, levels = levels, scale = 100,
                cmap = cmap_ignitions)
plotModComponet(model.potential_suppression(), 4, levels = levels, scale = 100,
                cmap = cmap_suppression)
plt.suptitle(title_output, fontsize=16)
if fig:
    plt.savefig(dir_fig + 'potential_limitation.png')

# #### Sensitivty

# In[21]:


plt.figure(figsize = (10, 7.5))  
levels = np.array([0, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10])
plotModComponet(model.sensitivity_fuel(), 1, scale = 100, levels = levels, 
                cmap = cmap_fuel, extend = 'max')
plotModComponet(model.sensitivity_moisture(), 2, scale = 100, levels = levels, 
                cmap = cmap_moisture, extend = 'max')
plotModComponet(model.sensitivity_ignitions(), 3, scale = 100, levels = levels, 
                cmap = cmap_ignitions, extend = 'max')
plotModComponet(model.sensitivity_suppression(), 4, scale = 100, levels = levels, 
                cmap = "Greys", extend = 'max')
plt.suptitle(title_output, fontsize=16)
if fig:
    plt.savefig(dir_fig + 'sensitivity.png')

# In[22]:
cubes = [model.burnt_area]
cubes = cubes + [model.fuel, model.moisture,
                 model.ignitions, model.suppression]

cubes = cubes + [model.standard_fuel, model.standard_moisture,
                 model.standard_ignitions, model.standard_suppression]

cubes = cubes + [model.potential_fuel(), model.potential_moisture(),
                 model.potential_ignitions(), model.potential_suppression()]

cubes = cubes + [model.sensitivity_fuel(), model.sensitivity_moisture(),
                 model.sensitivity_ignitions(), model.sensitivity_suppression()]

cubes = iris.cube.CubeList(cubes)
iris.save(cubes, outfile + File)

# In[ ]:




