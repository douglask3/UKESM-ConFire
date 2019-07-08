import iris
import numpy as np
import cartopy.crs as ccrs

import iris.plot as iplt
import matplotlib.pyplot as plt
from pdb import set_trace as browser


def grid_area(cube):
    if cube.coord('latitude').bounds is None:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)    
    #grid_areas[cube.data.mask] = 0.0
    return grid_areas

### Running mean/Moving average
def running_N_mean(l, N):
    sum = 0
    result = list( 0 for x in l)
    
    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = float('nan')#sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result

def cube_TS(cube, running_mean = False):
    grid_areas = grid_area(cube) 
    #cube = cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    
    cube = cube.collapsed(['longitude'], iris.analysis.MEAN).\
                collapsed(['latitude' ], iris.analysis.MEAN, weights = grid_areas[:,:,0])
    
    if (running_mean): cube.data = running_N_mean(cube.data, 12)
    return cube   

def plot_cube_TS(cubes, running_mean, xticksLabs = None, ylabel = '', ylim = None):    
    cubes = [cube_TS(cube, running_mean) for cube in cubes]    
    
    for cube in cubes:
        label = cube.name() if cube.var_name is None else cube.var_name
        plt.plot(cubes[0].data)
    
    tickMarks = range(0, len(cube.data), len(cube.data) / len(xticksLabs))
    plt.xticks(tickMarks, xticksLabs)
    if len(cube.data) == 12: plt.xlim([0,11])
    plt.ylabel(ylabel)
    plt.legend(ncol = 2, loc = 0)
    plt.grid(True)    
    plt.axis('tight')
    if ylim is not None: plt.gca().set_ylim([-0.25,0.05])
    plt.gca().set_ylabel(cubes[0].units, fontsize=16)
