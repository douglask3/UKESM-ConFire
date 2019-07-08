import numpy as np
import math

import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    try:
        average  = np.average(values, weights = weights)
        variance = np.average((values - average)**2, weights = weights)
        variance = math.sqrt(variance)
    except:
        average = variance = float('nan')
    
    return(average, variance)

def weight_array(ar, wts):
        zipped = zip(ar, wts)
        weighted = []
        for i in zipped:
            for j in range(i[1]):
                weighted.append(i[0])
        return weighted

def weightedBoxplot(data, weights = None, minW = 0.1, *args, **kw):
    def sampleDat(dat, weight):
        weight[weight < minW] = 0.0
        weight = (weight / minW)
        weight = np.around(weight)
        weight = weight.astype(int)
        return(weight_array(dat, weight))
    if (weights is not None):
        data = [sampleDat(d, w) for d, w in zip(data, weights)]

    return plt.boxplot(data, *args, **kw), data
