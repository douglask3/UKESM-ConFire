import iris
import numpy as np
from   libs.plot_temporals  import *
from   libs.ExtractLocation import *

regionNames = ['global', 'NA' , 'NA2', 'Asia', 'Asia2', 'Asia3', 'NA3', 'Asia4']
east        = [ None   , 260.0, 275.0,  80.0 ,  80.0  , 80     , 230.0, 0.1]
west        = [ None   , 310.0, 295.0, 105.0 , 105.0  , 105.0  , 310.0, 179.9 ]
south       = [ None   ,  50.0,  45.0,  35.0 ,  45.0  , 30.0   , 20.0,  5.0]
north       = [ None   ,  65.0,  55.0,  50.0 ,  55.0  , 40.0   , 75.0 , 80.0]

def plotRegion(dat, fign, regionName, jobID, levels, units = 'albedo', cmap = "pink",
              nms = None, tsYlim = None, *args, **kw):
    dat = ExtractLocation(dat, *args, **kw).cubes
    figN = fign + '-' + regionName + '-'

#plot_cubes_map(cubes, nms, cmap, levels, extend = 'neither',
#                   figName = None, units = '', nx = None, ny = None, 
##                   cbar_yoff = 0.0, figXscale = 1.0, figYscale = 1.0, 
 #                  totalMap = None, *args, **kw):
    if len(dat.coord('time').points) == 1:
        plot_cubes_map(dat, nms, cmap, levels, figName = 'figs/' + figN, projection = None,
                        figYscale = 0.5)
        return 
    if len(dat.coord('time').points) > 12:
        plotInterAnnual(dat, jobID, figN, mnthLength = 1,
                        timeCollapse = iris.analysis.MEAN,
                        levels = levels, cmap = cmap, units = units, tsYlim = tsYlim)

    plotClimatology(dat, jobID, figN, mnthLength = 1,
                    timeCollapse = iris.analysis.MEAN, nyrNormalise = False,
                    levels = levels, cmap = cmap, units = units, tsYlim = tsYlim)

def plotAllRegions(dat, nm, *args, **kw):
    for r, e, w, s, n in zip(regionNames, east, west, south, north):
        plotRegion(dat, nm, r, east = e, west = w, south = s, north = n, *args, **kw)



