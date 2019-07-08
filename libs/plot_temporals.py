import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from   libs.plot_maps       import *
from   libs.plot_TS         import *
import iris
import numpy as np

def plotInterAnnual(dat, jbName, figN, mnthLength = 30, timeCollapse = iris.analysis.SUM, 
                    *args, **kw):    
    if mnthLength > 1:
        mclim = dat[(mnthLength/2)::mnthLength].copy()

        for mn in range(mclim.shape[0]):
            md = mn * mnthLength
            mclim.data[mn] = dat[md:(md+mnthLength)].collapsed('time', timeCollapse).data
    else:
        mclim = dat
    aclim = mclim[6::12].copy()
    
    for yr in range(aclim.shape[0]):
        ym = yr * 12
        aclim.data[yr] = mclim[ym:(ym+12)].collapsed('time', timeCollapse).data    
        
    labels = [str(i)[10:14] for i in aclim.coord('time')]
    mclim = plotMapsTS(mclim, figN + 'IA', jbName,
                       labels = labels, mdat = aclim, totalMap = timeCollapse,
                       *args, **kw)
    return mclim, aclim, labels

def convert2Climatology(dat, mnthLength = 30,
                        timeCollapse = iris.analysis.SUM, nyrNormalise = True):
    
    yrLength = 12 * mnthLength
    dclim = dat[0:yrLength].copy()
    nyrs = np.floor(dat.shape[0] / yrLength) if nyrNormalise else 1.0
    
    for t in range(0, yrLength):
        dclim.data[t] = dat[t::yrLength].collapsed('time', timeCollapse).data / nyrs
    
    if mnthLength > 1:
        mclim = dat[(mnthLength/2):yrLength:mnthLength].copy()
        for mn in range(0, 12):
            md = mn * 30
            mclim.data[mn] = dclim[md:(md+30)].collapsed('time', timeCollapse).data
    else:
        mclim = dclim
    
    return mclim

def plotClimatology(dat, jbName, figN, mnthLength = 30,
                    timeCollapse = iris.analysis.SUM,  nyrNormalise = True, *args, **kw):
    mclim = convert2Climatology(dat,  mnthLength, timeCollapse, nyrNormalise)
    mclim = plotMapsTS(mclim, figN + 'clim', jbName,
                       *args, **kw)
    return mclim, mclim, 'JFMAMJJASOND'


def plotMapsTS(dat, figN, jbName, levels, cmap = 'pink', labels = 'JFMAMJJASOND',
               units = 'days', mdat = None, nx = 6, ny = 3, running_mean = False,
               tsYlim = None, *args, **kw):
    if mdat is None: mdat = dat
    print(figN)
    try: dat.units = units
    except: pass

    plot_cubes_map(mdat, labels, cmap, levels, nx = 6, ny = 3,
                   cbar_yoff = 0.25, projection = None, *args, **kw)
    
    plt.subplot(4, 1, 4)
    
    tdat = dat if type(dat) == list else [dat]
    
    plot_cube_TS(tdat, running_mean, xticksLabs = labels, ylabel = units, ylim = tsYlim)
    
    plt.title(jbName)  
    
    fig_name = 'figs/' + figN + '-' +jbName + '.png'
    plt.savefig(fig_name)

    try: dat.var_name = dat.long_name = jbName
    except: pass
    return dat
