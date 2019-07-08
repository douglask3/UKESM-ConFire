import iris
import numpy as np
from   libs.load_stash      import *
from   libs.plotRegions     import *
from pdb import set_trace as browser

def plotSWoverSW(mod_dir, fign,
                 SWd__code = 'm01s01i210', SWu__code = 'm01s01i211', 
                 data_dir = 'data/', *args, **kw):
    albedo = openSWoverSW(mod_dir, SWd__code, SWu__code, data_dir)
    plotAllRegions(albedo, fign + 'all', *args, **kw)
    return albedo

def openSWoverSW(mod_dir,
                 SWd__code = 'm01s01i210', SWu__code = 'm01s01i211', 
                 data_dir = 'data/'):

    files   = sort(listdir_path(data_dir + mod_dir))
    print(mod_dir)
    if len(files) > 24: files = files[0:120]

    SWd     = load_stash(files, SWd__code, 'SWdown' )
    SWu     = load_stash(files, SWu__code, 'SWup')
    albedo  = SWu/SWd
    
    return albedo
