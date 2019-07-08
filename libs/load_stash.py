import iris
from   pdb import set_trace as browser
from   pylab import sort 
from   libs.ExtractLocation import *
from   libs.listdir_path import *

def load_stash_dir(dir, *args, **kw):
    files = listdir_path(dir)[0:120]
    return load_stash(files, *args, **kw)

def load_stash(files, code, name = None, units = None):
    print name
    print code
    
    stash_constraint = iris.AttributeConstraint(STASH = code)
    cube = iris.load(files, stash_constraint)[0]

    if name  is not None: cube.var_name = name
    cube.standard_name = None
    if units is not None: cube.units = units
    return cube   

def loadCube(dir, data_dir, code = None, *args, **kw):
    files = sort(listdir_path(data_dir + dir))
    files = files[0:120]
    
    dat = iris.load_cube(files)
     
    dat = ExtractLocation(dat, *args, **kw).cubes

    dat.data = (dat.data > 0.00001) / 1.0
    return dat 
