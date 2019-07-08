import numpy as np
import iris
from   pdb               import set_trace as browser

import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from   libs.ExtractLocation import *

from scipy.optimize import minimize
from   pdb               import set_trace as browser

class Albedo(object):
    def __init__(self, frac, lai, alpha_0, alpha_inf, k):
        self.frac = frac
        self.lai  = lai
        self.alpha_inf = alpha_inf 
        self.alpha_0 = alpha_0    
        self.k = k 

        self.frac_has_time = len(which_coord(self.frac, 'time')) > 0
    
    def tile(self, tile, alpha_inf = None, k = None):
        
        lai = self.extraPseudo_level(self.lai, tile)
        
        if k         is None: k         = self.k
        if alpha_inf is None: alpha_inf = self.alpha_inf

        #try:
        k         = [i[tile] for i in k        ]
        #except:
        #    browser()
        alpha_inf = [i[tile] for i in alpha_inf]

        def tileMap(alpha_inf, k):
            tile0 = self.lai.coord('pseudo_level').points[0]
            alpha = self.extraPseudo_level(self.lai, tile0).copy()
            alpha.coord('pseudo_level').points[0] = tile
            alpha.attributes = None
        
            if lai is None:            
                alpha.data[:,:,:] = self.alpha_0.data if alpha_inf < 0.0 else alpha_inf   
            else:
                F = alpha.copy()
                F.data[:,:] = 1.0 if k is None else 1.0 - np.exp(-k * lai.data)
                F = F * alpha_inf + (F * (-1) + 1) * self.alpha_0
                alpha.data.data[:] = F.data.data[:]
        
            for i in range(0, alpha.shape[0]):
                alpha.data[i][self.alpha_0.data.mask] = None
            alpha.standard_name = None
            alpha.varn_name = 'albedo'
        
            return alpha
        
        return [tileMap(i, j) for i,j in zip(alpha_inf, k)]
        
    def tiles(self, annual = True, alpha_inf = None, k = None):
        #if alpha_inf is None: 
        alphas =  [self.tile(tile, alpha_inf, k) for tile in self.frac.coord('pseudo_level').points]
        #else:            
        #    alphas =  [self.tile(tile, aInf, k) for tile, aInf, k in zip(self.frac.coord('pseudo_level').points, alpha_inf[0], k[0])]
        #browser()    
        ## switch alphas time and frac dimensions round
        
        alphas = [[i[map] for i in alphas] for map in range(len(alphas[0]))]
        alphas = [iris.cube.CubeList(i).merge_cube() for i in alphas]
        
        if annual: alphas = [i.collapsed('time', iris.analysis.MEAN) for i in alphas]
        return alphas

    def cell(self, annual = False, *arg, **kw):
        alphas = self.tiles(annual, *arg, **kw)
        
        if annual:
            frac = self.frac.collapsed('time', iris.analysis.MEAN)  \
                if self.frac_has_time else self.frac
            alphas = [i * frac for i in alphas]
        else:
            if self.frac_has_time:
                for i in alphas: i.data *= self.frac.data
            else:                
                for tile in range(alphas[0].shape[1]):                   
                    for i in alphas: i.data[:,tile,:] = i.data[:,tile,:] * self.frac.data 

        return [i.collapsed('pseudo_level', iris.analysis.SUM) for i in alphas]


    def extraPseudo_level(self, cube, x):
        return cube.extract(iris.Constraint(pseudo_level = x))


    def Initials(self, x, index = None):
            if index is None: index = self.frac.coord('pseudo_level').points
            return [x[i] for i in index]

    def addBounds(self, cube):
        coords = ('longitude', 'latitude') ## Unhardcode
        try: 
            cube.coord('latitude').guess_bounds()
        except:
            pass
        try: 
            cube.coord('longitude').guess_bounds()
        except:
            pass
        return cube

    def makeParamsForMod(self, params):  
        ln = len(self.antIndex) + len(self.vegIndex)
        nVegP = len(self.alpha_inf) * len(self.vegIndex)
        def Param4Map(mapN):                    
            param = [0 for i in range(ln)]
            sc = mapN * len(self.vegIndex)
                    
            for i in range(len(self.vegIndex)): param[self.vegIndex[i]] = params[sc + i]
            for i in range(len(self.antIndex)): param[self.antIndex[i]] = params[nVegP + i]
            return param
   
        params = [Param4Map(i) for i in range(len(self.alpha_inf))]
        params = [self.Initials(i, self.indexInverse) for i in params]
        params = [dict(zip(self.frac.coord('pseudo_level').points, i)) for i in params]
        return(params)

    def splitOptimizationParams(self, params):

        aInf = params[:self.phalf]
        aK   = params[self.phalf:]
        aK = [None if i < 0 else i for i in aK]            

        print('----')
        print(aInf)
        print(aK)
            
        alpha_inf = self.makeParamsForMod(aInf)
        k         = self.makeParamsForMod(aK  )      
                
        return alpha_inf, k

    def optimize(self, observed, para_index, *args, **kw):
        ##########################
        ## prepare obs          ##
        ##########################
        coords = ('longitude', 'latitude') ## Unhardcode
        self.observed = [self.addBounds(i) for i in observed]

        ##########################
        ## prepare initals      ##
        ##########################        
        para_index0 = para_index
        para_index = self.Initials(para_index)
        indicies   = np.unique(para_index)
        self.indexInverse = [np.where(indicies == i)[0][0] for i in para_index]
        
        vegIndex = self.Initials(self.k[0] , index = indicies)
    
        self.antIndex = np.where([i is     None for i in vegIndex])[0]
        self.vegIndex = np.where([i is not None for i in vegIndex])[0]
    
        def startParams(params):
            start0 = [self.Initials(i , index = indicies) for i in params]
            start  = [[k[i] for i in self.vegIndex] for k in start0]
            start.append([start0[0][i] for i in self.antIndex])
            start  = [i for l in start for i in l]
            start  = [-1.0 if i is None else i for i in start]
            bounds = [(0.01, 1) if i >= 0.0 else (-1.0, -1.0) for i in start]
            
            return start, bounds

        start_k  , bounds_k   = startParams(self.k        )
        start_inf, bounds_inf = startParams(self.alpha_inf)

        start = start_inf + start_k
        bounds = bounds_inf + bounds_k

        self.indexM1 = np.where(np.array(start) == -1.0)[0]
        
        self.phalf = len(start)/2
        
        self.observed = [ExtractLocation(i, *args, **kw).cubes for i in self.observed]
        grid_areas = iris.analysis.cartography.area_weights(self.observed[0])
 
        def minFun(params):
            for i in self.indexM1: params[i] = -1.0


            alpha_inf, k = self.splitOptimizationParams(params)
            modelled = self.cell(annual = False, alpha_inf = alpha_inf, k = k)

            def diffPerMap(mod, obs):
                mod = self.addBounds(mod)
                mod = ExtractLocation(mod, *args, **kw).cubes
                mod.data[i.data < 0.0] = 0.0
                diff = obs.copy()
                diff.data = abs(mod.data - obs.data)
                diff.data = np.nan_to_num(diff.data) 
                
                collapsed_cube = diff.collapsed(coords, iris.analysis.MEAN, weights = grid_areas)
                return collapsed_cube.data.mean()
        
            diff = [diffPerMap(mod, obs) for mod, obs in zip(modelled, self.observed)]
            diff = np.mean(diff)
            print(diff)
            return diff
        
        
        res = minimize(minFun, start, bounds = bounds, method='SLSQP',  options={'xtol': 1e-3, 'disp': True}).x
        #res = [i + 0.1 for i in start]
        
        self.alpha_inf, self.k = self.splitOptimizationParams(res)
        
        #alpha_inf = self.Initials(res[:self.phalf], self.indexInverse)
        #k         = self.Initials(res[self.phalf:], self.indexInverse)
        #self.k   = [None if i < 0 else i for i in self.k]   

        #tileID = self.frac.coord('pseudo_level').points
        #self.alpha_inf = dict(zip( tileID, alpha_inf))
        #self.k         = dict(zip( tileID, k))     

        return self.alpha_inf, self.k


def coord_names(cube):
    return [coord.name() for coord in cube.coords()]

def which_coord(cube, nm):
    return np.where([coord.name() == nm for coord in cube.coords()])[0]


    
