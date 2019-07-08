import numpy as np
import math

def nanRound(vs, *args, **kw):
    def fun(v):
        if math.isnan(v): return(v)
        return np.around(v, *args, **kw)
    return [fun(i) for i in vs]
