import os
def listdir_path(path):    
    return [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(path)) for f in fn]
    
