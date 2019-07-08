import sys
paths = ['/users/global/rjel/.local/lib/python2.7/site-packages/Iris-1.9.1-py2.7-linux-x86_64.egg/', 
         '/users/global/rjel/.local/lib/python2.7/site-packages/cf_units-1.0.0-py2.7.egg', 
         '/users/global/rjel/.local/lib/python2.7/site-packages/Biggus-0.12.0-py2.7.egg', 
         '/users/global/rjel/.local/lib/python2.7/site-packages/pyugrid-0.1.6-py2.7.egg',
         '/users/global/rjel/.local/lib/python2.7/site-packages/Shapely-1.5.13-py2.7-linux-x86_64.egg',
         '/users/global/rjel/.local/lib/python2.7/site-packages/mo_pack-0.2.0-py2.7-linux-x86_64.egg',
         '/users/global/rjel/.local/lib/python2.7/site-packages/Cartopy-0.13.1-py2.7-linux-x86_64.egg',
        '/users/global/rjel/.local/lib/python2.7/site-packages/Iris-1.9.1-py2.7-linux-x86_64.egg',
        '/users/global/rjel/codeAndScripts/python',
        '/users/global/rjel/bin/Jinja/lib/python2.7/site-packages',
        '/users/global/rjel/.local/lib/python2.7/site-packages']

paths.reverse()
for p in paths: sys.path.insert(0, p)
