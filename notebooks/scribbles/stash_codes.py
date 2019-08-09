# ---Printing stash codes---#

# UKESM ensemble 1, historical data

import iris
file = "bc179a.p41930oct.pp"

# cube list
cl = iris.load(file)
list_of_stash = []

stash_conFIRE = {'vegcover'           : 'm01s03i317',
                 'alpha'              : 'm01s08i223',
                 'emc'                : 'm01s03i245',
                 'treeCover'          : 'm01s03i317',
                 'lightning'          : 'm01s50i082',
                 'pasture'            : 'm01s00i458',
#                'population_density' : 'population_density2000-2014.nc',
                 'relative_humidity'  : 'm01s03i245',
                 'cropland'           : 'm01s00i448'}

for c in cl:
    c = cl[0]
    print(c.name())
    list_of_stash.append(c.attributes["STASH"])
	#print(list_of_stash.head())
