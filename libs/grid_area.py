import iris

def grid_area(cube):
    if cube.coord('latitude').bounds is None:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    return iris.analysis.cartography.area_weights(cube)
