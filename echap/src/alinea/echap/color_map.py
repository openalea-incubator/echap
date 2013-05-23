""" Adaptation of the colormap function from openalea.mtg and openalea.alep to the needs of echap. """
# Imports #########################################################################
from matplotlib.colors import LinearSegmentedColormap

def green_yellow_red(levels=10):
    """ Generate a colormap from green to yellow then red.
    """
    return LinearSegmentedColormap.from_list(name='green_yellow_red', 
                                             colors =[(0., 1., 0.), 
                                                      (1., 1., 0.), 
                                                      (1, 0., 0.)],
                                             N=levels)

def green_lightblue_blue(levels=10):
    """ Generate a colormap from green to light blue then dark blue.
    """
    return LinearSegmentedColormap.from_list(name='green_lightblue_blue', 
                                             colors =[(0., 1., 0.), 
                                                      (0., 1., 1.), 
                                                      (0., 0., 1.)],
                                             N=levels)

