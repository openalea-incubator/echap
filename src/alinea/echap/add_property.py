""" utilities for creating the differents properties used by echap module (usefull for testing)
"""
import random
from alinea.echap.microclimate_leaf import local_climate


def add_surfacic_doses(g, dose=1):
    g.add_property('surfacic_doses')
    surfacic_doses = g.property('surfacic_doses')
    for vid in g:
        if g.label(vid).startswith('LeafElement'):
            surfacic_doses[vid] = {'defaultCompound':dose}           
    return g
    
def add_microclimate(g, weather_data):
    g.add_property('microclimate')
    g.add_property('rain_star')
    g.add_property('light_star')
    microclimate = g.property('microclimate')
    rain_star = g.property('rain_star')
    light_star = g.property('light_star')
    for vid in g:
        if g.label(vid).startswith('LeafElement'):
            rain_star[vid]=random.random()
            light_star[vid]=random.random()
            microclimate[vid] = local_climate(weather_data, rain_star[vid], light_star[vid])           
    return g
