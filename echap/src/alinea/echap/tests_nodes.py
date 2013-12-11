# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 17:15:32 2013

@author: lepse
"""
from openalea.color import colormap

import numpy
from numpy import recfromcsv
import matplotlib.pyplot as plt
from pandas import *
import pylab
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm

from datetime import datetime, timedelta

from alinea.echap.wheat_mtg import *
from alinea.echap.color_map import *

# Color map
from alinea.echap.color_map import green_lightblue_blue
green_lightblue_blue = green_lightblue_blue(levels=10)

###############################################

def plot_pesticide(g, compound_name='Epoxiconazole', colmap=cm.winter_r):
    """ plot the plant with pesticide doses """
    from matplotlib import mpl
    cmap = mpl.cm.get_cmap(colmap)
    green = (0,180,0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'surfacic_doses' in n.properties():
            r,gg,b,s= cmap(n.surfacic_doses[compound_name]*100)
            n.color = (int(r*255),int(gg*255),int(b*255))           
        else : 
            n.color = green
    scene = plot3d(g)
    Viewer.display(scene)
    return scene


def dose_norm(dose, dose_max_ha):
    """ normalise doses(g.m-2) """
    dn = 0
    if dose_max_ha > 0 :
        dn = float(dose * 1e4) / dose_max_ha
    return dn

def plot_pesticide_norm(g, property_name='surfacic_doses', compound_name='Epoxiconazole', dose_max_ha=125, cmap=green_lightblue_blue):
    """ plot the plant with pesticide doses """
    prop = g.property(property_name)
    keys = prop.keys()
    value = []
    for k, val in prop.iteritems():
        value.append(val[compound_name])
        val = numpy.array(value)
    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap())
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap()
    green = (0,180,0)
    for i in range(0,len(value)):
        val[i] = dose_norm(value[i], dose_max_ha)
    colors = (_cmap(val)[:,0:3])*255
    colors = numpy.array(colors,dtype=numpy.int).tolist()
    for vid in g.vertices(scale=g.max_scale()): 
        n = g.node(vid)
        if 'surfacic_doses' in n.properties():
            n.color = tuple(dict(zip(keys,colors))[vid])
        else : 
            n.color = green
    scene = plot3d(g)
    Viewer.display(scene)
    return g

    
def products_from_csv(csvname, delimiter = ';'):
    """ 
    Read a csv of products parameters and import them in a dict.
    Expected columns are :
        - 'product' : commercial name of the product
        - 'compound' : name of the active compound of the product
        - 'dose' : dose of active compound in the product (g.l-1)    
    """
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = {}
    for i in range(0,len(tab['compound'])):
        d[tab[i][0]] = {tab[i][1]:tab[i][2]}
    return d


def compounds_from_csv(csvname, delimiter = ';'):
    """ 
    Read a csv of compounds parameters and import them in a dict.
    Expected columns are :
        - 'compound'
        - 'dose_max_ha'
        - 'type_code'
        - 'Ap'
        - 'Kp'
        - 'Ae'
        - 'Ke'
        - 'decay_rate'    
    """
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = [dict(zip(tab.dtype.names, data)) for data in tab]
    return d


def sum_temp_global(g, globalclimate):    
    sum_temp = sum(globalclimate['temperature_air'])
    if not 'sum_temp' in g.property_names():
        for v in g.vertices(scale=1): 
            n = g.node(v)
            n.sum_temp = sum_temp
    else:
        for v in g.vertices(scale=1): 
            n = g.node(v)
            n.sum_temp += n.sum_temp
    return g


def generate_scene(g):
    scene = plot3d(g)
    return g, scene


def plot_scene(g):
	scene = plot3d(g)
	Viewer.display(scene)

def wheat_mtg(nb_sect=1):
    g = adel_mtg2(nb_sect=nb_sect)
    return g

def wheat_mtg3(nb_sect=1):
    g = adel_mtg3(nb_sect=nb_sect)
    return g

def update_no_doses(g, label = 'LeafElement'):
    """ Read weather data for a step of simulation and apply it to each leaf.
    """        
    vids = [vid for vid in g if g.label(vid).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.temp = 12
        n.rain_intensity = 0
        n.relative_humidity = 100 
        n.wetness = True
    return g
    
def plot_DU(g):
    """ plot the plant with elements carrying dispersal units in yellow """
    green = (0,180,0)
    yellow = (247, 220, 17)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties() and n.dispersal_units:
            n.color = yellow
        else : 
            n.color = green
    scene = plot3d(g)
    Viewer.display(scene)


def plot_lesions(g):
    """ plot the plant with infected elements in red """
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            n.color = red
        else : 
            n.color = green
    scene = plot3d(g)
    Viewer.display(scene)


def set_initial_properties_g(g, surface_leaf_element=5., position_senescence=None, label = 'LeafElement'):
    """ Give initial values for plant properties of each LeafElement. 
    
    :Parameters:
    ----------
    - 'g': MTG
        MTG representing the canopy
    - 'surface': float
        Initial surface of each leaf element
    - 'position_senescence': float
        Position of senescence on blade axis
    - 'label': str
        Label of the part of the MTG concerned by the calculation
        
    :Returns:
    -------
    - 'g': MTG
        Updated MTG representing the canopy
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.surface = surface_leaf_element
        n.healthy_surface = surface_leaf_element # TODO : Manage properly
        n.position_senescence = position_senescence
    return g