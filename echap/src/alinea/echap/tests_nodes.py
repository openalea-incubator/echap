# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 17:15:32 2013

@author: lepse
"""
from openalea.color import colormap

import numpy as np
from numpy import recfromcsv
import matplotlib.pyplot as plt
from pandas import *
import pylab
from matplotlib import cm

from datetime import datetime, timedelta

from alinea.echap.wheat_mtg import *
from alinea.echap.color_map import *

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


def plot_pesticide_norm(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue, lognorm=False):
    """ plot the plant with normalized pesticide doses """
    prop = g.property(property_name)
    keys = prop.keys()
    value = []
    for k, val in prop.iteritems():
        value.append(val[compound_name])
        v = np.array(value)

    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap())
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap()

    green = (0,180,0)

    norm = Normalize(vmin=0, vmax=max(v)) if not lognorm else LogNorm(vmin=0, vmax=max(v)) 
    values = norm(v)

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    for vid in g.vertices(scale=g.max_scale()): 
        n = g.node(vid)
        if 'surfacic_doses' in n.properties():
            n.properties()['color'] = dict(zip(keys,colors))
        else : 
            n.color = green

    scene = plot3d(g)
    Viewer.display(scene)
    return scene

    
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