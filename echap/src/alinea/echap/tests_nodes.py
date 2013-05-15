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

from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe


def plot_pesticide(g, compound_name='Epoxiconazole', colmap=cm.winter_r):
    """ plot the plant with pesticide doses """
    from matplotlib import mpl
    cmap = mpl.cm.get_cmap(colmap)
    green = (0,180,0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'surfacic_doses' in n.properties():
            r,gg,b,s= cmap(n.surfacic_doses[compound_name])
            n.color = (int(r*255),int(gg*255),int(b*255))           
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


def update_meteo_date(timestep=1, t_deb='2000-10-01 01:00:00'):
    date_object = datetime.strptime(t_deb, '%Y-%m-%d %H:%M:%S')
    d = date_object + timedelta(hours=timestep)
    t_deb = datetime.strftime(d, '%Y-%m-%d %H:%M:%S')
    return t_deb


def interface_meteo(meteo_reader, timestep=1, t_deb='2000-10-01 01:00:00'):
    mean_globalclimate, globalclimate, t_deb = meteo_reader.get_meteo_file(timestep, t_deb)
    return mean_globalclimate, globalclimate, t_deb


def sum_temp_global(g, globalclimate):    
    sum_temp = sum(globalclimate['Tair'])
    if not 'sum_temp' in g.property_names():
        for v in g.vertices(scale=1): 
            n = g.node(v)
            n.sum_temp = sum_temp
    else:
        for v in g.vertices(scale=1): 
            n = g.node(v)
            n.sum_temp += n.sum_temp
    return g