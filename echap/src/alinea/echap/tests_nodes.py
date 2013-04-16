# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 17:15:32 2013

@author: lepse
"""
from openalea.color import colormap

import numpy as np
import matplotlib.pyplot as plt
from pandas import *
import pylab
from matplotlib import cm

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
    return g



