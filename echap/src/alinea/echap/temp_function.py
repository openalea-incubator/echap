# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:10:31 2013

@author: lepse
"""


def som_temp_micro(g, dt):    
    microclimate = g.property('microclimate')
    temp = g.property('temp')
    temp_leaf={}
    for vid,clim in microclimate.iteritems():
        temp_leaf = dict({vid:clim['Tair']}, **temp_leaf)
    mean_temp_leaf = sum(temp_leaf.values()) / len(temp_leaf.values())
    mean_temp_air = sum(temp.values()) / len(temp.values())
    for v in g.vertices(scale=1): 
        n = g.node(v)
        if not 'sum_temp_leaf' in g.property_names():
            n.sum_temp_leaf = mean_temp_leaf * dt
        else:
            n.sum_temp_leaf += mean_temp_leaf * dt
        if not 'sum_temp_air' in g.property_names():
            n.sum_temp_air = mean_temp_air * dt
        else:
            n.sum_temp_air += mean_temp_air * dt
    return g


def som_temp_global(g, globalclimate):    
    sum_temp = sum(globalclimate['Tair'])
    for v in g.vertices(scale=1): 
        n = g.node(v)
        if not 'sum_temp' in g.property_names():
            n.sum_temp = sum_temp
        else:
            n.sum_temp += sum_temp
    return g
