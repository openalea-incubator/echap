# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.CaribuScene import CaribuScene
import alinea.caribu.sky_tools.turtle as turtle


def climate_from_csv(csvname, delimiter = ';'):
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


def microclimate_leaf(sectors, mean_globalclimate, scene):
    """ 
    energy is the global radiation in W.m-2 (PAR=ppfd in micromol.m-2.sec-1) 0.2174
    """
    microclimate = {}
    energy, emission, direction, elevation, azimuth = turtle.turtle(sectors, energy=mean_globalclimate['Pluie']) 
    sources = zip(energy, direction)
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    rain_leaf = c_scene.output_by_id(output, idmap)['Einc']
    energy, emission, direction, elevation, azimuth = turtle.turtle(sectors, energy=mean_globalclimate['PAR']*0.48) 
    sources = zip(energy, direction)
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    EiInf = c_scene.output_by_id(output, idmap)['EiInf']
    EiSup = c_scene.output_by_id(output, idmap)['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                microclimate[Infid] = {'radiation': e + a, 'rain': rain_leaf[Infid]} 
    return microclimate


class CaribuMicroclimModel(object):
    """ Adaptor for Caribu model compliying echap microclimate_model protocol 
    """
    def __init__(self, sectors='46'):
        self.sectors = sectors
    def microclim(self, mean_globalclimate, scene):
        local_meteo = microclimate_leaf(self.sectors, mean_globalclimate, scene)
        return local_meteo





