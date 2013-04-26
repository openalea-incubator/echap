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
    microclimate = {}
    energy, emission, direction, elevation, azimuth = turtle.turtle(sectors, mean_globalclimate['PAR']) 
    sources = zip(energy, direction)
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    if mean_globalclimate['Pluie'] > 0:
        rain_leaf = c_scene.output_by_id(output, idmap)['Einc']
    EiInf = c_scene.output_by_id(output, idmap)['EiInf']
    EiSup = c_scene.output_by_id(output, idmap)['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                if mean_globalclimate['Pluie'] == 0:
                    microclimate[Infid] = {'radiation': e + a, 'rain': 0} 
                else:
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





