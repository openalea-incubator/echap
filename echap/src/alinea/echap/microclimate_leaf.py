# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.CaribuScene import CaribuScene
import alinea.caribu.sky_tools.turtle as turtle


def microclimate_leaf(sectors, energy, microclimate, rain, scene):
    energy, emission, direction, elevation, azimuth = turtle.turtle(sectors, energy) 
    sources = zip(energy, direction)
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    if rain >0:
        rain_leaf = c_scene.output_by_id(output, idmap)['Einc']
    EiInf = c_scene.output_by_id(output, idmap)['EiInf']
    EiSup = c_scene.output_by_id(output, idmap)['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                if rain == 0:
                    microclimate[Infid] = {'radiation': e + a, 'rain': 0} 
                else:
                    microclimate[Infid] = {'radiation': e + a, 'rain': rain_leaf[Infid]} 
    return microclimate


class CaribuMicroclimModel(object):
    """ Adaptor for Caribu model compliying echap microclimate_model protocol 
    """
    def __init__(self, sectors='46', energy=1):
        self.sectors = sectors
        self.energy = energy
    def microclim(self, microclimate, rain, scene):
        local_meteo = microclimate_leaf(self.sectors, self.energy, microclimate, rain, scene)
        return local_meteo





