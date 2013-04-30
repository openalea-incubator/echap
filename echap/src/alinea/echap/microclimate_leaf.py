# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.CaribuScene import CaribuScene
import alinea.caribu.sky_tools.turtle as turtle


def microclimate_leaf(sectors, mean_globalclimate, scene):
        """ Calculate the radiation and rain interception with the Caribu interception model and the Tair, humidity and wind for each leaf and stem element id

        :Parameters:
        ----------

        - sectors: 
        - mean_globalclimate: Dict of global climat averaged over an hour defined time step 
        - scene: Scene containing the simulated system
            
        :Returns:
        -------

        - microclimate: Dict of microclimate variables (radiation, rain, Tair, humidity and wind) for each id of leaf and stem element
        """            
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
                microclimate[Infid] = {'radiation': e + a, 'rain': rain_leaf[Infid], 'Tair': mean_globalclimate['Tair'], 'humidity':mean_globalclimate['HR'], 'wind':mean_globalclimate['Vent']} 
    return microclimate


class CaribuMicroclimModel(object):
    """ Adaptor for Caribu model compliying echap local_microclimate model protocol (climate_model)
    """
    def __init__(self, sectors='46'):
        self.sectors = sectors
    def microclim(self, mean_globalclimate, scene):
        """ Return the local meteo calculated with the microclimate_leaf function

        :Parameters:
        ----------

        - mean_globalclimate: Dict of global climat averaged over an hour defined time step 
        - scene: Scene containing the simulated system

        :Returns:
        -------

        - local_meteo: Dict of microclimate variables (radiation, rain, Tair, humidity and wind) for each id of leaf and stem element
        """            
        local_meteo = microclimate_leaf(self.sectors, mean_globalclimate, scene)
        return local_meteo





