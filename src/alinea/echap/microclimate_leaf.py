# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.caribu_star import rain_and_light_star

# new approach

def local_climate(global_climate, rain_star = 0.5, light_star = 0.5, light_cols = ['global_radiation', 'PPFD'], rain_cols = ['rain'], temp_cols = ['temperature_air']):
    """ compute a simplistic microclimate"""
    local = global_climate.copy()
    for col in light_cols:
        if col in local.columns:
            local[col] *= light_star
    for col in rain_cols:
        if col in local.columns:
            local[col] *= rain_star
    for col in temp_cols:
        if col in local.columns and 'PPFD' in local.columns:
            local[col] = local[col] + local['PPFD'] / 300
    return local


def microclimate_leaf(g, weather_data, light_sectors='16', domain = None, convUnit = 0.01, label='LeafElement'):
    
    if not ('rain_star' in g.properties() and 'light_star' in g.properties()):
        g = rain_and_light_star(g, light_sectors=light_sectors, output_by_triangle = False, domain = domain, convUnit = convUnit, trigger = 1)   
    rain_star = g.property('rain_star')
    light_star = g.property('light_star')
    
    if not 'microclimate' in g.properties():
        g.add_property('microclimate')
    microclimate = g.property('microclimate')
      
    for vid in g:
        if g.label(vid).startswith(label):
            if vid in light_star: # may not be here because of asynchrony of leaf formation and with caribulightstar
                if vid in rain_star:
                    microclimate[vid] = local_climate(weather_data, rain_star[vid], light_star[vid]) 
            
    return g

# deprecated
from alinea.astk.caribu_interface import *


class MicroclimateLeaf(object):
    """ Template class of microclimate complying with the interfaces of echap.
    """
    def __init__(self, sectors='16'):
        self.sectors=sectors


    def microclim(self, scene_geometry, weather_data):
        """ 
        Calculate the radiation and rain interception with the Caribu interception model and the Tair, humidity and wind for each leaf and stem element id

        :Parameters:
        ----------
        - `sectors` (int) 
        - `weather` (class)
            A class embending the weather reader
        - `scene_geometry`
        - 'timestep' (int)
            The timestep of the simulation
        - 't_deb' (datetime)

        :Returns:
        -------
        - `microclimate` (Dict)
            Dict of microclimate variables (radiation, rain, Tair, humidity, wind) for each id of leaf and stem element
                - `radiation`: The global radiation in kJ.m-2.h-1
                PAR = PPFD in micromol.m-2.sec-1 (1 PPFD = 0.2174 Watts.m-2.sec-1 PAR, 1 W.m-2.sec-1 global = 0.48 Watts.m-2.sec-1 PAR)
                - `Tair`: Temperature of air near the leaf (Celcius)
        """
        sectors = self.sectors
        sectors_rain = '1'

        microclimate = {}  
        PAR_leaf = {}

        mean_global_radiation = weather_data[['global_radiation']].mean().item(0)
        mean_vapor_pressure = weather_data[['vapor_pressure']].mean().item(0)

    # Temp model with PPFD
        id_out = turtle_interception(sectors, scene_geometry, energy = weather_data[['PPFD']].mean().item(0))
        EiInf = id_out['EiInf']
        EiSup = id_out['EiSup']
        for Infid, e in EiInf.iteritems():
            for Supid, a in EiSup.iteritems():
                if Infid == Supid:
                    PAR_leaf[Infid] = {'PAR': e + a}
    # Rain
        id_out = turtle_interception(sectors_rain, scene_geometry, energy = weather_data[['rain']].mean().item(0))
        rain_leaf = id_out['Einc']
    # PAR
        id_out = turtle_interception(sectors, scene_geometry, energy = mean_global_radiation)
        EiInf = id_out['EiInf']
        EiSup = id_out['EiSup']
        for Infid, e in EiInf.iteritems():
            for Supid, a in EiSup.iteritems():
                if Infid == Supid:
                    microclimate[Infid] = {'global_radiation': e + a, 
                                            'rain': rain_leaf[Infid],
                                            'relative_humidity':weather_data[['relative_humidity']].mean().item(0),
                                            'wind_speed':weather_data[['wind_speed']].mean().item(0),
                                            'vapor_pressure':mean_vapor_pressure} 
                    if PAR_leaf[Infid]['PAR'] == 0:
                        microclimate[Infid]['temperature_air'] = weather_data[['temperature_air']].mean().item(0)
                    else:
                        microclimate[Infid]['temperature_air'] = weather_data[['temperature_air']].mean().item(0) + (PAR_leaf[Infid]['PAR'] / 300)
        return microclimate





