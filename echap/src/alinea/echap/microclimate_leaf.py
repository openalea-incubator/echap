# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.exposition import rain_and_light

# new approach

def microclimate_leaf(g, weather_data, light_sectors='16', domain = None, convUnit = 0.01):
    
    if not ('rain_exposed_area' in g.properties() and 'light_exposed_fraction' in g.properties()):
        g = rain_and_light(g, light_sectors=light_sectors, output_by_triangle = False, domain = domain, convUnit = convUnit, dt = 1)   
    rain_exposed_area = g.property('rain_exposed_area')
    light_exposed_fraction = g.property('light_exposed_fraction')

    for var in ['global_radiation', 'vapor_pressure', 'relative_humidity', 'wind_speed', 'temperature_air', 'PPFD']:
        exec("%s = weather_data[['%s']].mean().item(0)"%(var,var))
    water_amount = weather_data[['rain']].sum().item(0)
    
    microclimate = {}
    
    for vid in light_exposed_fraction:
        microclimate[vid] = {'global_radiation': light_exposed_fraction[vid] * global_radiation,
                             'PPFD': light_exposed_fraction[vid] * PPFD,
                             'water': rain_exposed_area[vid] * water_amount,
                             'relative_humidity': relative_humidity,
                             'wind_speed': wind_speed,
                             'vapor_pressure': vapor_pressure} 
        if microclimate[vid]['PPFD'] == 0:
            microclimate[vid]['temperature_air'] = temperature_air
        else:
            microclimate[vid]['temperature_air'] = temperature_air + microclimate[vid]['PPFD'] / 300.
        
    if not 'microclimate' in g.properties():
        g.add_property('microclimate')
    g.property('microclimate').update(microclimate)    
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





