# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.weather.global_weather import *
from alinea.astk.caribu_interface import *



class MicroclimateLeaf(object):
    """ Template class of microclimate complying with the interfaces of echap.
    """
    def __init__(self, sectors='16'):
        self.sectors=sectors

    def microclim(self, weather, scene_geometry, timestep, t_deb):
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

        mean_globalclimate, gc = weather.get_weather(timestep, t_deb)
        mean_global_radiation, gc = weather.add_global_radiation(gc)
        mean_vapor_pressure, gc = weather.add_vapor_pressure(gc)

    # Temp
        id_out = turtle_interception(sectors, scene_geometry, energy = mean_globalclimate['PPFD'])
        EiInf = id_out['EiInf']
        EiSup = id_out['EiSup']
        for Infid, e in EiInf.iteritems():
            for Supid, a in EiSup.iteritems():
                if Infid == Supid:
                    PAR_leaf[Infid] = {'PAR': e + a}
    # Rain
        id_out = turtle_interception(sectors_rain, scene_geometry, energy = mean_globalclimate['rain'])
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
                                            'relative_humidity':mean_globalclimate['relative_humidity'],
                                            'wind_speed':mean_globalclimate['wind_speed'],
                                            'vapor_pressure':mean_vapor_pressure} 
                    if PAR_leaf[Infid]['PAR'] == 0:
                        microclimate[Infid]['temperature_air'] = mean_globalclimate['temperature_air']
                    else:
                        microclimate[Infid]['temperature_air'] = mean_globalclimate['temperature_air'] + (PAR_leaf[Infid]['PAR'] / 300)
        return microclimate





