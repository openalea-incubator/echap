# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.CaribuScene import CaribuScene
import alinea.caribu.sky_tools.turtle as turtle
from alinea.weather.global_weather import *

def runcaribu(sectors, scene, energy):
    """ 
    Calls Caribu for differents energy sources

    :Parameters:
    ----------

    - `sectors` (int)
    - `scene` - Scene containing the simulated system
    - `energy` (float) - Meteorological mean variable at the global scale. Could be:
        - 'PAR' : Quantum PAR (ppfd) in micromol.m-2.sec-1
        - 'Pluie' : Precipitation (mm)

    :Returns:
    ----------
        
    - 'id_out' (dict) - Meteorological variable at the leaf scale
    """
    energie, emission, direction, elevation, azimuth = turtle.turtle(sectors=sectors, energy=energy) 
    sources = zip(energie,direction)
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    id_out = c_scene.output_by_id(output, idmap)
    return id_out

class MicroclimateLeaf:
    """ Template class of microclimate complying with the interfaces of echap.
    """
    def __init__(self, sectors='16'):
        self.sectors=sectors

    def microclim(self, weather, scene, timestep, t_deb):
        """ Calculate the radiation and rain interception with the Caribu interception model and the Tair, humidity and wind for each leaf and stem element id

            :Parameters:
            ----------

            - `sectors` (int) 
            - `mean_globalclimate` (Dict) - Dict of global climate averaged over an hourly defined time step
            expected variables are:
                - 'PAR' : Quantum PAR (ppfd) in micromol.m-2.sec-1
                - 'Pluie' : Precipitation (mm)
                - 'Tair' : Temperature of air (Celcius)
                - 'HR': Humidity of air (kPa)
                - 'Vent' : Wind speed (m.s-1)
            - `scene` - Scene containing the simulated system

            :Returns:
            -------

            - `microclimate` (Dict) - Dict of microclimate variables (radiation, rain, Tair, humidity, wind) for each id of leaf and stem element
                - `radiation`: The global radiation in kJ.m-2.h-1
                PAR=ppfd in micromol.m-2.sec-1 (1 ppfd = 0.2174 Watts.m-2.sec-1 PAR, 1 W.m-2.sec-1 global = 0.48 Watts.m-2.sec-1 PAR)
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
        id_out = runcaribu(sectors, scene, energy = mean_globalclimate['PPFD'])
        EiInf = id_out['EiInf']
        EiSup = id_out['EiSup']
        for Infid, e in EiInf.iteritems():
            for Supid, a in EiSup.iteritems():
                if Infid == Supid:
                    PAR_leaf[Infid] = {'PAR': e + a}
    # Rain
        id_out = runcaribu(sectors_rain, scene, energy = mean_globalclimate['rain'])
        rain_leaf = id_out['Einc']
    # PAR
        id_out = runcaribu(sectors, scene, energy = mean_global_radiation)
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





