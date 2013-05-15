# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.CaribuScene import CaribuScene
import alinea.caribu.sky_tools.turtle as turtle
import pandas as pd
from datetime import datetime, timedelta


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


def microclimate_leaf(sectors, meteo_plugin, scene):
    microclimate = {}  
    PAR_leaf = {}
# Temp
    id_out = runcaribu(sectors, scene, energy = mean_globalclimate['PAR'])
    EiInf = id_out['EiInf']
    EiSup = id_out['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                PAR_leaf[Infid] = {'PAR': e + a}
# Rain
    id_out = runcaribu(sectors, scene, energy = mean_globalclimate['Pluie'])
    rain_leaf = id_out['Einc']
# PAR
    id_out = runcaribu(sectors, scene, energy = meteo_plugin.convert_par(mean_globalclimate['PAR']))
    EiInf = id_out['EiInf']
    EiSup = id_out['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                microclimate[Infid] = {'radiation': e + a, 'rain': rain_leaf[Infid], 'humidity':mean_globalclimate['HR'], 'wind':mean_globalclimate['Vent']} 
                if PAR_leaf[Infid]['PAR'] == 0:
                    microclimate[Infid]['Tair'] = mean_globalclimate['Tair']
                else:
                    microclimate[Infid]['Tair'] = meteo_plugin.temp_par(mean_globalclimate['Tair'], PAR_leaf[Infid]['PAR'])
    return microclimate


timestep=3
t_deb='2000-10-01 08:00:00'
meteo_plugin = Meteo()
mean_globalclimate, globalclimate = meteo_plugin.get_meteo_file(timestep, t_deb)
t_deb = meteo_plugin.next_date(timestep, t_deb)
local_meteo = microclimate_leaf(mean_globalclimate, scene)
