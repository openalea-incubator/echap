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


def microclimate_leaf(sectors, mean_globalclimate, scene):
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
    id_out = runcaribu(sectors, scene, energy = (((mean_globalclimate['PAR']*0.2174)/0.48)/1000)*3600)
    EiInf = id_out['EiInf']
    EiSup = id_out['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                microclimate[Infid] = {'radiation': e + a, 'rain': rain_leaf[Infid], 'humidity':mean_globalclimate['HR'], 'wind':mean_globalclimate['Vent']} 
                if PAR_leaf[Infid]['PAR'] == 0:
                    microclimate[Infid]['Tair'] = mean_globalclimate['Tair']
                else:
                    microclimate[Infid]['Tair'] = mean_globalclimate['Tair'] + (PAR_leaf[Infid]['PAR'] / 300)
    return microclimate


def parse(yr, doy, hr):
    an, jour, heure = [int(x) for x in [yr, doy, hr/100]]
    dt = datetime(an - 1, 12, 31)
    delta = timedelta(days=jour, hours=heure)
    return dt + delta


class CaribuMicroclimModel(object):
    """ Adaptor for Caribu model compliying echap local_microclimate model protocol (climate_model)
    """
    def __init__(self, sectors='46', data_file = 'E:/openaleapkg/echap/test/meteo01.csv'):
        self.sectors = sectors
        self.data = pd.read_csv(data_file, parse_dates={'datetime':['An','Jour','hhmm']},date_parser=parse,delimiter=';',usecols=['An','Jour','hhmm','PAR','Tair','HR','Vent','Pluie'])
    def microclim(self, scene, timestep, t_deb):
        """ Return the local meteo calculated with the microclimate_leaf function

        :Parameters:
        ----------

        - `scene` - Scene containing the simulated system
        - `timestep` - Timestep of the simulation
        - 't_deb' - Initial start t_deb for the simulation 

        :Returns:
        -------

        - `local_meteo` (Dict) - Dict of microclimate variables (radiation, rain, Tair, humidity, wind) for each id of leaf and stem element
        - `mean_globalclimate` (Dataframe) - Pandas dataframe with microclimate variables averaged for the timestep 
        - `globalclimate` (Dataframe) - Pandas dataframe with microclimate variables for the timestep 
        - `t_deb` (String) - New t_deb after the timestep 
        """            
        df_data = self.data
        index_deb = df_data.index[df_data['datetime']==t_deb][0]
        globalclimate = df_data.truncate(before = index_deb, after = index_deb + (timestep - 1))
        mean_globalclimate = globalclimate.set_index('datetime').mean()
        date_object = datetime.strptime(t_deb, '%Y-%m-%d %H:%M:%S')
        d = date_object + timedelta(hours=timestep)
        t_deb = datetime.strftime(d, '%Y-%m-%d %H:%M:%S')        
        local_meteo = microclimate_leaf(self.sectors, mean_globalclimate, scene)
        return local_meteo, mean_globalclimate, globalclimate, t_deb





