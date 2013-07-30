# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.astk.TimeControl import * 
from alinea.weather.global_weather import *
from alinea.astk.caribu_interface import *



class MicroclimateLeaf(object):
    """ Template class of microclimate complying with the interfaces of echap.
    """
    def __init__(self, sectors='16'):
        self.sectors=sectors


    def timing(self, delay, steps, weather, start_date):
        """ compute timing and time_control_sets for a simulation between start and stop. return False when there is no treatment
        """

        def str_to_datetime(t_deb):
            format = "%Y-%m-%d %H:%M:%S"
            if isinstance(t_deb, str):
                t_deb = datetime.strptime(t_deb,format)
            return t_deb

        istart = str_to_datetime(start_date)
        stop = istart + timedelta(hours=steps-1)
        i = istart
        step = []

        mean_globalclimate, gc = weather.get_weather(steps, istart)
        mean_global_radiation,gc = weather.add_global_radiation(gc)
        mean_vapor_pressure,gc = weather.add_vapor_pressure(gc)

        while i <= stop:
            step.append(i)
            i+= timedelta(hours=1)
        event = []
        id = 0
        for i in step:
            if i == str_to_datetime(gc['datetime'][id]):
                event.append([gc['PPFD'][id], gc['rain'][id], gc['global_radiation'][id], gc['vapor_pressure'][id], gc['relative_humidity'][id], gc['wind_speed'][id], gc['temperature_air'][id]])
                id += 1
            else:
                event.append(False)

        pp=[]
        ra=[]
        gr=[]
        vp=[]
        rh=[]
        ws=[]
        ta=[]
        for x in event:
            if x:
                pp.append(x[0])
                ra.append(x[1])
                gr.append(x[2])
                vp.append(x[3])
                rh.append(x[4])
                ws.append(x[5])
                ta.append(x[6])
                duration = len(pp)
            else:
                pp=None
                ra=None
                gr=None
                vp=None
                rh=None
                ws=None
                ta=None
                duration=0
        return (TimeControlSet(PPFD = pp, rain = ra, global_radiation = gr, vapor_pressure = vp, relative_humidity = rh, wind_speed = ws, temperature_air = ta, dt = duration))


    def microclim(self, scene_geometry, time_control):
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
        mean_global_radiation = np.mean(time_control.global_radiation)
        mean_vapor_pressure = np.mean(time_control.vapor_pressure)

    # Temp model with PPFD
        id_out = turtle_interception(sectors, scene_geometry, energy = np.mean(time_control.PPFD))
        EiInf = id_out['EiInf']
        EiSup = id_out['EiSup']
        for Infid, e in EiInf.iteritems():
            for Supid, a in EiSup.iteritems():
                if Infid == Supid:
                    PAR_leaf[Infid] = {'PAR': e + a}
    # Rain
        id_out = turtle_interception(sectors_rain, scene_geometry, energy = np.mean(time_control.rain))
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
                                            'relative_humidity':np.mean(time_control.relative_humidity),
                                            'wind_speed':np.mean(time_control.wind_speed),
                                            'vapor_pressure':mean_vapor_pressure} 
                    if PAR_leaf[Infid]['PAR'] == 0:
                        microclimate[Infid]['temperature_air'] = np.mean(time_control.temperature_air)
                    else:
                        microclimate[Infid]['temperature_air'] = np.mean(time_control.temperature_air) + (PAR_leaf[Infid]['PAR'] / 300)
        return microclimate





