# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:57:47 2013

@author: lepse
"""

from alinea.echap.global_meteo import *


def test_mean_meteo():
    t_deb = '2000-10-01 04:00:00'
    t = 0
    dt = 4
    nb_steps = 4
    reader_meteo = Meteo()
    for i in range(nb_steps):
        mean_globalclimate = reader_meteo.get_mean_meteo(t, t_deb, timestep=dt)
        print mean_globalclimate
        t += dt


def test_get_meteo():
    t_deb = '2000-10-01 04:00:00'
    t = 0
    dt = 1
    nb_steps = 4
    reader_meteo = Meteo()
    for i in range(nb_steps):
        globalclimate = reader_meteo.get_meteo_file(t, t_deb, timestep=dt)
        #mean_globalclimate = reader_meteo.get_mean_meteo(t, t_deb, timestep=dt)
        print globalclimate
        #print mean_globalclimate
        t += dt


def test_get_meteoT():
    t_deb = '2000-10-01 04:00:00'
    t = 0
    dt = 4
    nb_steps = 24
    reader_meteo = MeteoT()
    for i in range(nb_steps):
        globalclimate = reader_meteo.get_meteo_file(t_deb, timestep=dt)
        #mean_globalclimate = reader_meteo.get_mean_meteo(t, t_deb, timestep=dt)
        print globalclimate
        #print mean_globalclimate
        t += dt
        date_object = datetime.strptime(t_deb, '%Y-%m-%d %H:%M:%S')
        d = date_object + timedelta(hours=dt)
        t_deb = datetime.strftime(d, '%Y-%m-%d %H:%M:%S') # à mettre dans la fonction!!!!!!!!!
        


