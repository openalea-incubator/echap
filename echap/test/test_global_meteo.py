# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:57:47 2013

@author: lepse
"""

from alinea.echap.global_meteo import *
from datetime import datetime, timedelta

def update_meteo_date(t_deb, dt):
    date_object = datetime.strptime(t_deb, '%Y-%m-%d %H:%M:%S')
    d = date_object + timedelta(hours=dt)
    t_deb = datetime.strftime(d, '%Y-%m-%d %H:%M:%S')
    return t_deb


def test_get_meteo():
    t_deb = '2000-10-01 04:00:00'
    t = 0
    dt = 4
    nb_steps = 24
    reader_meteo = Meteo()
    for i in range(nb_steps):
        mean_globalclimate, globalclimate, t_deb = reader_meteo.get_meteo_file(timestep=dt, t_deb=t_deb)
        print t_deb
        t += dt






