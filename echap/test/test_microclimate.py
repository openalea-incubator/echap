from alinea.adel.astk_interface import AdelWheat
from alinea.astk.plant_interface import *

import pandas

#weather
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather

from alinea.echap.microclimate_leaf import microclimate_leaf

meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')
t_deb = "2000-10-01 01:00:00"
weather = Weather(data_file=meteo_path)
weather.check(['global_radiation','vapor_pressure'])
seq = pandas.date_range(start = "2000-10-02", periods=24, freq='H')

wdata = weather.get_weather(seq)

wheat = AdelWheat()
g,_ = new_canopy(wheat, age=100)
microclimate_leaf(g,wdata)
