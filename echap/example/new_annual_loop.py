import pandas

#from openalea.core.alea import load_package_manager, function
from openalea.deploy.shared_data import shared_data

from alinea.astk.Weather import Weather
from alinea.astk.TimeControl import *
from alinea.caribu.caribu_star import rain_and_light_star
from alinea.alep.protocol import infect
from alinea.echap.interception_leaf import pesticide_applications
from alinea.echap.microclimate_leaf import microclimate_leaf
from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record #to avoid confusion with numpy record


# pm = load_package_manager()
# node_factory = pm['alinea.echap.macros']['setup canopy']
# setup_canopy = function(node_factory)
#node_factory = pm['alinea.echap.macros']['update lesions']
#update_lesions = function(node_factory)
#node_factory = pm['alinea.echap.macros']['update pesticide']
#update_pesticides = function(node_factory)
#node_factory = pm['alinea.echap.macros']['dispersion']
#dispersion = function(node_factory)
#node_factory = pm['alinea.echap.macros']['contamination']
#contamination = function(node_factory)

from macros_annual_loop import setup_canopy, update_lesions, update_pesticides, dispersion, contamination, pesticide_intercept

from alinea.echap.tests_nodes import plot_pesticide

import alinea.septo3d
meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')

#setup
weather = Weather(data_file=meteo_path)
weather.check(['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
applications = """date,dose, product_name
2000-11-20 00:00:00, 0, bug
2000-11-20 1:00:00, 1, Opus
2000-11-20 2:00:00, 0, bug
2001-04-13 10:00:00, 1, Opus
2001-04-13 11:00:00, 0, bug
2001-05-13 10:00:00, 1, Opus
2001-05-13 11:00:00, 0, bug
"""
pest_calendar = pesticide_applications(applications)

seq = pandas.date_range(start = "2000-11-20", periods=100, freq='H')  # 5000 pour le cycle complet
g, adel, domain, domain_area, convUnit, nplants = setup_canopy()
SspoSol = 0.01
recorder = LeafElementRecorder()

TTmodel = DegreeDayModel(Tbase = 0)
every_rain = rain_filter(seq, weather)
every_h = time_filter(seq, delay=6)
every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 15)
every_pest = date_filter(seq, pest_calendar)

rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
canopy_timing = IterWithDelays(*time_control(seq, every_dd, weather.data))
doses_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
pest_timing = IterWithDelays(*time_control(seq, every_pest, pest_calendar))

for i,controls in enumerate(zip(canopy_timing, doses_timing, rain_timing, pest_timing)):
    canopy_iter, doses_iter, rain_iter, pest_iter = controls
    if canopy_iter:
        print canopy_iter.value.index[-1]
        g = adel.grow(g, canopy_iter.value)
        _=rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
        _=update_lesions(g, canopy_iter.dt, True)   
        _=do_record(g, canopy_iter.value, recorder)
    if pest_iter:
        _=pesticide_intercept(g, pest_iter.value)
        #plot_pesticide(g)
    if doses_iter:
        print 'updte microclimate / doses...'
        _=microclimate_leaf(g, doses_iter.value, domain = domain, convUnit = convUnit)
        _=update_pesticides(g, doses_iter.value)
        plot_pesticide(g)
        _=do_record(g, doses_iter.value, recorder, header={'iter':i, 'TT':adel.canopy_age})
    if rain_iter:
        wdata = rain_iter.value
        dispersion(g, wdata, domain, domain_area, convUnit)
        _=contamination(g,wdata, SspoSol, domain, domain_area, convUnit)
        _=infect(g, rain_iter.dt)
        
recorder.save_records('test.csv')     