import pandas
from openalea.deploy.shared_data import shared_data

from alinea.astk.Weather import Weather
from alinea.astk.TimeControl import *
from alinea.caribu.caribu_star import rain_and_light_star
from alinea.alep.protocol import infect
from alinea.echap.interception_leaf import pesticide_applications
from alinea.echap.microclimate_leaf import microclimate_leaf
from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record #to avoid confusion with numpy record

from macros_annual_loop import setup_canopy, update_lesions, update_pesticides, dispersion, contamination, pesticide_intercept
import alinea.septo3d

from alinea.echap.tests_nodes import plot_pesticide


def pesticide_loop(meteo_file='meteo00-01.txt', start="2000-04-25", periods=8, freq='H', TB=0, delayH=1, delayDD=15, applications=""" """):
    # Fichier meteo
    meteo_path = shared_data(alinea.septo3d, meteo_file)
    # Setup
    weather = Weather(data_file=meteo_path)
    weather.check(['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
    # Pesticide protocol
    pest_calendar = pesticide_applications(applications)
    # Periode
    seq = pandas.date_range(start, periods=periods, freq=freq)  # 5000 et "2000-11-20" pour le cycle complet
    g, adel, domain, domain_area, convUnit, nplants = setup_canopy()
    recorder = LeafElementRecorder()
    # Delay
    TTmodel = DegreeDayModel(Tbase=TB)
    every_h = time_filter(seq, delay=delayH)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay=delayDD)
    every_pest = date_filter(seq, pest_calendar)
    # Time_Control
    canopy_timing = IterWithDelays(*time_control(seq, every_dd, weather.data))
    doses_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    pest_timing = IterWithDelays(*time_control(seq, every_pest, pest_calendar))
    # Simul
    for i,controls in enumerate(zip(canopy_timing, doses_timing, pest_timing)):
        canopy_iter, doses_iter, pest_iter = controls
        if canopy_iter:
            print canopy_iter.value.index[-1]
            g = adel.grow(g, canopy_iter.value)
            _=rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            print('recording, iter %d, tt%f'%(i,adel.canopy_age))
            _=do_record(g, canopy_iter.value, recorder, header={'iter':i, 'TT':adel.canopy_age})
        if pest_iter:
            print pest_iter.value
            _=pesticide_intercept(g, pest_iter.value)
        if doses_iter:
            print 'updte microclimate / doses...'
            _=microclimate_leaf(g, doses_iter.value, domain = domain, convUnit = convUnit)
            _=update_pesticides(g, doses_iter.value)
            #_=do_record(g, doses_iter.value, recorder, header={'iter':i, 'TT':adel.canopy_age})
    return g, recorder
      

def repartition_at_application(appdate = '2011-04-19', dose = 0.5, age = 1166):
    from macros_annual_loop import setup_canopy
    from alinea.echap.recorder import LeafElementRecorder
    
    recorder = LeafElementRecorder()
    g, adel, domain, domain_area, convUnit, nplants = setup_canopy(age=1166)
    applications= """date,dose, product_name
    %s 10:00:00, %f, Opus
    """%(appdate, dose)
    application_data = pesticide_applications(applications)
    g,_=pesticide_intercept(g, application_data)
    do_record(g, application_data, recorder)
    df =  recorder.get_records()
    return df