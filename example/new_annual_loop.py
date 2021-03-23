import pandas
import matplotlib.pyplot as plt

plt.ion()

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

from macros_annual_loop import reconst, update_lesions, update_pesticides, dispersion, contamination, pesticide_intercept
from alinea.echap.tests_nodes import plot_pesticide

import alinea.septo3d
from alinea.astk.Weather import Weather
from alinea.echap.weather_data import *

##------------------------------------------------------------------------------------------------------------

def plot(dataframefile, what, t, by, axe, plant, style):
    from itertools import cycle, islice
    d = dataframefile
    d = d[d[what].notnull()]
    if plant == 'all':
        dm=d[d.axe==axe]
    else:
        dm=d[(d.axe==axe) & (d.plant==plant)]
    gr=dm.groupby(t)
    dmp=gr.agg(numpy.sum)
    dmp = dmp.reset_index()
    '''gr=dmp.groupby([by,t])
    dms=gr.agg(numpy.mean)
    dms=dms.reset_index()
    gr=dms.groupby(by)'''
    dmp.plot(t, what, style=style)
    
def date(seq, weather, tx):
    bid = thermal_time(seq, weather.data)
    bid = bid[seq]
    bid = pandas.DataFrame(list(bid.values), index=bid.index)
    bid = bid.reset_index()
    bid.columns = ['date','TT']
    #bid['HS'] = bid['TT']*tx
    return bid

def annual_loop(name='Mercia', dTT_stop=0, original=False, nplants=3):

    if name is 'Mercia' or 'Rht3':
        weather = Boigneville_2010_2011()
        periods = 5000 # 5000 pour le cycle complet
        seq = pandas.date_range(start = "2010-11-02", periods=periods, freq='H')  
        seq_bid = pandas.date_range(start = "2010-10-15", periods=periods+7000, freq='H')
        applications = """date,dose, product_name
        2010-11-02 00:00:00, 0, bug
        2010-11-02 1:00:00, 1, Opus
        2010-11-02 2:00:00, 0, bug
        2010-04-29 10:00:00, 1, Opus
        2010-04-29 11:00:00, 0, bug
        2010-05-22 10:00:00, 1, Opus
        2010-05-22 11:00:00, 0, bug
        """

    pest_calendar = pesticide_applications(applications)

    pgen, adel, domain, domain_area, convUnit, nplants = reconst(name=name, dTT_stop=dTT_stop, original=original, n=nplants)
    tx = pgen['dynT_user'].a_cohort[0]
    SspoSol = 0.01
    recorder = LeafElementRecorder()

    TTmodel = DegreeDayModel(Tbase = 150)
    every_rain = rain_filter(seq, weather)
    every_h = time_filter(seq, delay = 6)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 15)
    every_pest = date_filter(seq, pest_calendar)

    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    canopy_timing = IterWithDelays(*time_control(seq, every_dd, weather.data))
    doses_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    pest_timing = IterWithDelays(*time_control(seq, every_pest, pest_calendar))

    g = adel.setup_canopy(age=150)

    for i,controls in enumerate(zip(canopy_timing, doses_timing, rain_timing, pest_timing)):
        canopy_iter, doses_iter, rain_iter, pest_iter = controls
        
        if canopy_iter:
            print('--', canopy_iter.value.index[-1], '--')
            g = adel.grow(g, canopy_iter.value)
            _=rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            _=update_lesions(g, canopy_iter.dt, True)
            _=do_record(g, canopy_iter.value, recorder)
        if pest_iter:
            _=pesticide_intercept(g, pest_iter.value)
            #plot_pesticide(g)
        if doses_iter:
            print('-- update microclimate / doses --')
            _=microclimate_leaf(g, doses_iter.value, domain = domain, convUnit = convUnit)
            _=update_pesticides(g, doses_iter.value)
            #plot_pesticide(g)
            _=do_record(g, doses_iter.value, recorder, header={'iter':i, 'TT':adel.canopy_age})
        if rain_iter:
            print('-- rain --')
            wdata = rain_iter.value
            dispersion(g, wdata, domain, domain_area, convUnit)
            _=contamination(g,wdata, SspoSol, domain, domain_area, convUnit)
            _=infect(g, rain_iter.dt)

    # recuperer recorder sous forme csv    
    #recorder.save_records('test.csv')   

    # remplissage correct de la colonne TT
    rec = recorder.get_records()
    rec = rec[rec.TT >= 0] 
    
    rec['date'] = rec['date'].map(lambda x: x.strftime('%Y-%m-%d'))
    rec = rec.drop('TT', 1)
    bid = date(seq_bid,weather,tx)
    bid['date'] = bid['date'].map(lambda x: x.strftime('%Y-%m-%d'))
    bid = bid.groupby('date').mean()
    bid.to_csv('bid.csv')
    bid=bid.reset_index()
    
    rec = pandas.merge(bid,rec,on='date')
    rec.to_csv('rec_bid.csv')
    
    return rec, bid
    
    #plt.plot(rec, what='green_area', t = 'date', by = 'ntop', axe = 'MS', plant = 'all', style = 'ob')