"""

Run annual loop by decomposing sub_runs and using pickle to restore

"""
import pandas
import pickle

from alinea.adel.newmtg import move_properties
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *

from alinea.echap.architectural_reconstructions import reconst_db

from alinea.echap.pesticide_model import update_pesticides, pesticide_intercept, pesticide_applications

from alinea.caribu.caribu_star import rain_and_light_star
from alinea.echap.microclimate_leaf import microclimate_leaf

from alinea.septo3d.disease import disease
from alinea.alep.protocol import *

from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record

weather = Boigneville_2010_2011()
Mercia = reconst_db['Mercia']
pgen, adel, domain, domain_area, convUnit, nplants = Mercia(nplants = 3, nsect=5)
tx = pgen['dynT_user'].a_cohort[0]
TTmodel = DegreeDayModel(Tbase = 0)

#modele
seq = pandas.date_range(start = "2010-11-02", periods=100, freq='H')  # hours = 5000 pour le cycle complet
# tout le cycle hivernal
#seq = pandas.date_range(start = "2010-11-02", end = "2011-07-01 18:00:00", freq='H')
applications = """date,dose, product_name
        2010-11-02 00:00:00, 0, bug
        2010-04-29 10:00:00, 1, Opus
        2010-04-29 11:00:00, 0, bug
        2010-05-22 10:00:00, 1, Opus
        2010-05-22 11:00:00, 0, bug
        """

pest_calendar = pesticide_applications(applications)
every_pest = date_filter(seq, pest_calendar)
every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 10)
every_rain = rain_filter(seq, weather)

every_all = filter_or([every_dd, every_rain, every_pest])

# save wheat and microclimate for the whole season
def make_canopy(dir = './Tem'):   
    canopy_timing = IterWithDelays(*time_control(seq, every_all, weather.data))
    g = adel.setup_canopy(age=150)
    rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)

    it = 0
    adel.save(g,it, dir)
    
    for control in canopy_timing:
        if control:
            print(control.value.index[-1])
            it += 1
            g = adel.grow(g, control.value)
            rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            microclimate_leaf(g, control.value, domain = domain, convUnit = convUnit)
            adel.save(g,it, dir)
            
        
def run_disease(SspoSol = 0.01, compute_star=False, dir = './Nopest'):

    fungus, growth_controler, infection_controler, emitter, transporter, inoc, contaminator = disease(domain, domain_area, SspoSol, compute_star=compute_star)
    canopy_timing = IterWithDelays(*time_control(seq, every_all, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    
    recorder = LeafElementRecorder()
    it = 0
    g, TT = adel.load(it, './Tem')
    do_record(g, weather.get_weather_start(seq), recorder, header={'iter':it, 'TT':TT, 'HS': TT * tx})
    adel.save(g, it, dir)

    
    for i, controls in enumerate(zip(canopy_timing, rain_timing)):
        canopy_iter, rain_iter = controls
        if canopy_iter:
            it += 1
            print(('canopy iter %d ...'%(it)))
            newg,TT = adel.load(it, './Tem')
            move_properties(g,newg)
            g = newg
            update(g, canopy_iter.dt, growth_controler, senescence_model=None, label='LeafElement', weather_data = canopy_iter.value)
            do_record(g, canopy_iter.value, recorder, header={'iter':it, 'TT':TT, 'HS': TT * tx})
            adel.save(g, it, dir)
            
        if rain_iter:
            wdata = rain_iter.value
            if wdata.rain.sum() > 0:
                print('raining...')
                g = external_contamination(g, inoc, contaminator, wdata)
                g = disperse(g, emitter, transporter, fungus.name, label='LeafElement', weather_data=wdata)
                infect(g, rain_iter.dt, infection_controler, label='LeafElement')
                       
    recorder.save_records(dir+'/records.csv')
    return g, recorder
  

def run_disease_and_decay(start = 0, refdir = './Nopest', dir = './doses'):

    fungus, growth_controler, infection_controler, emitter, transporter, inoc, contaminator = disease(domain, domain_area, 0, compute_star=False)
    
    canopy_timing = IterWithDelays(*time_control(seq, every_all, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    pest_timing = IterWithDelays(*time_control(seq, every_pest, pest_calendar))
    
    recorder = LeafElementRecorder()
    it = 0
    g, TT = adel.load(it, refdir)
    do_record(g, weather.get_weather_start(seq), recorder, header={'iter':it, 'TT':TT, 'HS': TT * tx})
    adel.save(g, it, dir)

    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, pest_timing)):
        canopy_iter, rain_iter, pest_iter = controls
        if canopy_iter:
            it += 1
            print(('canopy iter %d ...'%(it)))
            newg,TT = adel.load(it, refdir)
            move_properties(g,newg)
            g = newg
            update(g, canopy_iter.dt, growth_controler, senescence_model=None, label='LeafElement', weather_data = canopy_iter.value)
            update_pesticides(g, canopy_iter.value)
            do_record(g, canopy_iter.value, recorder, header={'iter':it, 'TT':TT, 'HS': TT * tx})
            adel.save(g, it, dir)
            
        if rain_iter:
            wdata = rain_iter.value
            if wdata.rain.sum() > 0:
                print('raining...')
                g = external_contamination(g, inoc, contaminator, wdata)
                g = disperse(g, emitter, transporter, fungus.name, label='LeafElement', weather_data=wdata)
                infect(g, rain_iter.dt, infection_controler, label='LeafElement')
        if pest_iter:
            if pest_iter.value.dose.sum() > 0:
                print('pesticide application...')
                pesticide_intercept(g, pest_iter.value)
                       
    recorder.save_records(dir+'/records.csv')
    return g, recorder
