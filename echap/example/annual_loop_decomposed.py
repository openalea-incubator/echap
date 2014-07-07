"""

Run annual loop by decomposing sub_runs and using pickle to restore

"""
import pandas
import pickle

from alinea.adel.newmtg import move_properties

from alinea.echap.architectural_reconstructions import reconst_db
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *

from alinea.caribu.caribu_star import rain_and_light_star

from alinea.septo3d.disease import disease
from alinea.alep.protocol import *

from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record

weather = Boigneville_2010_2011()
Mercia = reconst_db['Mercia']
pgen, adel, domain, domain_area, convUnit, nplants = Mercia(nplants = 3, nsect=5)
tx = pgen['dynT_user'].a_cohort[0]
TTmodel = DegreeDayModel(Tbase = 0)



seq = pandas.date_range(start = "2010-11-02", periods=100, freq='H')  # hours = 5000 pour le cycle complet
every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 10)
every_rain = rain_filter(seq, weather)
every_dd_or_rain = filter_or([every_dd, every_rain])

def make_canopy():   
    canopy_timing = IterWithDelays(*time_control(seq, every_dd_or_rain, weather.data))
    g = adel.setup_canopy(age=150)
    rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)

    it = 0
    adel.save(g, it)
    
    for control in canopy_timing:
        if control:
            print control.value.index[-1]
            it += 1
            g = adel.grow(g, control.value)
            rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            adel.save(g,it)

def lesions_to_dict(g):
    lesions = g.property('lesions')
    for vid, les in lesions.iteritems():
        lesions[vid] = lesions.[vid].c_to_dict()
        
        
def dict_to_lesions(g):
    lesions = g.property('lesions')
    for vid, les in lesions.iteritems():
        lesions[vid] = lesions.[vid].dict_to_c()
        
def run_disease(SspoSol = 0.01, compute_star=False):

    fungus, growth_controler, infection_controler, emitter, transporter, inoc, contaminator = disease(domain, domain_area, SspoSol, compute_star=compute_star)
    canopy_timing = IterWithDelays(*time_control(seq, every_dd_or_rain, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    
    recorder = LeafElementRecorder()
    it = 0
    g,TT = adel.load(it)
    do_record(g, weather.get_weather_start(seq), recorder, header={'iter':it, 'TT':TT, 'HS': TT * tx})
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing)):
        canopy_iter, rain_iter = controls
        if canopy_iter:
            it += 1
            newg,TT = adel.load(it)
            lesions_to_dict(g)
            move_properties(g,newg)
            g = newg
            dict_to_lesions(g)
            update(g, canopy_iter.dt, growth_controler, senescence_model=None, label='LeafElement', weather_data = canopy_iter.value)
            do_record(g, canopy_iter.value, recorder, header={'iter':it, 'TT':TT, 'HS': TT * tx})  
        if rain_iter:
            print 'raining...'
            wdata = rain_iter.value
            g = external_contamination(g, inoc, contaminator, wdata)
            g = disperse(g, emitter, transporter, fungus.name, label='LeafElement', weather_data=wdata)
            infect(g, rain_iter.dt, infection_controler, label='LeafElement')
            
    return g,recorder
    
    