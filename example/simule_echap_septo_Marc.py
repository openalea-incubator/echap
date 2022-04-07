""" Simulation echap peticide eample with mercia"""

# import echap architectural reconstruction and setup canopy
from alinea.echap.architectural_reconstructions import EchapReconstructions

echap = EchapReconstructions()

mercia_arch = echap.get_reconstruction(name="Mercia",nplants=30)
g_mercia = mercia_arch.setup_canopy(age=1000)

# pesticide interception
from alinea.echap.interception_leaf import pesticide_applications , InterceptModel
from alinea.echap.interfaces import pesticide_interception, pesticide_efficacy, pesticide_penetrated_decay
from alinea.echap.pesticide_efficacy import PesticideEfficacyModel
from alinea.echap.milne_leaf import PenetratedDecayModel
from alinea.echap.tests_nodes import plot_pesticide

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
productsDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}
interception = InterceptModel(productsDB)

pesticide_interception(g=g_mercia, interception_model=interception,application_data=pest_calendar , label='LeafElement')
pgl.Viewer.display(plot_pesticide(g=g_mercia))

# pesticide penetration decay

milne_parameters = {'Epoxiconazole': {'Ae': 0.5 ,
                                      'Ke': 7.01,
                                      'decay_rate': 0.069,#d-1
                                      'dose_max_ha': 125,
                                      'Ap': 0.71,
                                      'Kp': 6.,
                                      'type_code': 2},
                    'Chlorothalonil':{'Ae': 1., 
                                      'Ke': 6.49, 
                                      'decay_rate': 0.01, 
                                      'dose_max_ha': 1000, 
                                      'Ap': 0, 
                                      'Kp': 0, 
                                      'type_code': 3}
                                
                    }

decayModel = PenetratedDecayModel(milne_parameters)
pesticide_penetrated_decay(g=g_mercia,decay_model=decayModel,weather_data=[10,20,30])

# peciticide efficacy

pest_efficacy_model= PesticideEfficacyModel()
pesticide_efficacy(g_mercia,efficacy_model=pest_efficacy_model,weather_data=[10,20,30],label='LeafElement')



# loop echap with weather data
from alinea.echap.architectural_reconstructions import EchapReconstructions
from alinea.echap.weather_data import get_weather

from alinea.alep.mini_models import linear_degree_days
from alinea.alep.protocol import update
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.disease_outputs import AdelWheatRecorder

from alinea.astk.TimeControl import DegreeDayModel, thermal_time_filter
from alinea.astk.TimeControl import IterWithDelays , time_control




import pandas

# setup weatherdata
weather = get_weather(variety="Mercia") # get weather data from echap for Mercia variety

## check weather data and add degree day variable using linear degree day model
weather.check(["degree_days",'temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'],models={'degree_days':linear_degree_days})

## weather index selection
weather_time_index = pandas.date_range(start ="2010-11-02", end='2011-05-22', freq='H')  # Range date index. among weather data select weather period for simulation

## canopy growth timing
TTmodel = DegreeDayModel(Tbase = 0) # Thermal times model
every_dd = thermal_time_filter(time_sequence=weather_time_index, weather=weather, model=TTmodel, delay = 15) # allows to filter thermal time data here each 15 degree day
canopy_timing = IterWithDelays(*time_control(time_sequence=weather_time_index, eval_filter=every_dd, data=weather.data)) # create data from weatherdata according to filter


#Setup pesticide application data

from alinea.echap.interception_leaf import pesticide_applications , InterceptModel
from alinea.echap.interfaces import pesticide_interception,pesticide_efficacy, pesticide_penetrated_decay
from alinea.astk.TimeControl import date_filter

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
productsDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}
interception = InterceptModel(productsDB)

## pesticide interception timing
every_pest = date_filter(weather_time_index, pest_calendar)
pest_timing = IterWithDelays(*time_control(weather_time_index, every_pest, pest_calendar))

# pesticide penetrated setup
from alinea.echap.milne_leaf import PenetratedDecayModel
decayModel = PenetratedDecayModel()

# pesticide efficacy setup
from alinea.echap.pesticide_efficacy import PesticideEfficacyModel
pest_efficacy_model = PesticideEfficacyModel()

# growth control
growth_controler = NoPriorityGrowthControl()

# recorder
from alinea.echap.recorder import LeafElementRecorder
recorder = LeafElementRecorder()
#recorder.record(g,date=0)

from alinea.echap.interfaces import record
# simulation
from alinea.echap.tests_nodes import plot_pesticide
import openalea.plantgl.all as pgl

echap = EchapReconstructions()
wheat = echap.get_reconstruction(name="Mercia", nplants=2)
g = wheat.setup_canopy(age=500)

for i,controls in enumerate(zip(canopy_timing, pest_timing)):

    canopy_iter, pest_iter = controls

    if canopy_iter.eval:
        print('-- update microclimate / canopy --')
        print('--', canopy_iter.value.index[-1], '--')
        g = wheat.grow(g, canopy_iter.value)
        _=update(g, canopy_iter.dt, growth_controler)
        _= wheat.plot(g)
        #recorder.record(g=g, date=i)
    
    if pest_iter.eval:
        print('-- pesticide interception --')
        print('--', pest_iter.value.index[-1], '--')
        _=pesticide_interception(g, interception, pest_iter.value, label='LeafElement')
        _=pesticide_penetrated_decay(g=g,decay_model=decayModel,weather_data=pest_iter.value)
        _=pesticide_efficacy(g=g,efficacy_model=pest_efficacy_model,weather_data=pest_iter.value,label='LeafElement')
        #_= pgl.Viewer.display(plot_pesticide(g))
    
    record(g, pest_iter.value, recorder, header={'iter':i, 'dd':wheat.canopy_age})
    
recorder.save_records(path='C:/Users/mlabadie/Desktop/recoder_test.csv')    

        
