#from alinea.echap.imports_echap import *
from alinea.echap.wheat_mtg import *
from openalea.deploy.shared_data import shared_data
from alinea.echap.interception_leaf import CaribuInterceptModel
from alinea.echap.TimeControl import TimeControl
from alinea.astk.Weather import Weather
from alinea.popdrops.PesticideInterception import PesticideInterceptionModel
from alinea.echap.interfaces import pesticide_interception, rain_interception
import alinea.echap

############# Pesticide interception
def test_intercept_caribu():
    g = adel_mtg()
    pest_calendar = {'datetime':["2000-10-01 01:00:00"], 'dose':[1.5], 'product_name':['Opus']}
    interception_model = CaribuInterceptModel(pest_calendar=pest_calendar)
    tc = TimeControl(steps = 1, model = interception_model, start_date = "2000-10-01 01:00:00")
    time_control = next(tc)
    g = pesticide_interception(g, interception_model, time_control)
    return g


def test_intercept_popdrops():    
    g = adel_mtg()
    pest_calendar = {'datetime':["2000-10-01 01:00:00"], 'dose':[1.5], 'product_name':['Opus']}
    interception_model = PesticideInterceptionModel(pest_calendar=pest_calendar)
    tc = TimeControl(steps = 1, model = interception_model, start_date = "2000-10-01 01:00:00")
    time_control = next(tc)
    g = pesticide_interception(g, interception_model, time_control)
    return g


############# Rain interception
def test_intercept_rain():    
    g = adel_mtg()
    meteo01_filepath = shared_data(['alinea/echap'])/ 'meteo01.csv'
    t_deb = "2000-10-01 01:00:00"
    weather = Weather(data_file=meteo01_filepath)
    interception_model = RainInterceptionModel()
    rain_timing = TimeControl(delay = 1, steps = 1, model = interception_model, weather = weather, start_date = t_deb)
    time_control = next(rain_timing)
    g = rain_interception(g, interception_model, time_control)
    return g