from alinea.echap.imports_echap import *


def test_intercept_caribu():
    g = adel_mtg()
    pest_calendar = {'datetime':["2000-10-01 01:00:00"], 'dose':[1.5], 'product_name':['Opus']}
    interception_model = CaribuInterceptModel(pest_calendar=pest_calendar)
    tc = TimeControl(steps = 1, model = interception_model, start_date = "2000-10-01 01:00:00")
    time_control = tc.next()
    g = pesticide_interception(g, interception_model, time_control)
    return g


def test_intercept_popdrops():    
    g = adel_mtg()
    pest_calendar = {'datetime':["2000-10-01 01:00:00"], 'dose':[1.5], 'product_name':['Opus']}
    interception_model = PesticideInterceptionModel(pest_calendar=pest_calendar)
    tc = TimeControl(steps = 1, model = interception_model, start_date = "2000-10-01 01:00:00")
    time_control = tc.next()
    g = pesticide_interception(g, interception_model, time_control)
    return g