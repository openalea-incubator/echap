# -*- coding: latin1 -*- 
import pandas
import numpy
import datetime
from openalea.deploy.shared_data import shared_data
import alinea.echap

from alinea.astk.Weather import Weather, linear_degree_days
from alinea.astk.TimeControl import DegreeDayModel
from openalea.deploy.shared_data import shared_data

def arvalis_reader(data_file):
    """ reads meteorological data from arvalis database 
    
    Caution : missing data labeled as 'Abs!' should be removed before parsing
    
    To Check  if time is UTC ?
    """
    
    data = pandas.read_csv(data_file, parse_dates={'datetime':[0,1]}, dayfirst=True, names=['date','h','temperature_air','relative_humidity','rain','wind_speed','global_radiation'], sep=';', index_col=0,skiprows=2, decimal=',')
    #convert Rg J/cm2 -> J.m-2.s-1
    data['global_radiation'] *= (10000. / 3600)
    #convert wind km/h -> m.s-1
    data['wind_speed'] *= (1000. / 3600)
    return data
  
#to do : one function with year as argument
  
def Boigneville_2010_2011():
    meteo_path = shared_data(alinea.echap, 'Boigneville_01092010_31082011_h.csv')
    weather = Weather(meteo_path, reader = arvalis_reader)
    weather.check(['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
    return weather
    
def Boigneville_2011_2012():
    meteo_path = shared_data(alinea.echap, 'Boigneville_01092011_31082012_h.csv')
    weather = Weather(meteo_path, reader = arvalis_reader)
    weather.check(['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
    return weather
    
def Boigneville_2012_2013():
    meteo_path = shared_data(alinea.echap, 'Boigneville_01092012_31082013_h.csv')
    weather = Weather(meteo_path, reader = arvalis_reader)
    weather.check(['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
    return weather
    
def format_date(start_date, end_date):
    """ Turn string start and end dates in datetime format to read weather """
    start_date+=" 12:00:00"
    end_date+=" 01:00:00"
    start = datetime.datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end = datetime.datetime.strptime(end_date, '%Y-%m-%d %H:%M:%S')
    return start, end

def read_weather(start, end):
    """ Read weather and add variables relevant for epidemiological analysis """
    if start.year >= 2010:
        filename = 'Boigneville_0109'+str(start.year)+'_3108'+str(end.year)+'_h.csv'
        meteo_path = shared_data(alinea.echap, filename)
        weather = Weather(meteo_path, reader = arvalis_reader)
        weather.check(['temperature_air', 'PPFD', 'relative_humidity',
                       'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
    weather.check(varnames=['degree_days'], models={'degree_days':linear_degree_days},
                    start_date = start)
    return weather

def read_weather_year(year = 2011):
    if year == 2011:
        start_date = "2010-10-15"
        end_date = "2011-06-20"
    elif year == 2012:
        start_date = "2011-10-21"
        end_date = "2012-07-18"
    elif year == 2013:
        start_date = "2012-10-29"
        end_date = "2013-08-01"
    start, end = format_date(start_date, end_date)
    return read_weather(start, end)

def mat_ray_obs(name = 'Tremie12'):
    import re
    if name is 'Mercia' or name is 'Rht3':
        data_file = 'METEO_stationINRA_20102011.csv'
    elif name is 'Tremie12': 
        data_file = 'METEO_stationINRA_20112012.csv'
    elif name is 'Tremie13':
        data_file = 'METEO_stationINRA_20122013.csv'
    
    header_row = ['DATE','H','RECORD','QR_PAR_Avg','SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4','SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8']
    #data = pandas.read_csv(data_file, parse_dates={'datetime':[0,1]}, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')
    data_path = shared_data(alinea.echap, data_file)
    df = pandas.read_csv(data_file, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')

    #filtre seulement heure entre 9h et 19h pour chaque journee
    df = df[df['H']>9]; df = df[df['H']<19]
    
    #tableau non traité
    dfa=df; dfa = dfa.reset_index()
    dfa.columns = ['datetime','H','RECORD','QR_PAR_Avg','SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4','SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8']
    dfa = dfa.groupby(dfa['datetime'], axis=0).mean(); dfa = dfa.reset_index()
    da=0
    while da<len(dfa):
        dfa['datetime'][da] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", dfa['datetime'][da])
        da = da+1
    dfa['%SE1'] = dfa['SE_PAR_Avg_1']/dfa['QR_PAR_Avg']; dfa['%SE2'] = dfa['SE_PAR_Avg_2']/dfa['QR_PAR_Avg']
    dfa['%SE3'] = dfa['SE_PAR_Avg_3']/dfa['QR_PAR_Avg']; dfa['%SE4'] = dfa['SE_PAR_Avg_4']/dfa['QR_PAR_Avg']
    dfa['%SE5'] = dfa['SE_PAR_Avg_5']/dfa['QR_PAR_Avg']; dfa['%SE6'] = dfa['SE_PAR_Avg_6']/dfa['QR_PAR_Avg']
    dfa['%SE7'] = dfa['SE_PAR_Avg_7']/dfa['QR_PAR_Avg']; dfa['%SE8'] = dfa['SE_PAR_Avg_8']/dfa['QR_PAR_Avg']
    
    # traitement tableau de données obs
    #avg PAR niveau du sol ('SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4')
    dat0 = df.ix[:,'SE_PAR_Avg_1':'SE_PAR_Avg_4']; dat0 = dat0.apply(numpy.mean,1)
    #avg PAR env 20 cm du sol ('SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8')
    dat20 = df.ix[:,'SE_PAR_Avg_5':'SE_PAR_Avg_8']; dat20 = dat20.apply(numpy.mean,1)
    # tableau % dat0/dat200
    tab_prc0 = dat0/df.QR_PAR_Avg
    # tableau % dat20/dat200
    tab_prc20 = dat20/df.QR_PAR_Avg
    
    # chaine de traitement pour hauteur au niveau du sol
    tab0 = pandas.DataFrame(tab_prc0)
    tab0 = tab0.reset_index(); tab0.columns = ['datetime','%']
    tab0 = tab0.groupby(tab0['datetime'], axis=0).mean()
    tab0 = tab0.reset_index()
    # changement format date JJ/MM/AAAA en JJ-MM-AAAA
    dt=0
    while dt<len(tab0):
        tab0['datetime'][dt] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", tab0['datetime'][dt])
        dt = dt+1

    # chaine de traitement pour hauteur a 20 cm du sol
    tab20 = pandas.DataFrame(tab_prc20)
    tab20 = tab20.reset_index(); tab20.columns = ['datetime','%']
    tab20 = tab20.groupby(tab20['datetime'], axis=0).mean()
    tab20 = tab20.reset_index()
    # changement format date JJ/MM/AAAA en JJ-MM-AAAA
    dte=0
    while dte<len(tab20):
        tab20['datetime'][dte] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", tab20['datetime'][dte])
        dte = dte+1
    
    return dfa, tab0, tab20