import pandas
from openalea.deploy.shared_data import shared_data
import alinea.echap

from alinea.astk.Weather import Weather

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