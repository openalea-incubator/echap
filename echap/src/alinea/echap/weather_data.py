# -*- coding: latin1 -*-
""" Readers to weathers recors"""
import pandas
import numpy
import os
import scipy.stats

from alinea.astk.Weather import Weather
from openalea.deploy.shared_data import shared_data
import alinea.echap


def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor can be looked up from a table, for 95pcent confidence
    on a reasonable size sample (>=500) 1.96 is used.
    """

    n, v = len(lst), numpy.var(lst, ddof=1)
    c = scipy.stats.t.interval(perc_conf * 1.0 / 100, n - 1)[1]

    return numpy.sqrt(v / n) * c


boigneville = {'city': 'Boigneville', 'latitude': 48.3, 'longitude': 2.3,
               'altitude': 100}


def arvalis_reader(data_file):
    """ reads meteorological data from arvalis database 
    
    Caution : missing data labeled as 'Abs!' should be removed before parsing
    
    To Check  if time is UTC ?
    """

    data = pandas.read_csv(data_file, names=['date', 'h', 'temperature_air',
                                             'relative_humidity', 'rain',
                                             'wind_speed', 'global_radiation'],
                           sep=';', skiprows=2, decimal=',')
    # create datetime index
    data.index = pandas.to_datetime(data['date'].map(str) + ' ' + data['h'],
                                    dayfirst=True)
    data['date'] = data.index
    # convert Rg J/cm2 -> J.m-2.s-1
    data['global_radiation'] *= (10000. / 3600)
    # convert wind km/h -> m.s-1
    data['wind_speed'] *= (1000. / 3600)
    return data


def get_weather(variety='Mercia'):
    if variety in ['Mercia', 'Rht3']:
        year = 2010
    elif variety == 'Tremie12':
        year = 2011
    elif variety == 'Tremie13':
        year = 2012
    else:
        raise ValueError('No year associated to this variety: ' + variety)

    filename = 'Boigneville_' + str(year) + '_' + str(year + 1) + '_h.csv'
    meteo_path = shared_data(alinea.echap) / 'weather_data' / filename
    weather = Weather(meteo_path, reader=arvalis_reader,
                      localisation=boigneville)
    weather.check(
        ['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain',
         'global_radiation', 'vapor_pressure'])

    return weather


def par_sensors_hourly(variety='Tremie12'):
    if variety is 'Mercia' or variety is 'Rht3':
        data_file = 'METEO_stationINRA_20102011.csv'
    elif variety is 'Tremie12':
        data_file = 'METEO_stationINRA_20112012.csv'
    elif variety is 'Tremie13':
        data_file = 'METEO_stationINRA_20122013.csv'

    header_row = ['date', 'h', 'record', 'QR_PAR_Avg', 'SE_PAR_Avg_1',
                  'SE_PAR_Avg_2', 'SE_PAR_Avg_3', 'SE_PAR_Avg_4',
                  'SE_PAR_Avg_5', 'SE_PAR_Avg_6', 'SE_PAR_Avg_7',
                  'SE_PAR_Avg_8']
    data_path = shared_data(alinea.echap) / 'weather_data' / data_file
    df = pandas.read_csv(data_path, names=header_row, sep=';',
                         skiprows=1, decimal='.')
    df['daydate'] = pandas.to_datetime(df['date'].values,
                                       dayfirst=True).strftime('%Y-%m-%d')
    df.pop('date')
    return df.groupby(('daydate', 'h')).aggregate(numpy.mean).reset_index()


def par_sensors_daily(variety='Tremie12'):
    df = par_sensors_hourly(variety)
    df.pop('record'); df.pop('h')
    return df.groupby('daydate').agg(sum).reset_index()


def par_transmittance(variety='Tremie12'):
    df = par_sensors_daily(variety)
    df0 = df.ix[:, 'SE_PAR_Avg_1':'SE_PAR_Avg_4'].apply(lambda x: x / df['QR_PAR_Avg'])
    df20 = df.ix[:, 'SE_PAR_Avg_5':'SE_PAR_Avg_8'].apply(lambda x: x / df['QR_PAR_Avg'])
    aggregated = df.loc[:,('daydate', 'QR_PAR_Avg')]
    aggregated['po_0_mean'] = df0.mean(axis=1)
    aggregated['po_0_conf_int'] = df0.apply(conf_int, axis=1)
    aggregated['po_20_mean'] = df20.mean(axis=1)
    aggregated['po_20_conf_int'] = df20.apply(conf_int, axis=1)
    return aggregated


# def mat_ray_obs(name='Tremie12'):
#     import re
#     if name is 'Mercia' or name is 'Rht3':
#         data_file = 'METEO_stationINRA_20102011.csv'
#     elif name is 'Tremie12':
#         data_file = 'METEO_stationINRA_20112012.csv'
#     elif name is 'Tremie13':
#         data_file = 'METEO_stationINRA_20122013.csv'
#
#     header_row = ['DATE', 'H', 'RECORD', 'QR_PAR_Avg', 'SE_PAR_Avg_1',
#                   'SE_PAR_Avg_2', 'SE_PAR_Avg_3', 'SE_PAR_Avg_4',
#                   'SE_PAR_Avg_5', 'SE_PAR_Avg_6', 'SE_PAR_Avg_7',
#                   'SE_PAR_Avg_8']
#     # data = pandas.read_csv(data_file, parse_dates={'datetime':[0,1]}, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')
#     data_path = shared_data(alinea.echap, data_file)
#     df = pandas.read_csv(data_path, dayfirst=True, names=header_row, sep=';',
#                          index_col=0, skiprows=1, decimal='.')
#
#     # filtre seulement heure entre 9h et 19h pour chaque journee
#     df = df[df['H'] > 9];
#     df = df[df['H'] < 19]
#
#     # tableau non traité
#     dfa = df;
#     dfa = dfa.reset_index()
#     dfa.columns = ['datetime', 'H', 'RECORD', 'QR_PAR_Avg', 'SE_PAR_Avg_1',
#                    'SE_PAR_Avg_2', 'SE_PAR_Avg_3', 'SE_PAR_Avg_4',
#                    'SE_PAR_Avg_5', 'SE_PAR_Avg_6', 'SE_PAR_Avg_7',
#                    'SE_PAR_Avg_8']
#     dfa = dfa.groupby(dfa['datetime'], axis=0).mean();
#     dfa = dfa.reset_index()
#     da = 0
#     while da < len(dfa):
#         dfa['datetime'][da] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})",
#                                      "\\1-\\2-\\3", dfa['datetime'][da])
#         da = da + 1
#     dfa['%SE1'] = dfa['SE_PAR_Avg_1'] / dfa['QR_PAR_Avg'];
#     dfa['%SE2'] = dfa['SE_PAR_Avg_2'] / dfa['QR_PAR_Avg']
#     dfa['%SE3'] = dfa['SE_PAR_Avg_3'] / dfa['QR_PAR_Avg'];
#     dfa['%SE4'] = dfa['SE_PAR_Avg_4'] / dfa['QR_PAR_Avg']
#     dfa['%SE5'] = dfa['SE_PAR_Avg_5'] / dfa['QR_PAR_Avg'];
#     dfa['%SE6'] = dfa['SE_PAR_Avg_6'] / dfa['QR_PAR_Avg']
#     dfa['%SE7'] = dfa['SE_PAR_Avg_7'] / dfa['QR_PAR_Avg'];
#     dfa['%SE8'] = dfa['SE_PAR_Avg_8'] / dfa['QR_PAR_Avg']
#
#     # traitement tableau de données obs
#     # avg PAR niveau du sol ('SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4')
#     dat0 = df.ix[:, 'SE_PAR_Avg_1':'SE_PAR_Avg_4'];
#     dat0 = dat0.apply(numpy.mean, 1)
#     # avg PAR env 20 cm du sol ('SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8')
#     dat20 = df.ix[:, 'SE_PAR_Avg_5':'SE_PAR_Avg_8'];
#     dat20 = dat20.apply(numpy.mean, 1)
#     # tableau % dat0/dat200
#     tab_prc0 = dat0 / df.QR_PAR_Avg
#     # tableau % dat20/dat200
#     tab_prc20 = dat20 / df.QR_PAR_Avg
#
#     # chaine de traitement pour hauteur au niveau du sol
#     tab0 = pandas.DataFrame(tab_prc0)
#     tab0 = tab0.reset_index();
#     tab0.columns = ['datetime', '%']
#     tab0 = tab0.groupby(tab0['datetime'], axis=0).mean()
#     tab0 = tab0.reset_index()
#     # changement format date JJ/MM/AAAA en JJ-MM-AAAA
#     dt = 0
#     while dt < len(tab0):
#         tab0['datetime'][dt] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})",
#                                       "\\1-\\2-\\3", tab0['datetime'][dt])
#         dt = dt + 1
#
#     # chaine de traitement pour hauteur a 20 cm du sol
#     tab20 = pandas.DataFrame(tab_prc20)
#     tab20 = tab20.reset_index();
#     tab20.columns = ['datetime', '%']
#     tab20 = tab20.groupby(tab20['datetime'], axis=0).mean()
#     tab20 = tab20.reset_index()
#     # changement format date JJ/MM/AAAA en JJ-MM-AAAA
#     dte = 0
#     while dte < len(tab20):
#         tab20['datetime'][dte] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})",
#                                         "\\1-\\2-\\3", tab20['datetime'][dte])
#         dte = dte + 1
#
#     return dfa, tab0, tab20
