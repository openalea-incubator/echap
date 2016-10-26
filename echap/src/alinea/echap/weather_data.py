# -*- coding: latin1 -*- 
import pandas
import numpy
try:
    import cPickle as pickle
except ImportError:
    import pickle
from alinea.astk.Weather import Weather, linear_degree_days
from openalea.deploy.shared_data import shared_data
import alinea.echap


boigneville={'city': 'Boigneville', 'latitude': 48.3,
                               'longitude': 2.3,
             'altitude':100}

def arvalis_reader(data_file):
    """ reads meteorological data from arvalis database 
    
    Caution : missing data labeled as 'Abs!' should be removed before parsing
    
    To Check  if time is UTC ?
    """
    
    data = pandas.read_csv(data_file, names=['date','h','temperature_air','relative_humidity','rain','wind_speed','global_radiation'], sep=';',skiprows=2, decimal=',')
    # create datetime index
    data.index = pandas.to_datetime(data['date'].map(str) + ' ' + data['h'], dayfirst=True)
    data['date'] = data.index
    #convert Rg J/cm2 -> J.m-2.s-1
    data['global_radiation'] *= (10000. / 3600)
    #convert wind km/h -> m.s-1
    data['wind_speed'] *= (1000. / 3600)
    return data


def get_weather(variety='Mercia'):
    if variety in ['Mercia', 'Rht3']:
        year = 2010
    elif variety is 'Tremie12':
        year =  2011
    elif variety is 'Tremie13':
        year = 2012

    filename = 'Boigneville_' + str(year) +  '_' + str(year + 1) + '_h.csv'
    meteo_path = shared_data(alinea.echap) / 'weather_data' / filename
    weather = Weather(meteo_path, reader=arvalis_reader, localisation=boigneville)
    weather.check(
        ['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain',
         'global_radiation', 'vapor_pressure'])

    return weather


def tt_lin(variety='Mercia'):
    """Build or read  in the cache the daydate/TT time scales"""

    sowing = {'Mercia': '2010-10-05', 'Rht3': '2010-10-05', 'Tremie12': '2011-10-21',
              'Tremie13': '2012-10-29'}
    harvest = {'Mercia': '2011-06-20', 'Rht3': '2011-06-20', 'Tremie12': '2012-06-19',
              'Tremie13': '2013-06-09'}
    filename = {'Mercia': 'MerciaRht3_TTlin_sowing_harvest.csv',
               'Rht3': 'MerciaRht3_TTlin_sowing_harvest.csv',
               'Tremie12': 'Tremie12_TTlin_sowing_harvest.csv',
               'Tremie13': 'Tremie13_TTlin_sowing_harvest.csv'}

    df = None
    path = shared_data(alinea.echap) / 'cache' / filename[variety]
    try:
        df = pandas.read_csv(path)
    except IOError:
        w = get_weather(variety)
        tt = linear_degree_days(w.data, start_date=sowing[variety])
        seq = pandas.date_range(sowing[variety], harvest[variety], tz='UTC')
        df = pandas.DataFrame({'daydate': seq.strftime('%Y-%m-%d'),
                               'TT': tt.reindex(seq).values})
        df.to_csv(path, index=False)

    return df


def application_tag(variety, daydate, which='T1'):
    """ generate tag relative to pesticide application date

    :param variety: variety name
    :param daydate: a list of daydate strings
    :return: tags
    """

    sowing = {'Mercia': '2010-10-05', 'Rht3': '2010-10-05', 'Tremie12': '2011-10-21',
              'Tremie13': '2012-10-29'}
    t1 = {'Mercia': '2011-04-19','Rht3':'2011-04-19', 'Tremie12': '2012-04-11',
          'Tremie13': '2013-04-25'}
    t2 = {'Mercia': '2011-05-11', 'Rht3': '2011-05-11',
          'Tremie12': '2012-05-09', 'Tremie13': '2013-05-17'}


    if which == 'T1':
        origin = t1[variety]
    elif which == 'T2':
        origin = t2[variety]
    else:
        origin = sowing[variety]

    delta = (pandas.to_datetime(daydate) - pandas.to_datetime(origin)).days
    tags = numpy.where(delta == 0, which,
                       numpy.where(delta > 0, map(lambda x: which + '+' + str(x), delta),
                                   map(lambda x: which + '-' + str(abs(x)), delta)))
    return tags


def tt_hs_tag(variety='Mercia'):
    """A multiple time index"""
    filename = str(shared_data(alinea.echap)/'cache'/'HS_fit.pckl')
    hsfit = None
    try:
        with open(filename) as input:
             hsfit = pickle.load(input)
    except IOError:
        print('HaunStage needs to be fitted before calling this function')

    df = None
    path = shared_data(alinea.echap) / 'cache' / variety + '_TT_HS_tag.csv'
    try:
        df = pandas.read_csv(path)
    except IOError:
        df = tt_lin(variety)
        df['HS'] = hsfit[variety](df['TT'])
        for t in ('T1', 'T2'):
            df['tag' + t] = application_tag(variety, df['daydate'].values, t)
        df.to_csv(path, index=False)

    return df






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
    df = pandas.read_csv(data_path, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')

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