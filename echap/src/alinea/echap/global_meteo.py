# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 14:29:15 2013

@author: lepse
"""

import pandas as pd
from datetime import datetime, timedelta


def parse(yr, doy, hr):
    an, jour, heure = [int(x) for x in [yr, doy, hr/100]]
    dt = datetime(an - 1, 12, 31)
    delta = timedelta(days=jour, hours=heure)
    return dt + delta


class Meteo(object):
    """ 
    """
    def __init__(self, data_file = 'E:/openaleapkg/alep/test/meteo01.csv'):
        self.data = pd.read_csv(data_file, parse_dates={'datetime':['An','Jour','hhmm']},date_parser=parse,delimiter=';',usecols=['An','Jour','hhmm','PAR','Tair','HR','Vent','Pluie'])
    def get_meteo_file(self, timestep, t_deb):
        df_data = self.data
        index_deb = df_data.index[df_data['datetime']==t_deb][0]
        globalclimate = df_data.truncate(before = index_deb, after = index_deb + (timestep - 1))
        mean_globalclimate = globalclimate.set_index('datetime').mean()
        return mean_globalclimate, globalclimate, t_deb




