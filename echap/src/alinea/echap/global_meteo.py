# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 14:29:15 2013

@author: lepse
"""

import pandas as pd
from datetime import datetime, timedelta


def parse(yr, doy, hr):
    """ Convert the 'An', 'Jour' and 'hhmm' variables of the meteo dataframe in a datetime object (%Y-%m-%d %H:%M:%S format)
    """
    an, jour, heure = [int(x) for x in [yr, doy, hr/100]]
    dt = datetime(an - 1, 12, 31)
    delta = timedelta(days=jour, hours=heure)
    return dt + delta


class Meteo(object):
    """ Class compliying echap local_microclimate model protocol (meteo_reader).
    """
    def __init__(self, data_file = 'E:/openaleapkg/echap/test/meteo01.csv'):
        self.data = pd.read_csv(data_file, parse_dates={'datetime':['An','Jour','hhmm']},date_parser=parse,delimiter=';',usecols=['An','Jour','hhmm','PAR','Tair','HR','Vent','Pluie'])
    def get_meteo_file(self, timestep, t_deb):
        """ Read an hourly meteo file and return the global climate averaged and the global climate detail dataframe for an hour time step and the start date of each time step.

        :Parameters:
        ----------

        - `timestep` - The time step on which the averages are performed and for which the data frame is return
        - `t_deb` - The start date for reading the meteo file 

        :Returns:
        -------

        - `mean_globalclimate` - Mean variables of the global climate dataframe
        - `globalclimate` - Pandas dataframe with hourly meteo for the time step from t_deb
        - `t_deb` - The new start date for the new read of the meteo file 
        """    
        df_data = self.data
        index_deb = df_data.index[df_data['datetime']==t_deb][0]
        globalclimate = df_data.truncate(before = index_deb, after = index_deb + (timestep - 1))
        mean_globalclimate = globalclimate.set_index('datetime').mean()
        date_object = datetime.strptime(t_deb, '%Y-%m-%d %H:%M:%S')
        d = date_object + timedelta(hours=timestep)
        t_deb = datetime.strftime(d, '%Y-%m-%d %H:%M:%S')
        return mean_globalclimate, globalclimate, t_deb




