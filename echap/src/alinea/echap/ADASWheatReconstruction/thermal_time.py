'''
Created on 25 oct. 2011

@author: cchambon
'''

import numpy
import pandas
from pandas.core.datetools import Hour
import math
import datetime

def parton_logan(dailyMinMaxTemp, latitude=55, param="air150cm"):
    '''Estimate hourly temperature from daily temperature. Based on Parton & Logan formula (AFM 1981) 
  
    :Parameters:
        - `dailyMinMaxTemp` : the min and max temperatures for each consecutive julian day number.
        - `latitude` : latitude in degrees
        - `param` : a dictionary of 3 coefficients 'a', 'b' and 'c' (e.g.: {'a': 1.86, 'b': 2.2, 'c': - 0.17}) 
        OR a string indicating the type of temperature recorded ("air150cm", "air10cm", "soilsurface", or "soil-10cm").
        a,b,c are for :
        a = time lag coeff for max Temperature after noon (h)
        b = coefficient for temperature decrease at night
        c = timeLag for the minimum temperature after sunrise (h)
        
    :Types:
        - `dailyMinMaxTemp` : pandas.DataFrame
        - `latitude` : float
        - `param` : dictionary
  
    :return: the air temperature at hourly time steps.
    :rtype: pandas.core.series.TimeSeries
    
    '''

    paramref = {'air150cm': {'a': 1.86, 'b': 2.2, 'c': - 0.17},
                'air10cm': {'a': 1.52, 'b': 2, 'c': - 0.18},
                'soilsurface': {'a': 0.5, 'b': 1.81, 'c': 0.49},
                'soil-10cm': {'a': 0.45, 'b': 2.28, 'c': 1.83}}
                 
    if type(param) is str:
        p = paramref[param]
    else:
        p = param
        
    a = p['a']
    b = p['b']
    c = p['c']

    phi = latitude / 180.0 * math.pi
        
    Tmin = dailyMinMaxTemp['Tmin']
    Tmax = dailyMinMaxTemp['Tmax']
    
    daily_series = pandas.Series(range(Tmin.size), index=Tmin.index)    
        
    def calc_daylength(i):
        dayNumber = dailyMinMaxTemp.index[i].toordinal() + 1 - datetime.datetime(dailyMinMaxTemp.index[i].year,1,1).toordinal()
        delta_i = 0.4014 * math.sin(2 * math.pi * (dayNumber - 77) / 365.0)
        t2_i = (-math.tan(phi)) * math.tan(delta_i)
        t1_i = math.sqrt(1 - t2_i ** 2)
        # angular amplitude
        ahr_i = math.atan2(t1_i, t2_i)
        return ahr_i / math.pi * 24
        
    daylength = daily_series.map(calc_daylength)
    
    def calc_sunrise(i):
        return 12 - daylength[i] / 2.0
        
    sunrise = daily_series.map(calc_sunrise)
    
    def calc_sunset(i):
        return 12 + daylength[i] / 2.0
        
    sunset = daily_series.map(calc_sunset)
    
    # minimal temperature hour
    def calc_hmin(i):
        return sunrise[i] + c
        
    hmin = daily_series.map(calc_hmin)
    
    # sunset temperature
    def calc_Tsunset(i):
        return Tmin[i] + (Tmax[i] - Tmin[i]) * math.sin(math.pi * (sunset[i] - hmin[i]) / (daylength[i] + 2 * a))
        
    Tsunset = daily_series.map(calc_Tsunset)    
    
    # sunset temperature at the day before
    Tsunsetb = pandas.Series(numpy.concatenate(([Tsunset.values[0]], Tsunset.values[:-1])), index=Tsunset.index)
    # minimal temperature at the day after
    Tmina = pandas.Series(numpy.concatenate((Tmin.values[1:], [Tmin.values[-1]])), index=Tsunset.index)
    
    hourly_idx = pandas.DateRange(Tmin.index[0], Tmin.index[-1] + 23 * Hour(), offset=Hour())

    hourly_series = pandas.Series(range(hourly_idx.size), index=hourly_idx)                               
                                          
    def calc_hourlyTemp(abs_hour):
        abs_day = abs_hour / 24
        rel_hour = abs_hour % 24 + 1
        if rel_hour < sunrise[abs_day]:
            return Tmin[abs_day] + (Tsunsetb[abs_day] - Tmin[abs_day]) * math.exp(-b * (24 - sunset[abs_day] + rel_hour) / (24 - daylength[abs_day]))
        elif rel_hour > sunset[abs_day]:
            return Tmina[abs_day] + (Tsunset[abs_day] - Tmina[abs_day]) * math.exp(-b * (rel_hour - sunset[abs_day]) / (24 - daylength[abs_day]))
        else:
            return Tmin[abs_day] + (Tmax[abs_day] - Tmin[abs_day]) * math.sin(math.pi * (rel_hour - hmin[abs_day]) / (daylength[abs_day] + 2 * a))
    
    temp_series = hourly_series.map(calc_hourlyTemp)
    return pandas.DataFrame(temp_series, index=temp_series.index, columns=['Tair'])


def dsT(Th, Tb=0.0):
    '''Hourly thermal time increment.
    
    :Parameters:
        - `Th` : ???
        - `Tb` : ???

    :Types:
        - `Th` : numpy.array
        - `Tb` : float
  
    :return: Hourly thermal time increment.
    :rtype: numpy.array
    
    '''
    return numpy.maximum(numpy.zeros_like(Th), Th - Tb)


def loiT(TK, k, EaR, DSR, DHR): 
    '''???
    
    :Parameters:
        - `TK` : temperature in Kelvin
        - `k` : ???
        - `EaR` : ???
        - `DSR` : ???
        - `DHR` : ???

    :Types:
        - `TK` : float
        - `k` : float
        - `EaR` : float
        - `DSR` : float
        - `DHR` : float
  
    :return: ???
    :rtype: float

    '''
    return k * TK * numpy.power(math.e, -EaR / TK) / (1 + numpy.power(math.e, DSR - DHR / TK))


def dsTc(Tair, TCref=12.0, k=3.8088e10, EaR=8899.6, DSR=68.0, DHR=20736.0):
    '''Hourly compensated thermal time increment.
    
    :Parameters:
        - `Tair` : temperature in Kelvin
        - `TCref` : ???
        - `k` : ???
        - `EaR` : ???
        - `DSR` : ???
        - `DSR` : ???

    :Types:
        - `Tair` : numpy.array
        - `TCref` : float
        - `k` : float
        - `EaR` : float
        - `DSR` : float
        - `DSR` : float
  
    :return: Hourly compensated thermal time increment.
    :rtype: numpy.array

    '''
    TKref = 273.0 + TCref
    TKref_mod = loiT(TKref, k, EaR, DSR, DHR)
    TK_array = Tair + 273.0
    TK_array_mod = loiT(TK_array, k, EaR, DSR, DHR)
    return TK_array_mod * TCref / TKref_mod / 24.0

        