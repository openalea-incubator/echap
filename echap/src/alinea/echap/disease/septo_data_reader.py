# -*- coding: latin1 -*- 

""" Read data files from disease measurements in Grignon on:
    - Mercia & Rht3 wheat in 2010/11
    - Tremie12 wheat in 2011/12
    - Tremie13 wheat in 2012/13 
    
    Files are stored in 'echap>src>data>disease_measurements>...'
    Files with disease measurements are all in the same format as 'tremie_12_ctrl'
   
    Suffix in filenames stand for:
    - 'ctrl' control treatment with no fungicide (tagged plants followed over the campaign)
    - 'ref' reference treatment with fungicide (tagged plants followed over the campaign)
    - 'global' global notations on control treatment with (occasional and destructive
    with bigger samples)
    - 'nec' total necrosis (apical+disease on tagged plants of reference treatment)
    - 'ctrl_ligulation' 
"""

# Generic imports
import pandas as pd
import numpy as np
import collections
import re

# Imports for weather
import datetime
import alinea.echap
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.echap.weather_data import arvalis_reader
from alinea.astk.Weather import Weather
from alinea.alep.alep_weather import (wetness_rapilly, add_septoria_infection_risk, 
                                    add_septoria_risk_with_event, linear_degree_days)
from alinea.alep.alep_time_control import add_notation_dates

# Imports for wheat reconstruction
from alinea.echap.architectural_reconstructions import EchapReconstructions

# Reading of weather files #########################################################################
def format_date(start_date, end_date):
    """ Turn string start and end dates in datetime format to read weather """
    start_date+=" 12:00:00"
    end_date+=" 01:00:00"
    start = datetime.datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end = datetime.datetime.strptime(end_date, '%Y-%m-%d %H:%M:%S')
    return start, end

def read_weather(start, end, wetness_duration_min = 10., temp_min = 0., temp_max = 25.):
    """ Read weather and add variables relevant for epidemiological analysis """
    if start.year >= 2010:
        filename = 'Boigneville_0109'+str(start.year)+'_3108'+str(end.year)+'_h.csv'
        meteo_path = shared_data(alinea.echap, filename)
        weather = Weather(meteo_path, reader = arvalis_reader)
        weather.check(['temperature_air', 'PPFD', 'relative_humidity',
                       'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
        notation_dates_file = shared_data(alinea.alep, 'notation_dates/notation_dates_'+str(end.year)+'.csv')
        weather.check(varnames=['notation_dates'], models={'notation_dates':add_notation_dates}, notation_dates_file = notation_dates_file)
    else:
        start_yr = str(start.year)[2:4]
        end_yr = str(end.year)[2:4]
        filename = 'meteo'+ start_yr + '-' + end_yr + '.txt'
        meteo_path = shared_data(alinea.septo3d, filename)
        weather = Weather(data_file=meteo_path)
    weather.check(varnames=['degree_days'], models={'degree_days':linear_degree_days}, start_date=start, base_temp=0., max_temp=30.)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    weather.check(varnames=['septo_infection_risk'], 
                  models={'septo_infection_risk':add_septoria_infection_risk}, 
                  wetness_duration_min = wetness_duration_min, temp_min = temp_min, temp_max = temp_max)
    weather.check(varnames=['septo_infection_risk_with_event'], 
                  models={'septo_infection_risk_with_event':add_septoria_risk_with_event}, 
                  wetness_duration_min = wetness_duration_min, temp_min = temp_min, temp_max = temp_max)
    return weather
    
# Reading of disease files and manipulation of generated pandas DataFrame ##########################
def str2int(x):
    """ Turn string format into integer when possible """
    try:
        return int(x)
    except:
        return x

def is_iterable(obj):
    return isinstance(obj, collections.Iterable)
    
def change_index(df, new_index = 'date'):
    """ Set new index on dataframe """
    df.reset_index(inplace=True)
    df.set_index(new_index, inplace = True)

def get_missing(numbers, start, end):
    """ Find missing numbers in list compared to range(start,end) """
    all_numbers = set(range(start, end))
    missing = all_numbers - set(numbers)
    return np.array(list(missing))
    
def get_fnl_by_plant(data):
    """ Group plants from data by final leaf number (fnl) 
    
        :Returns:
        - dict([(fnl, [plant_ids])])
    """
    # Get data for last date and top leaves
    df = data.xs(max(set(data.index)))
    df = df[df['num_leaf_top']==1]
    # Group plant ids by final number of leaves (FNL)
    grps = {fnl:[plant for plant in df[df['num_leaf_bottom']==fnl]['num_plant']] for fnl in set(df['num_leaf_bottom'])}
    # Find fnl if missing value for a leaf on last notation
    if len(sum(grps.itervalues(),[]))<max(data['num_plant']):
        miss = get_missing(np.sort(sum(grps.itervalues(),[])), start=1, end=max(data['num_plant']))
        for m in miss:
            if len(data[data['num_plant']==m])>0:
                fnl = max(data['num_leaf_bottom'][data['num_plant']==m]) + min(data['num_leaf_top'][data['num_plant']==m]) - 1
                grps[fnl].append(m)
    return grps               

def group_by_fnl(data, from_top = True):
    """ Split dataframe by final leaf number (fnl) of plants
    
        :Returns:
        - dict([(fnl, [dataframe_fnl])])
    """
    grps = get_fnl_by_plant(data)
    # Get data by FNL
    if from_top == True:
        return {fnl:data[data['num_plant'].isin(grps[fnl])] for fnl in grps.iterkeys()}
    else:
        groups = {}
        for fnl in grps.iterkeys():
            group = data[data['num_plant'].isin(grps[fnl])]
            groups[fnl] = group
        return groups
        
def get_fnls(adel):
    """ Get final leaf numbers over all plants in wheat canopy simulated with adel """
    return set(adel.phenT().set_index('num_plant').groupby(level=0).count()['n'] - 1)
    
def data_reader_2011(data_file='mercia_2011.csv'):       
    """ Read data of septoria measurements from Grignon 2011 (Mercia or Rht3 wheat) """
    file_path = shared_data(alinea.echap, 'disease_measurements/'+data_file)
    df = pd.read_csv(file_path, parse_dates={'date':[1]}, sep=';', index_col=[0])
    df['num_plant'] += 30*(df['bloc']-1)
    df['septo_necro'] = np.nan*df.septo_green
    df['num_leaf_bottom'] = 13 - df['num_leaf_top']
    return df
    
def data_reader_tremis(data_file='tremis_2012_tnt.csv', green_on_green = True):
    """ Read data of septoria measurements from Grignon on Tremie wheat (2012 or 2013) """
    
    def severity(septo_green_tot, necro, ratio_septo_necro):
        return septo_green_tot + necro*ratio_septo_necro/100
    
    def septo_necro(necro, ratio_septo_necro):
        return necro*ratio_septo_necro/100
    
    def septo_green_on_green(septo_green, necro):
        septo_green_on_green = map(lambda (x,y): 100*x / (100 - y) if y<100 else 0., zip(septo_green, necro))
        #septo_green_on_green[septo_green_on_green == np.inf] = 0.
        return septo_green_on_green
    
    file_path = shared_data(alinea.echap, 'disease_measurements/'+data_file)
    df = pd.read_csv(file_path, parse_dates={'date':[1]}, sep=';', index_col=[0])
    
    df['septo_green_tot'] = df['septo_green'].copy()
    if green_on_green == True:
        df['septo_green'] = septo_green_on_green(df.septo_green_tot, df.necro)
        
    if isinstance(df['ratio_septo_necro'][0], str):
        #df = df[(df['ratio_septo_necro']!='<') & (df['ratio_septo_necro']!='>') & (df['ratio_septo_necro']!='=')
         #       & (df['ratio_septo_necro']!='?') & (df['ratio_septo_necro']!='_')]
        indx = (df['ratio_septo_necro']!='<') & (df['ratio_septo_necro']!='>') & (df['ratio_septo_necro']!='=') & (df['ratio_septo_necro']!='?') & (df['ratio_septo_necro']!='_')
        if green_on_green == True:
            df['ratio_septo_necro'][~indx] = df['septo_green'][~indx]
        else:
            df['ratio_septo_necro'][~indx] = np.nan
        df.ratio_septo_necro = df.ratio_septo_necro.map(str2int)
    df['severity'] = severity(df.septo_green_tot, df.necro, df.ratio_septo_necro)
    df['septo_necro'] = septo_necro(df.necro, df.ratio_septo_necro)
    
    # Correct leaf number from bottom in tremie 2012
    # if int(re.search(r'\d+', data_file).group()) == 2012 and 'tnt' in data_file:
        # df['num_leaf_bottom'] -= 1
    return df
    
def zeros_green_to_nan(data):
    """ Delete notations of septoria on green returning to zero on last notations dates """
    data_pl = data.copy()
    df = pd.DataFrame()
    for lf in set(data_pl['num_leaf_top']):
        data_lf = data_pl[data_pl['num_leaf_top']==lf]
        for pl in set(data_lf['num_plant']):
            df_ = data_lf[data_lf['num_plant']==pl]
            if is_iterable(df_['septo_green']):
                mask = df_['septo_green'].diff()<0
                if mask.iloc[-1] == True:
                    nb_false = [i for i,j in enumerate(mask.values) if j == False][-1]+1
                elif mask.iloc[-1] == False and not all(mask.values == False):
                    nb_false = [i for i,j in enumerate(mask.values) if j == True][-1]
                else:
                    nb_false = len(df_)
                nb_true = len(df_) - nb_false
                mask = [False for f in range(nb_false)] + [True for t in range(nb_true)]
                df_['septo_green'][mask] = np.nan
                df_['septo_green'][df_['necro']==100] = np.nan
                #df = df.append(df_)
            if is_iterable(df_['severity']) and df_['severity'].iloc[-1]==0.:
                mask = df_['severity'].diff()<0
                if mask.iloc[-1] == True:
                    nb_false = [i for i,j in enumerate(mask.values) if j == False][-1]+1
                elif mask.iloc[-1] == False and not all(mask.values == False):
                    nb_false = [i for i,j in enumerate(mask.values) if j == True][-1]
                else:
                    nb_false = len(df_)
                nb_true = len(df_) - nb_false
                mask = [False for f in range(nb_false)] + [True for t in range(nb_true)]
                df_['severity'][mask] = np.nan
            df = df.append(df_)
    return df

def get_leaf_dates(adel, event='tip', nff = None):
    """ Get dates of leaf developmental stage in wheat canopy simulated with adel """
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('num_plant').groupby(level=0).count()['n'] - 1
    df_phenT['nff'] =  [nffs[v] for v in df_phenT['num_plant']]
    df_phenT['n_from_top'] = df_phenT['nff'] - df_phenT['n']
    if nff != None:
        df_phenT = df_phenT[df_phenT['nff']==nff]
    if len(df_phenT)>0:
        dates = [np.mean(df_phenT[df_phenT['n_from_top']==lf][event]) for lf in range(1, int(max(df_phenT['n'])+1))]
        df_dates = pd.DataFrame(dates, index = range(1, len(dates)+1), columns = ['date'])
        df_dates.index.name = 'leaf_number'
        if event == 'ssi':
            df_dates.ix[1:4, 'date'] = df_dates.ix[5, 'date']
            df_dates.ix[5:len(df_dates)-1, 'date'] = df_dates.ix[6:, 'date'].values
        return df_dates
    else:
        raise ValueError("No plant with given NFF in simulation")

def get_leaf_emergence_dates(adel, nff = None):
    """ Get dates of leaf emergence in wheat canopy simulated with adel """
    return get_leaf_dates(adel, event='tip', nff = nff)

def get_leaf_ligulation_dates(adel, nff = None):
    """ Get dates of leaf ligulation in wheat canopy simulated with adel """
    return get_leaf_dates(adel, event='col', nff = nff)
        
def get_leaf_senescence_dates(adel, nff = None):
    """ Get dates of leaf senescence in wheat canopy simulated with adel """
    return get_leaf_dates(adel, event='ssi', nff = nff)

def get_leaf_dates_from_file(filename = 'TTleaf_Tremie13.csv', event = 'emergence', fnl = 12):
    """ Get dates of leaf developmental stage in filename in argument """
    if event == 'emergence':
        event = 'TTem'
    elif event == 'ligulation':
        event = 'TTlig'
    file_path = shared_data(alinea.echap, 'disease_measurements/'+filename)
    dates = pd.read_csv(file_path, sep = ';')
    dates.set_index('ranktop', inplace = True)
    dates.index.names = ['leaf_number']
    if fnl is not None:
        dates = dates[dates['nff'] == fnl][event].to_frame()
    else:
        df = pd.DataFrame(index = set(dates.index))
        for fnl in set(dates['nff']):
            df[fnl] = dates[dates['nff'] == fnl][event]
        dates = df.mean(axis = 1).to_frame()
        dates.index.names = ['leaf_number']
    dates.columns = ['date']
    return dates

def get_fnls_from_file(filename = 'TTleaf_Tremie13.csv'):
    """ Get nffs from file of architecture notations """
    file_path = shared_data(alinea.echap, 'disease_measurements/'+filename)
    dates = pd.read_csv(file_path, sep = ';')
    return set(dates['nff'])

def find_dates_filename(weather):
    """ Find filenames for files with dates of leaf stage for Tremie in 2012 or 2013 """
    if weather.data.index[-1].year == 2012:
        return 'TTleaf_Tremie12.csv'
    elif weather.data.index[-1].year == 2013:
        return 'TTleaf_Tremie13.csv'
    else:
        return None

def get_date_name_event(xaxis):
    if xaxis == 'age_leaf':
        return 'date_emergence_leaf'
    elif xaxis == 'age_leaf_lig':
        return 'date_ligulation_leaf'
    elif xaxis == 'age_leaf_vs_flag_lig':
        return 'date_ligulation_flag_leaf'
    elif xaxis == 'age_leaf_vs_flag_emg':
        return 'date_emergence_flag_leaf'
        
def get_df_dates_xaxis(data, xaxis):
    """ Get dates of leaf stage in argument xaxis from data """
    date_name = get_date_name_event(xaxis)
    df_dates = pd.DataFrame(index = set(data['num_leaf_top']))
    df_dates[xaxis] = map(lambda x: np.mean(data[date_name][data['num_leaf_top']==x]), df_dates.index)
    return df_dates
    
def get_count_one_leaf(df, variable = 'severity', xaxis = 'degree_days',
                          num_leaf = 1, from_top = True):
    """ Get number of notations of argument variable on given leaf over all plants in canopy """
    if from_top == True:
        df = df[df['num_leaf_top'] == num_leaf]
    else:
        df = df[df['num_leaf_bottom'] == num_leaf]
    if xaxis in ['date', 'degree_days']:
        return df.groupby(xaxis).count()[variable]
    elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
        df_count = df.groupby('degree_days').count()[variable]
        df_dates = get_df_dates_xaxis(df, xaxis)
        df_count.index -= df_dates.loc[num_leaf, xaxis].astype(int)
        return df_count
        
def table_count_notations(data, weather, variable = 'septo_green', xaxis = 'date',
                            add_ddays=True, from_top = True):
    """ Get number of notations of argument variable on all leaves over all plants in canopy """
    df = data.copy()
    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'
    if xaxis in ['date', 'degree_days']:
        change_index(df, new_index = [xaxis, num_leaf])
        df.index.rename(['Date', 'Leaf number'], inplace=True)
        df = df.groupby(level=[0,1]).count()
        df = df[variable].unstack()
        df[(np.isnan(df)) | (df==0)]='-'
        if xaxis=='date' and add_ddays == True:
            df = add_index_ddays(df, weather)
        return df
    elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
        dfs = []
        for lf in set(df[num_leaf]):
            df_count_lf = get_count_one_leaf(df, variable = variable, xaxis = xaxis, 
                                                num_leaf = lf, from_top = from_top)
            df_count_lf = df_count_lf.to_frame()
            df_count_lf.columns = [lf]
            df_count_lf[(np.isnan(df_count_lf)) | (df_count_lf==0)]='-'
            dfs.append(df_count_lf)
        return pd.concat(dfs, axis = 1)
        
def age_leaf_notations(data, weather, adel, nff=None):
    """ Get age of leaves at the moment of notations """
    df_notations = table_count_notations(data, weather)
    df_emergence = get_leaf_emergence_dates(adel, nff = nff)
    df_emergence = df_emergence.ix[df_notations.columns,:].T
    for i, row in df_notations.iterrows():
        if '-' in row.values:
            ind_na = row[row == '-'].index
        else:
            ind_na = None
        df_notations.loc[i,:] = np.ones(len(row))*i[1] - df_emergence.values
        if ind_na!=None:
            df_notations.loc[i, ind_na] = '-'
    return df_notations

def add_leaf_dates_to_data(df, adel = None, filename = 'TTleaf_Tremie12.csv',):
    
    def complete_df_with_nearest_fnl(df, fnl, nearest_fnl):
        df['date'] = map(lambda x: x[1] + 10*(x[0]-1)*(fnl - nearest_fnl), 
                                        zip(df.index, df['date']))
        return df
    
    grp_fnls = get_fnl_by_plant(df)

    df['FNL'] = map(lambda x: {vv:k for k,v in grp_fnls.iteritems() for vv in v }[x], df['num_plant'])

    dict_emerg = {}
    dict_lig = {}
    dict_flag_lig = {}
    dict_flag_emg = {}
    #TODO use plantgen extensions HS_fit : TTem_leaf, TTlig_leaf
    for fnl in grp_fnls.iterkeys():
        if filename is None:
            assert adel is not None, ("Provide either an adel canopy mock up or file with emergence and ligulation dates")
            fnls = get_fnls(adel)
            if fnl in fnls:
                dict_emerg[fnl] = get_leaf_emergence_dates(adel, fnl)
                dict_lig[fnl] = get_leaf_ligulation_dates(adel, fnl)
            else:
                nearest_fnl = [i for i in fnls][np.argmin([abs(i-fnl) for i in fnls])]
                df_emg = get_leaf_emergence_dates(adel, nearest_fnl)
                df_lig = get_leaf_ligulation_dates(adel, nearest_fnl)
                dict_emerg[fnl]  = complete_df_with_nearest_fnl(df_emg, fnl, nearest_fnl)
                dict_lig[fnl] = complete_df_with_nearest_fnl(df_lig, fnl, nearest_fnl)
        else:
            fnls = get_fnls_from_file(filename)
            if fnl in fnls:
                dict_emerg[fnl] = get_leaf_dates_from_file(filename = filename, 
                                                            event = 'emergence', fnl = fnl)
                dict_lig[fnl] = get_leaf_dates_from_file(filename = filename, 
                                                            event = 'ligulation', fnl = fnl)
            else:
                fnls = get_fnls_from_file(filename)
                nearest_fnl = [i for i in fnls][np.argmin([abs(i-fnl) for i in fnls])]
                df_emg = get_leaf_dates_from_file(filename = filename, 
                                                    event = 'emergence', fnl = nearest_fnl)
                df_lig = get_leaf_dates_from_file(filename = filename, 
                                                    event = 'ligulation', fnl = nearest_fnl)
                dict_emerg[fnl] = complete_df_with_nearest_fnl(df_emg, fnl, nearest_fnl)
                dict_lig[fnl] = complete_df_with_nearest_fnl(df_lig, fnl, nearest_fnl)

        dict_flag_lig[fnl] = dict_lig[fnl].loc[1, 'date']
        dict_flag_emg[fnl] = dict_emerg[fnl].loc[1, 'date']

    df['date_emergence_leaf'] = map(lambda x : dict_emerg[x[0]].loc[x[1], 'date'],
                                        zip(df['FNL'], df['num_leaf_top']))
    df['age_leaf'] = df['degree_days'] - df['date_emergence_leaf']

    df['date_ligulation_leaf'] = map(lambda x : dict_lig[x[0]].loc[x[1], 'date'],
                                        zip(df['FNL'], df['num_leaf_top']))
    df['age_leaf_lig'] = df['degree_days'] - df['date_ligulation_leaf']

    df['date_ligulation_flag_leaf'] = map(lambda x: dict_flag_lig[x], df['FNL'])
    df['age_leaf_vs_flag_lig'] = df['degree_days'] - df['date_ligulation_flag_leaf']
    
    df['date_emergence_flag_leaf'] = map(lambda x: dict_flag_emg[x], df['FNL'])
    df['age_leaf_vs_flag_emg'] = df['degree_days'] - df['date_emergence_flag_leaf']
    return df
    
def adapt_data(data, weather, adel = None, hide_old_data = False, correct_leaf_number = False):
    """ Complete data after reading """
    df = data.copy()
    
    # Add Degree days
    df['degree_days'] = [int(weather.data['degree_days'][t]) for t in df.index]
    
    # Add age of leaves
    if adel is None:
        filename = find_dates_filename(weather)
    else:
        filename = None
    df = add_leaf_dates_to_data(df, adel=adel, filename = filename)
    
    if correct_leaf_number == True:
        df['num_leaf_bottom'][df['FNL']==14] -= 1
        df['FNL'][df['FNL']==14] -= 1
    
    # If variable is septoria on green : transform zeros to NaN at the end of season
    if weather.data.index[-1].year in [2012, 2013]:
        df = zeros_green_to_nan(df)

    if hide_old_data == True:
        df['septo_green'][df['green']<50] = np.nan
    return df
    
def add_index_ddays(df, weather):
    """ Customize data frames with degree days in multiple index level 1"""
    df['degree_days'] = [int(weather.data['degree_days'][t]) for t in df.index]
    df.set_index('degree_days', append=True, inplace=True)
    df.index.rename(['Date', 'Degree days'], inplace=True)
    return df
    
def get_only_numeric_data(df):
    """ Drop non numeric data in given columns """
    data = df.copy()
    data.drop([col for col in data.columns if sum(isinstance(s,str) for s in data[col])>0],
                inplace=True,axis=1)
    return data
    
def data_reader(year = 2012, variety = 'Tremie12', from_file = 'control', green_on_green = True, 
                wetness_duration_min = 10., temp_min = 0., temp_max = 25., hide_old_data = False,
                with_adel = False, correct_lf_nb_2012_ctrl = True):
    """ Specific reader for the analysis of septoria measurements in Grignon 2011, 2012, 2013 """
    correct_leaf_number = False
    if year == 2011:
        start_date = "2010-10-15"
        end_date = "2011-06-20"
        
        # Read septoria data
        if variety == 'Mercia':
            data_file='mercia_2011.csv'
        elif variety == 'Rht3':
            data_file='rht3_2011.csv'
        else:
            raise ValueError('Variety unknown for this year: try Mercia or Rht3')
        data = data_reader_2011(data_file = data_file)
        
    elif year == 2012:
        start_date = "2011-10-21"
        end_date = "2012-07-18"
        
        # Read septoria data        
        if variety!='Tremie12':
            raise ValueError('Variety unknown for this year: try Tremie12')
        if from_file == 'control':
            file_suffix = 'ctrl'
            correct_leaf_number = correct_lf_nb_2012_ctrl
        elif from_file == 'reference':
            file_suffix = 'ref'
        elif from_file == 'global':
            file_suffix = 'global'
        elif from_file == 'senescence':
            file_suffix = 'nec'
        else:
            raise ValueError("Modality unknown for this year: try 'control' or 'reference' or 'global' or 'senescence'")
        data = data_reader_tremis('tremis_2012_'+file_suffix+'.csv', 
                                  green_on_green = green_on_green)
            
    elif year == 2013:
        start_date = "2012-10-29"
        end_date = "2013-08-01"
        
        # Read septoria data        
        if variety!='Tremie13':
            raise ValueError('Variety unknown for this year: try Tremie13')
        if from_file == 'control':
            file_suffix = 'ctrl'
        elif from_file == 'reference':
            file_suffix = 'ref'
        elif from_file == 'global':
            file_suffix = 'global'
        elif from_file == 'senescence':
            file_suffix = 'nec'
        else:
            raise ValueError("Modality unknown for this year: try 'control' or 'reference' or 'global' or 'senescence'")
        data = data_reader_tremis('tremis_2013_'+file_suffix+'.csv', 
                                  green_on_green = green_on_green)
    
    if year in [2012, 2013] and from_file == 'global':
        data['num_plant'] += 30*(data['bloc']-1)
    elif year in [2012, 2013] and from_file == 'senescence':
        data[data['data_set']=='tagged']['num_plant'] *= 2
    
    # Read weather data
    start, end = format_date(start_date, end_date)
    weather = read_weather(start, end, wetness_duration_min = wetness_duration_min, 
                        temp_min = temp_min, temp_max = temp_max)

    # Generate wheat stand for the year
    if with_adel == True:
        reconst = EchapReconstructions()
        adel = reconst.get_reconstruction(name=variety, nplants = 100, nsect = 7, 
                                            disc_level = 5, aspect = 'line')
    else:
        adel = None
        
    # Adapt septoria data
    data = adapt_data(data, weather, adel, hide_old_data = hide_old_data, 
                        correct_leaf_number = correct_leaf_number)
    
    if with_adel == True:
        return data, weather, adel
    else:
        return data, weather
    
def compare_contingency_by_fnl(data, weather, variable = 'septo_green'):
    """ Generate comparative dataframe with number of notations for given fnl """
    data_grouped = group_by_fnl(data)
    dict_df_count = {}
    for fnl, df in data_grouped.iteritems():
        df_count = table_count_notations(df, weather, variable = variable)
        df_count.columns = pd.MultiIndex.from_product([['FNL %d' %fnl], df_count.columns])
        dict_df_count[fnl] = df_count
    df_empty = pd.DataFrame(['/' for i in range(len(df_count))], index = df_count.index, 
                            columns=pd.MultiIndex.from_product([['/'], ['/']]))
    return pd.concat([pd.concat([df, df_empty], axis=1) for df in dict_df_count.itervalues()], axis = 1)