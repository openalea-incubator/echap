'''
Created on 25 oct. 2011
Check test results in the created temp directory. 

@author: cchambon
'''

import tempfile
import datetime

import pandas
from openalea.core.path import path

import thermal_time


def setup_module():
    module_vars = globals()
    module_vars['latitude'] = {'HM': 54.11,'RM': 52.13}
    # read data to fit from csv files
    HMj_df_raw = pandas.read_csv(path()/'src/alinea/echap/ADASWheatReconstruction/data/HMj.csv', index_col=0, na_values=['NA'], parse_dates=True)
    RMj_df_raw = pandas.read_csv(path()/'src/alinea/echap/ADASWheatReconstruction/data/RMj.csv', index_col=0, na_values=['NA'], parse_dates=True)
    RMh_00_df_raw = pandas.read_csv(path()/'src/alinea/echap/ADASWheatReconstruction/data/RMh.csv', index_col=0, na_values=['NA'], parse_dates=True)

    start = datetime.datetime(1998, 9, 1)
    end = datetime.datetime(1999, 8, 31)

    def crit_year(index_i):
        # Permits to extract data from `start` to `end`
        if index_i >= start and index_i <= end:
            return True
        return False
    
    module_vars['HMj_99_df'] = HMj_df_raw.select(crit_year)
    module_vars['RMj_99_df'] = RMj_df_raw.select(crit_year)
    
    start = datetime.datetime(1999, 9, 1)
    end = datetime.datetime(2000, 8, 31)
        
    module_vars['HMj_00_df'] = HMj_df_raw.select(crit_year)
    
    # RMj.csv already contains hourly min and max air temperature. 
    # Thus we just have to calculate the mean for each (Tmin,Tmax) tuple.
    Tair_mean = (RMh_00_df_raw['Tmin'] + RMh_00_df_raw['Tmax']) / 2.0
    module_vars['RMh_00_df'] = pandas.DataFrame(Tair_mean, index=RMh_00_df_raw.index, columns=['Tair'])
    

def test_parton_logan():
    module_vars = globals()
    # convert daily data to hourly data
    module_vars['HMh_99_df'] = thermal_time.parton_logan(module_vars['HMj_99_df'], module_vars['latitude']['HM'])
    module_vars['HMh_00_df'] = thermal_time.parton_logan(module_vars['HMj_00_df'], module_vars['latitude']['HM'])
    module_vars['RMh_99_df'] = thermal_time.parton_logan(module_vars['RMj_99_df'], module_vars['latitude']['RM'])
    

def test_dsT():
    module_vars = globals()
    # calculate dsT for hourly data
    module_vars['HMh_99_df']['dsT'] = thermal_time.dsT(module_vars['HMh_99_df']['Tair'].values)
    module_vars['HMh_00_df']['dsT'] = thermal_time.dsT(module_vars['HMh_00_df']['Tair'].values)
    module_vars['RMh_99_df']['dsT'] = thermal_time.dsT(module_vars['RMh_99_df']['Tair'].values)
    module_vars['RMh_00_df']['dsT'] = thermal_time.dsT(module_vars['RMh_00_df']['Tair'].values)


def test_dsTc():
    module_vars = globals()
    # calculate dsTc for hourly data
    module_vars['HMh_99_df']['dsTc'] = thermal_time.dsTc(module_vars['HMh_99_df']['Tair'].values)
    module_vars['HMh_00_df']['dsTc'] = thermal_time.dsTc(module_vars['HMh_00_df']['Tair'].values)
    module_vars['RMh_99_df']['dsTc'] = thermal_time.dsTc(module_vars['RMh_99_df']['Tair'].values)
    module_vars['RMh_00_df']['dsTc'] = thermal_time.dsTc(module_vars['RMh_00_df']['Tair'].values)
    
    
def teardown_module():
    module_vars = globals()
    results_directory = path(tempfile.mkdtemp(suffix='_thermal_time_tests'))
    # write fitted data to csv files in `results_directory`
    module_vars['HMh_99_df'].to_csv(results_directory/'HMh_99.csv', na_rep='NA', index_label='date')
    module_vars['HMh_00_df'].to_csv(results_directory/'HMh_00.csv', na_rep='NA', index_label='date')
    module_vars['RMh_99_df'].to_csv(results_directory/'RMh_99.csv', na_rep='NA', index_label='date')
    module_vars['RMh_00_df'].to_csv(results_directory/'RMh_00.csv', na_rep='NA', index_label='date')
    # you can now check test results in results_directory 

    