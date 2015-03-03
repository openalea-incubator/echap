""" Compare results of simulation with disease notations. """

from alinea.alep.disease_outputs import get_recorder, get_mean_by_leaf, mean_by_leaf
from septo_data_reader import *
from weather_reader import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy
plt.ion()

def add_leaf_dates_obs(data, weather):
    pass

def recorder_to_dataframe(recorder, weather = None):
    recos = []
    for pl, rec_pl in recorder.iteritems():
        for lf, rec_lf in rec_pl.iteritems():
            rec_lf.data['plant'] = int(pl[1:])
            rec_lf.data['num_leaf_top'] = int(lf[1:])
            rec_lf.data['num_leaf_ref'] = int(lf[1:])
            rec_lf.data['variety'] = 'tremie'
            rec_lf.data = rec_lf.data.rename(columns = {'date_sequence':'datetime'})
            # Temp : to move in disease outputs
            rec_lf.data['pycnidia_coverage'] = rec_lf.data['ratio_spo'] + rec_lf.data['ratio_empty']
            rec_lf.data['pycnidia_coverage_on_green'] = rec_lf.data['ratio_spo'] + rec_lf.data['ratio_empty']
            # TODO link with weather and date ligulation/emergence to add columns
            if weather is not None:
                rec_lf.data = add_leaf_dates_obs(rec_lf.data, weather)
            recos.append(rec_lf.data)
    return pd.concat(recos)

def get_mean_one_leaf_sim(df, variable = 'severity', xaxis = 'degree_days',
                          num_leaf = 1, from_top = True):
    df=df.rename(columns = {'date_sequence':'date'})
    if from_top == True:
        df = df[df['num_leaf_top'] == num_leaf]
    else:
        df = df[df['num_leaf_bottom'] == num_leaf]
    return df.groupby(xaxis).mean()[variable]

def get_data_obs(year = 2012, variety = 'Tremie12', 
                 from_file = 'control', renumber_leaves = False, 
                 wetness_duration_min = 10., temp_min = 0., temp_max = 25.):
    """ Return data_obs, weather and wheat mock up """
    return data_reader(year = year, variety = variety, 
                        from_file = from_file, renumber_leaves = renumber_leaves, 
                        wetness_duration_min = wetness_duration_min,
                        temp_min = temp_min, temp_max = temp_max)
                        
def get_data_sim(filename = 'recorder_2012_1pl_7sect_frac2e-4.pckl'):
    recorder = get_recorder('../alep/example/tremie/' + filename)
    return recorder_to_dataframe(recorder)

def get_best_rmse(df_mean_obs, df_mean_sim):
    notnull_obs = df_mean_obs[i+1][df_mean_obs[i+1].notnull()]
    nearest_sim = [df_mean_sim.iloc[np.abs(np.array(df_mean_sim.index) - ind).argmin(), i+1] for ind in notnull_obs.index]
    return numpy.sqrt(((nearest_sim - notnull_obs) ** 2).mean())
    
def plot_one_leaf_sim(df_sim, variable = 'severity', xaxis = 'degree_days',
                      num_leaf = 1, from_top = True, ax = None, return_df_mean = True,
                      color = 'r', linestyle = '-', marker = ''):
    if ax is None:
        fig, ax = plt.subplots(1)
    df_mean_sim = get_mean_one_leaf_sim(df = df_sim, variable = variable, xaxis = xaxis,
                                        num_leaf = num_leaf, from_top = from_top)
    ax.plot(df_mean_sim.index, df_mean_sim['F%d' % num_leaf], 
            color = color, linestyle = linestyle, marker = marker)
    if return_df_mean == True:
        return df_mean_sim

def plot_comparison_sim_obs(data_obs, weather, variable='severity', xaxis = 'degree_days',
                            leaves = range(1,6), xlims = [1000, 2200]):
    df_mean_obs, df_low, df_high, fig, axs = plot_confidence_and_boxplot(data_obs, weather, 
                                                                         leaves = leaves,
                                                                         variable = variable,
                                                                         xaxis = xaxis,
                                                                         xlims = xlims,
                                                                         return_fig = True)
    if variable == 'severity':
        variable = 'pycnidia_coverage'
    
    for i, ax in enumerate(axs.flat):
        df_mean_sim = plot_one_leaf_sim(df_sim, variable = variable, 
                                        xaxis = xaxis, num_leaf = i+1, ax = ax)
        ax.annotate('RMSE : %.2f' %get_best_rmse(df_mean_obs, df_mean_sim), xy=(0.05, 0.75), xycoords='axes fraction', fontsize=14)

        
df_mean_sim, df_low, df_high, fig, axs = plot_confidence_and_boxplot(data_sim, weather, leaves = range(1,6), variable = 'pycnidia_coverage_on_green', xaxis = 'degree_days', xlims = [1000, 2200], return_fig = True)
# Get simulation results
# filename = '../alep/example/tremie/recorder_2012_1pl_7sect_frac2e-4.pckl'
# filename = '../alep/example/tremie/recorder_2012_30pl_7sect_frac2e-4.pckl'

# TODO : rename /!\

filename = '../alep/example/tremie/recorder_2012_30pl_7sect_frac2e-4_lat330.pckl'
recorder = get_recorder(filename)
df_mean_sim = get_mean_by_leaf('necrosis_percentage', recorder)*100

for i, ax in enumerate(axs.flat):
    ax.plot(df_mean_sim.index, df_mean_sim['F%d' % (i+1)], 'r')
    notnull_obs = df_mean_obs[i+1][df_mean_obs[i+1].notnull()]
    nearest_sim = [df_mean_sim.iloc[np.abs(np.array(df_mean_sim.index) - ind).argmin(), i+1] for ind in notnull_obs.index]
    rmse = numpy.sqrt(((nearest_sim - notnull_obs) ** 2).mean())
    ax.annotate('RMSE : %.2f' %rmse, xy=(0.05, 0.75), xycoords='axes fraction', fontsize=14)
    
df_mean_sim = get_mean_by_leaf('nb_lesions', recorder)

for i, ax in enumerate(axs.flat):
    ax.plot(df_mean_sim.index, df_mean_sim['F%d' % (i+1)], 'b')
    notnull_obs = df_mean_obs[i+1][df_mean_obs[i+1].notnull()]
    nearest_sim = [df_mean_sim.iloc[np.abs(np.array(df_mean_sim.index) - ind).argmin(), i+1] for ind in notnull_obs.index]
    rmse = numpy.sqrt(((nearest_sim - notnull_obs) ** 2).mean())
    ax.annotate('RMSE : %.2f' %rmse, xy=(0.05, 0.75), xycoords='axes fraction', fontsize=14)
    
df_mean_sim = get_mean_by_leaf('nb_dispersal_units', recorder)

for i, ax in enumerate(axs.flat):
    ax.plot(df_mean_sim.index, df_mean_sim['F%d' % (i+1)], 'k')
    notnull_obs = df_mean_obs[i+1][df_mean_obs[i+1].notnull()]
    nearest_sim = [df_mean_sim.iloc[np.abs(np.array(df_mean_sim.index) - ind).argmin(), i+1] for ind in notnull_obs.index]
    rmse = numpy.sqrt(((nearest_sim - notnull_obs) ** 2).mean())
    ax.annotate('RMSE : %.2f' %rmse, xy=(0.05, 0.75), xycoords='axes fraction', fontsize=14)