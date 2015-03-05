# -*- coding: latin1 -*- 

""" Get simulation outputs of Alep septoria as pandas DataFrames complying with the format generated 
    by readers in 'septo_data_reader' and define comparison methods """

import pandas as pd
import numpy as np
from alinea.alep.disease_outputs import *
from alinea.echap.disease.septo_data_reader import *
from alinea.echap.disease.septo_data_treatment import *
    
# Get data as DataFrames ###########################################################################
def recorder_to_dataframe(recorder, weather = None, adel = None, skipna = True):
    """ Translate recorder object in DataFrame with the same format as disease notations """
    recos = []
    for pl, rec_pl in recorder.iteritems():
        for lf, rec_lf in rec_pl.iteritems():
            rec_lf.data['plant'] = int(pl[1:])
            rec_lf.data['num_leaf_top'] = int(lf[1:])
            rec_lf.data['num_leaf_bottom'] = len(rec_pl) - int(lf[1:]) + 1
            rec_lf.data['variety'] = 'tremie'
            rec_lf.data = rec_lf.data.rename(columns = {'date_sequence':'datetime'})
            # Temp : to move in disease outputs
            rec_lf.data['pycnidia_coverage'] = (rec_lf.data['ratio_spo'] + rec_lf.data['ratio_empty'])*100
            rec_lf.data['pycnidia_coverage_on_green'] = (rec_lf.data['ratio_spo'] + rec_lf.data['ratio_empty'])*100
            recos.append(rec_lf.data)
    data_sim = pd.concat(recos)
    
    # Convert ratios in percentage
    list_of_ratios = [var for var in data_sim.columns if var.startswith('ratio')]+['severity', 'necrosis_percentage']
    data_sim[list_of_ratios] = data_sim[list_of_ratios].apply(lambda x: x*100.)
    
    # Add leaf dates
    if weather is not None:
        filename = find_dates_filename(weather)
        data_sim = add_leaf_dates_to_data(data_sim, adel, filename = filename)
    
    # Ignore data from dates with dead leaves in each layer
    if skipna == True:
        df_count = table_count_notations(data_sim, weather, variable = 'severity', add_ddays = True)
        for lf in df_count.columns:
            df_lf = df_count[lf][map(lambda x: isinstance(x, (int, float)), df_count[lf])]
            nan_dates = df_lf[df_lf<df_lf.max()].reset_index().loc[:,'Date']
            if len(nan_dates)>0:
                for variable in data_sim.columns:
                    if variable not in ['datetime', 'degree_days', 'date_death',
                                        'variety', 'plant', 'num_leaf_top', 'num_leaf_bottom']:
                        data_sim[variable][(data_sim['num_leaf_top']==lf) & 
                                (data_sim['datetime'].isin([d for d in nan_dates]))] = np.nan
    return data_sim
    
def recorder_to_dataframe_single_variable(recorder, weather = None, adel = None, skipna = True):
    """ Translate recorder object in DataFrame with the same format as disease notations """
    to_keep = ['datetime', 'degree_days', 'date_death','variety',
                'plant', 'num_leaf_top', 'num_leaf_bottom'] + ['pycnidia_coverage']
    recos = []
    for pl, rec_pl in recorder.iteritems():
        for lf, rec_lf in rec_pl.iteritems():
            rec_lf.data['plant'] = int(pl[1:])
            rec_lf.data['num_leaf_top'] = int(lf[1:])
            rec_lf.data['num_leaf_bottom'] = len(rec_pl) - int(lf[1:]) + 1
            rec_lf.data['variety'] = 'tremie'
            rec_lf.data = rec_lf.data.rename(columns = {'date_sequence':'datetime'})
            # Temp : to move in disease outputs
            rec_lf.data['pycnidia_coverage'] = (rec_lf.data['ratio_spo'] + rec_lf.data['ratio_empty'])*100
            rec_lf.data.drop([c for c in rec_lf.data.columns if not c in to_keep],inplace=True,axis=1)
            recos.append(rec_lf.data)
    data_sim = pd.concat(recos)
    
    # Convert ratios in percentage
    # list_of_ratios = [var for var in data_sim.columns if var.startswith('ratio')]+['severity', 'necrosis_percentage']
    # data_sim[list_of_ratios] = data_sim[list_of_ratios].apply(lambda x: x*100.)
    
    # Add leaf dates
    if weather is not None:
        filename = find_dates_filename(weather)
        data_sim = add_leaf_dates_to_data(data_sim, adel, filename = filename)
    
    # Ignore data from dates with dead leaves in each layer
    if skipna == True:
        df_count = table_count_notations(data_sim, weather, variable = 'pycnidia_coverage', add_ddays = True)
        for lf in df_count.columns:
            df_lf = df_count[lf][map(lambda x: isinstance(x, (int, float, numpy.integer)), df_count[lf])]
            nan_dates = df_lf[df_lf<df_lf.max()].reset_index().loc[:,'Date']
            if len(nan_dates)>0:
                for variable in data_sim.columns:
                    if variable not in ['datetime', 'degree_days', 'date_death',
                                        'variety', 'plant', 'num_leaf_top', 'num_leaf_bottom']:
                        data_sim[variable][(data_sim['num_leaf_top']==lf) & 
                                (data_sim['datetime'].isin([d for d in nan_dates]))] = np.nan
    return data_sim
    
def get_data_sim(filename = '../alep/example/tremie/recorder_2012_30pl_7sect_frac2e-4.pckl',
                 weather = None, adel =None, skipna = True):
    recorder = get_recorder(filename)
    return recorder_to_dataframe(recorder, weather = weather, adel = adel, skipna = skipna)

# Data manipulation on results of simulation #######################################################
def get_mean_one_leaf_sim(df, variable = 'severity', xaxis = 'degree_days',
                          num_leaf = 1, from_top = True):
    """ Get average of input variable on input leaf over all plants in simulated canopy """
    if from_top == True:
        df = df[df['num_leaf_top'] == num_leaf]
    else:
        df = df[df['num_leaf_bottom'] == num_leaf]
    if xaxis in ['datetime', 'degree_days']:
        return df.groupby(xaxis).mean()[variable]
    elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig']:
        df_mean = df.groupby('degree_days').mean()[variable]
        df_dates = get_df_dates_xaxis(df, xaxis)
        df_mean.index -= df_dates.loc[num_leaf, xaxis]
        return df_mean

def get_df_mean_sim(data_sim, variable = 'severity', xaxis = 'degree_days', leaves = range(1,14)):
    """ Get average for given variable over xaxis for all given leaves """
    dfs = []
    for lf in leaves:
        df_mean_lf = get_mean_one_leaf_sim(data_sim, variable = variable, 
                                            xaxis = xaxis, num_leaf = lf, from_top = True)
        df_mean_lf = df_mean_lf.to_frame()
        df_mean_lf.columns = [lf]
        dfs.append(df_mean_lf)
    return pd.concat(dfs, axis = 1)
    
def get_mean_rmse(df_mean_obs, df_mean_sim, num_leaf):
    """ Get mean RMSE between observation data and simulation results """
    notnull_obs = df_mean_obs[num_leaf][df_mean_obs[num_leaf].notnull()]
    nearest_sim = [df_mean_sim[num_leaf].iloc[np.abs(np.array(df_mean_sim.index) - ind).argmin()] 
                    for ind in notnull_obs.index]
    return numpy.sqrt(((nearest_sim - notnull_obs) ** 2).mean())
    
# Compare observations and simulations #############################################################
def plot_one_leaf_sim(data_sim, variable = 'severity', xaxis = 'degree_days',
                      num_leaf = 1, from_top = True, ax = None, return_df_mean = True,
                      color = 'r', linestyle = '-', marker = ''):
    """ Plot input variable versus input axis """
    if ax is None:
        fig, ax = plt.subplots(1)
    df_mean_sim = get_mean_one_leaf_sim(df = data_sim, variable = variable, xaxis = xaxis,
                                        num_leaf = num_leaf, from_top = from_top)
    ax.plot(df_mean_sim.index, df_mean_sim, 
            color = color, linestyle = linestyle, marker = marker)
    if return_df_mean == True:
        return df_mean_sim

def plot_mean_by_leaf_sim(data_sim, variable = 'severity', leaves = range(1,14), 
                          xaxis = 'degree_days', ax = None, xlims = None,
                          ylims = None, ylabel = 'Severity', xlabel = 'Degree Days',
                          fig_size = (8,6), grid = True):
    """ Plot average of given variable on given separate leaves for all plants """    
    df_mean = get_df_mean_sim(data_sim, variable = variable, xaxis = xaxis, leaves = leaves)
    if ax == None:
        fig, ax = plt.subplots(1, figsize = fig_size)
    for leaf in df_mean.columns:
        x_data = df_mean.index[~np.isnan(df_mean.ix[:,leaf])]
        ax.plot(x_data, df_mean.loc[:, leaf][~np.isnan(df_mean.ix[:,leaf])])
    
    # Customize plot
    if ylims is not None:
        ax.set_ylim(ylims)
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize = 18)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize = 18)
    if grid == True:
        ax.grid()
    ax.legend(leaves, title = 'Leaf\nnumber', loc='center left', bbox_to_anchor=(1, 0.5))
    
def plot_comparison_confidence_and_boxplot_sim_obs(data_obs, data_sim,
                                                    weather, variable='severity', 
                                                    xaxis = 'degree_days', 
                                                    leaves = range(1,7), 
                                                    xlims = [1000, 2200], 
                                                    display_confidence = False, 
                                                    display_rmse = False):
    df_mean_obs, df_low, df_high, fig, axs = plot_confidence_and_boxplot(data_obs, weather, leaves = leaves,
                                                                         variable = variable, xaxis = xaxis,
                                                                         xlims = xlims, return_fig = True)
    if variable == 'severity':
        variable = 'pycnidia_coverage'
    
    if display_confidence == True:
        df_mean_sim, df_low, df_high, fig, axs = plot_confidence_and_boxplot(data_sim, weather, leaves = leaves,
                                                                             variable = variable, xaxis = xaxis,
                                                                             xlims = xlims, marker ='', fixed_color = 'r',
                                                                             fig = fig, axs = axs,
                                                                             display_box = False, return_fig = True)
    else:
        df_mean_sim = get_df_mean_sim(data_sim, variable = variable, xaxis = xaxis, leaves = leaves)
        for i, ax in enumerate(axs.flat):
            leaf = i+1
            x_data = df_mean_sim.index[~np.isnan(df_mean_sim.ix[:,leaf])]
            ax.plot(x_data, df_mean_sim.loc[:, leaf][~np.isnan(df_mean_sim.ix[:,leaf])], 'r')
            if display_rmse == True:
                ax.annotate('RMSE : %.2f' %get_mean_rmse(df_mean_obs, df_mean_sim, num_leaf = leaf), xy=(0.05, 0.75),
                            xycoords='axes fraction', fontsize=14)