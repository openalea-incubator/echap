""" Reconstruct a MTG wheat canopy at different ages and save leaf properties (area green, 
    area senesced) and compare it with field data """

import pandas
import numpy
import datetime
import matplotlib.pyplot as plt

from alinea.echap.architectural_reconstructions import echap_reconstructions, HS_fit
import alinea.echap.architectural_data as archidb
from alinea.adel.newmtg import adel_labels
from alinea.astk.plantgl_utils import get_height
from alinea.echap.weather_data import *
from alinea.astk.TimeControl import *
from alinea.echap.evaluation_canopy_properties import plot_mean, plot_sum

# Run simu and save leaf properties ################################################################
class AdelLeafRecorder:
    """ Record simulation output in a dataframe during simulation """
    def __init__(self, variety = 'Tremie12'):
        self.variety = variety
        self.data = pandas.DataFrame(columns = ['date',
                                                'degree_days',
                                                'num_plant',
                                                'axis',
                                                'num_leaf_bottom',
                                                'area', 
                                                'green_area',
                                                'senesced_area',
                                                'length',
                                                'green_length',
                                                'senesced_length',
                                                'height',
                                                'min_height',
                                                'max_height'])
    
    def get_values_single_leaf(self, g, date, degree_days, id_list):
        areas = g.property('area')
        green_areas = g.property('green_area')
        senesced_areas = g.property('senesced_area')
        lengths = g.property('length')
        green_lengths = g.property('green_length')
        senesced_lengths = g.property('senesced_length')
        geometries = g.property('geometry')
        
        dict_lf = {}
        dict_lf['date'] = date
        dict_lf['degree_days'] = degree_days
        a_label_splitted = self.a_labels[id_list[0]].split('_')
        dict_lf['num_plant'] = int(a_label_splitted[0].split('plant')[1])
        dict_lf['axis'] = a_label_splitted[1]
        dict_lf['num_leaf_bottom'] = int(a_label_splitted[2].split('metamer')[1])
        dict_lf['area'] = sum([areas[id] for id in id_list])
        dict_lf['green_area'] = sum([green_areas[id] for id in id_list])
        dict_lf['senesced_area'] = sum([senesced_areas[id] for id in id_list])
        dict_lf['length'] = sum([lengths[id] for id in id_list])
        dict_lf['green_length'] = sum([green_lengths[id] for id in id_list])
        dict_lf['senesced_length'] = sum([senesced_lengths[id] for id in id_list])
        heights = [numpy.mean(get_height({id:geometries[id]}).values()) for id in id_list]
        dict_lf['height'] = numpy.mean(heights) if len(heights)>0. else 0.
        dict_lf['min_height'] = min(heights) if len(heights)>0. else 0.
        dict_lf['max_height'] = max(heights) if len(heights)>0. else 0.
        return dict_lf
        
    def record(self, g, date = None, degree_days = None):
        self.a_labels = {vid:lab for vid, lab in adel_labels(g, scale = 5).iteritems() 
                            if 'LeafElement' in lab}
        v_length = g.property('visible_length')
        labels = g.property('label')
        geometries = g.property('geometry')
        areas = g.property('area')
        blades = [id for id,lb in labels.iteritems() if lb.startswith('blade') and v_length[id]>0]
        for blade in blades:
            id_list = [id for id in g.components(blade) if geometries.get(id) is not None 
                                                        and areas.get(id) is not None
                                                        and labels[id].startswith('LeafElement')]
            if len(id_list)>0:
                dict_lf = self.get_values_single_leaf(g = g, date = date, 
                                                      degree_days = degree_days, 
                                                      id_list = id_list)
                self.data = self.data.append(dict_lf, ignore_index = True)
                
    def add_leaf_numbers(self):
        for pl in set(self.data['num_plant']):
            df_pl = self.data[self.data['num_plant'] == pl]
            for ax in set(df_pl['axis']):
                df_ax = df_pl[df_pl['axis'] == ax]
                fnl = df_ax['num_leaf_bottom'].max()
                df_ax.loc[:, 'fnl'] = fnl
                df_ax.loc[:, 'num_leaf_top'] = fnl - df_ax['num_leaf_bottom'] + 1
                self.data.loc[df_ax.index, 'fnl'] = fnl
                self.data.loc[df_ax.index, 'num_leaf_top'] = fnl - df_ax['num_leaf_bottom'] + 1
                for date in set(df_ax['degree_days']):
                    df_date = df_ax[df_ax['degree_days'] == date]
                    current_max_bottom = df_date['num_leaf_bottom'].max()
                    cur_max_leaf_top = df_date['fnl'] - current_max_bottom + 1
                    self.data.loc[df_date.index, 'cur_max_leaf_top'] = cur_max_leaf_top
                    self.data.loc[df_date.index, 'cur_num_leaf_top'] = cur_max_leaf_top - df_date['num_leaf_top'] + 1
    
    def add_haun_stage(self):
        HSconv = HS_fit()[self.variety]
        self.data['HS'] = HSconv(self.data['degree_days'])
    
    def post_treatment(self):
        self.add_leaf_numbers()
        self.add_haun_stage()
    
    def save_mean_and_sum_by_axis(self, filename):
        def get_mean_data(variable):
            df_sum.loc[:, variable] = df_mean.loc[:, variable]
        
        df_mean = self.data.groupby(['date', 'num_plant', 'axis']).mean()
        df_mean = df_mean.reset_index()
        df_mean.loc[:, 'degree_days'] = map(lambda x: numpy.round(x),
                                            df_mean.loc[:, 'degree_days'])
        df_mean.loc[:, 'HS'] = map(lambda x: numpy.round(x, 2), df_mean.loc[:, 'HS'])
        df_mean.to_csv(filename[:-4]+'_mean'+filename[-4:], index = False)
        
        df_sum = self.data.groupby(['date', 'num_plant', 'axis']).sum()
        df_sum = df_sum.reset_index()
        for var in ['degree_days', 'HS', 'num_leaf_bottom', 'num_leaf_top', 'fnl', 
                    'cur_max_leaf_top', 'cur_num_leaf_top']:
            get_mean_data(var)
        df_sum.to_csv(filename[:-4]+'_sum'+filename[-4:], index = False)
        
    def save(self, filename):
        self.post_treatment()
        self.data.to_csv(filename, index = False)
        self.save_mean_and_sum_by_axis(filename)

def get_file_name(variety = 'Tremie12', nplants = 30, group = None):
    if group is None:
        return 'leaf_'+variety.lower() + '_' + str(nplants) + 'pl.csv'
    else:
        assert group in ['sum', 'mean'], ValueError("group unknown, try 'sum' or 'mean'")
        return 'leaf_'+variety.lower() + '_' + str(nplants) + 'pl_' + group + '.csv'
    
def get_file_path(variety = 'Tremie12', nplants = 30, group = None):
    filename = get_file_name(variety = variety, nplants = nplants, group = group)
    return str(shared_data(alinea.echap)/'architectural_simulations'/filename)

def add_notation_dates(weather_data, variety):
    if variety in ['Mercia', 'Rht3']:
        df = pandas.DataFrame(numpy.zeros(len(weather_data)), index = weather_data.index)
    else:
        dates = map(lambda x: pandas.to_datetime(x, dayfirst=True), archidb.scan_dates(variety = variety))
        df = pandas.DataFrame(numpy.zeros(len(weather_data)), index = weather_data.index)
        df.loc[dates, :] = 1
    return df.values == 1
    
def ddays_and_notations_filter(seq, weather, delay = 20.):
    notations = weather.data.notation_dates[seq]
    ddays = weather.data.degree_days[seq].diff()
    ddays.ix[0] = 0.
    df = notations.copy()
    count_ddays = 0.
    for i, row in df.iteritems():
        if row == True:
            count_ddays = 0.
        count_ddays += ddays[i]
        if count_ddays >= delay and i > df.index[0]:
            df[i - 1] = True
            count_ddays = 0.
    df[df.index[0]] = True
    return df
    
def set_adel_canopy(variety = 'Tremie12', nplants = 30, age = 0., 
                    reset = False, reset_data = False):
    reconst = echap_reconstructions(reset=reset, reset_data=reset_data)
    return reconst.get_reconstruction(name = variety, nplants = nplants)
    
def get_weather(variety = 'Tremie12'):
    year = get_year_for_variety(variety = variety)
    start, end = get_start_end_dates(year = year)
    weather = read_weather(start, end)
    weather.check(varnames=['notation_dates'], models={'notation_dates':add_notation_dates}, variety = variety)
    seq = pandas.date_range(start = start, end = end, freq='H')
    return weather, seq
    
def run_and_save(variety = 'Tremie12', nplants = 30, delay = 20., 
                reset = False, reset_data = False):
    # Initialize wheat canopy
    adel = set_adel_canopy(variety = variety, nplants = nplants, age = 0., 
                            reset = reset, reset_data = reset_data)
    recorder = AdelLeafRecorder()
    
    # Get weather
    weather, seq = get_weather(variety = variety)
    
    # Schedule of calls for each model
    time_filter = ddays_and_notations_filter(seq, weather, delay = delay)
    canopy_timing = IterWithDelays(*time_control(seq, time_filter, weather.data))
    
    # Run annual loop and save outputs
    for i, canopy_iter in enumerate(canopy_timing):
        if canopy_iter:
            print canopy_iter.value.index[0]
            g = adel.setup_canopy(age = canopy_iter.value.degree_days[0])
            recorder.record(g, date = canopy_iter.value.index[0], 
                            degree_days = canopy_iter.value.degree_days[0])
    
    filename = get_file_path(variety = variety, nplants = nplants)
    recorder.save(filename = filename)

def run_and_save_all_varieties(nplants = 30, delay = 20., 
                                reset = False, reset_data = False):
    for variety in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']:
        run_and_save(variety=variety, nplants=nplants, delay=delay,
                      reset=reset, reset_data=reset_data)
    
def get_simu_results(variety = 'Tremie12', nplants = 30, group = None):
    file_path = get_file_path(variety = variety, nplants = nplants, group = group)
    df = pandas.read_csv(file_path)
    df['necro'] = 100*df['senesced_area']/df['area'].replace({ 0 : numpy.inf })
    df['necro_tot'] = 100*df['senesced_area']/df['area'].replace({ 0 : numpy.inf })
    return df
    
def plot_by_leaf(data, variable = 'green_area', xaxis = 'degree_days', 
                  leaves = range(1, 14), from_top = True, plant_axis = ['MS'],
                  error_bars = False, error_method = 'confidence_interval', 
                  marker = '', empty_marker = False, linestyle = '-', fixed_color = None, 
                  title = None, legend = True, xlabel = None, ylabel = None,
                  xlims = None, ylims = None, ax = None, return_ax = False, fig_size = (10,8)):
    df = data.copy()
    if ax == None:
        fig, ax = plt.subplots(figsize = fig_size)
    colors = ax._get_lines.set_color_cycle()
    colors = ax._get_lines.color_cycle
        
    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'
    
    proxy = []
    labels = []
    for lf in leaves:
        df_lf = df[(df['axis'].isin(plant_axis)) & (df[num_leaf]==lf)]
        if fixed_color == None:
            color = next(colors)
        else:
            color = fixed_color
        plot_mean(df_lf, variable = variable, xaxis = xaxis, 
                  error_bars = error_bars, error_method = error_method, 
                  marker = marker, empty_marker = empty_marker, linestyle = linestyle, 
                  color = color, title = title, xlabel = xlabel, ylabel = ylabel,
                  xlims = xlims, ylims = ylims, ax = ax)
        proxy += [plt.Line2D((0,1),(0,0), color = color, linestyle ='-')]
        labels += ['L%d' %lf]
    
    if legend == True:
        colors = ax._get_lines.set_color_cycle()
        colors = ax._get_lines.color_cycle
        ax.legend(proxy, labels, title = 'Leaf\nNumber',
                    loc='center left', bbox_to_anchor=(1, 0.5))
    if return_ax == True:
        return ax
    
def plot_by_leaf_by_fnl(data, variable = 'green_area', xaxis = 'degree_days', 
                  leaves = range(1, 14), from_top = True,
                  error_bars = False, error_method = 'confidence_interval', 
                  markers = ['', '', ''], linestyle = ['-', '--', ':'], colors = [None, None, None], 
                  title = None, xlabel = None, ylabel = None,
                  xlims = None, ylims = None, ax = None, return_ax = False, fig_size = (10,8)):
    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'
    df = data[data['axis'] == 'MS']
    if ax is None:
        fig, ax = plt.subplots(figsize = fig_size)
    colors, markers, linestyles = map(lambda x: iter(x), (colors, markers, linestyles))
    markers = iter()
    labels = []
    for fnl in set(df['fnl']):
        df_fnl = df[df['fnl'] == fnl]
        (color, marker, linestyle) = map(lambda x: next(x), (colors, markers, linestyles))
        plot_by_leaf(df_fnl, variable = variable, xaxis = xaxis, 
                      error_bars = error_bars, error_method = error_method, 
                      marker = marker, linestyle = linestyle, fixed_color = color,
                      title = title, xlabel = xlabel, ylabel = ylabel,
                      xlims = xlims, ylims = ylims, ax = ax)
        labels += ['L '+str(lf)+': FNL '+str(fnl) for lf in leaves]
                  
    # Customize legend
    color_list = [next(ax._get_lines.color_cycle) for lf in leaves]
    leg = ax.legend(ax._get_legend_handles(), labels, loc='center left', bbox_to_anchor=(1, 0.5))
    labels, handles = zip(*sorted(zip(labels, ax._get_legend_handles()), 
                            key=lambda t: color_list.index(t[1].get_c())))
    leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    
# Read scan data ###################################################################################
def add_hs_and_degree_days(df, variety = 'Tremie12'):
    weather = read_weather_variety(variety = variety)
    df['degree_days'] = weather.data.degree_days[df['date']].values
    HSconv = HS_fit()[variety]
    df['HS'] = HSconv(df['degree_days'])

def get_archi_tagged_data(variety = 'Tremie12'):
    df = archidb.treated_archi_tagged_data(variety = variety)
    if variety == 'Tremie13':
        df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x, dayfirst = True))
    else:
        df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x))
    add_hs_and_degree_days(df, variety = variety)
    return df
    
def get_symptom_tagged_data(variety = 'Tremie12'):
    df = archidb.treated_symptom_tagged_data(variety = variety)
    df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x))
    add_hs_and_degree_days(df, variety = variety)
    return df

def get_scan_dimensions_single_date(variety = 'Tremie12', date = '09/05/2012'):
    df_all = archidb.scan_dimensions_single_date(variety = variety, date = date)
    df = df_all.loc[:, ['prelevement', 'plant', 'id_Axe', 'rank', 'lmax', 'wmax',
                        'A_bl', 'A_bl_green', 'stat']]
    df = df.rename(columns = {'prelevement':'date', 'plant':'num_plant', 'id_Axe':'axis', 
                              'rank':'num_leaf_bottom', 'lmax':'length', 'wmax':'width',
                              'A_bl':'area', 'A_bl_green':'area_green'})
    df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x, dayfirst=True))
    df['axis'] = df['axis'].apply(lambda x: 'MS' if x == 'MB' else 'T')
    return df
    
def get_scan_dimensions(variety = 'Tremie12'):
    df = pandas.concat([get_scan_dimensions_single_date(variety = variety, date = d) 
                        for d in archidb.scan_dates(variety)])
    add_hs_and_degree_days(df, variety = variety)
    return df.reset_index(drop = True)

def plot_sim_obs_MS(df_sim, df_obs, variable = 'length', xaxis = 'degree_days', 
                     leaves = range(7, 14), from_top = False, error_bars = [False, True],
                     markers = ['', 'o'], linestyles = ['-', '--'], 
                     xlims = None, ylims = None, legend = True, 
                     title = None, ax = None, fig_size = None):
    if ax is None:
        if fig_size is not None:
            fig, ax = plt.subplots(figsize = fig_size)
        else:
            fig, ax = plt.subplots()
    plot_by_leaf(df_sim, variable = variable, xaxis = xaxis, leaves = leaves, 
                      plant_axis = ['MS'], from_top = from_top, error_bars = error_bars[0], 
                      marker = markers[0], linestyle = linestyles[0], title = title,
                      xlabel = xaxis, ylabel = variable, xlims = xlims, 
                      ylims = ylims, legend = legend, ax = ax)
    plot_by_leaf(df_obs, variable = variable, xaxis = xaxis, leaves = leaves, 
                 plant_axis = ['MS'], from_top = from_top, error_bars = error_bars[1],
                 marker = markers[1], linestyle = linestyles[1], title = title,
                 xlabel = xaxis, ylabel = variable, xlims = xlims, ylims = ylims,
                 legend = legend, ax = ax)
    
def compare_sim_scan(variety = 'Tremie12', nplants = 30, variable = 'length', xaxis = 'degree_days', 
                     leaves = range(7, 14), from_top = False, exclude_stat = None, 
                     xlims = None, ylims = None, legend = True, ax = None, df_sim = None, 
                     fig_size = (10, 8)):
    if df_sim is None:
        df_sim = get_simu_results(variety = variety, nplants = nplants)
    df_obs = get_scan_dimensions(variety = variety)
    if exclude_stat is not None:
        df_obs = df_obs[df_obs['stat'] < exclude_stat]
    plot_sim_obs_MS(df_sim, df_obs, variable = variable, xaxis = xaxis, leaves = leaves, 
                      from_top = from_top, title = 'Simu vs. scans',
                      xlims = xlims, ylims = ylims, legend = legend, ax = ax, fig_size = fig_size)

def compare_sim_archi_tagged(variety = 'Tremie12', nplants = 30,
                            variable = 'length', xaxis = 'degree_days', 
                            leaves = range(8, 14), from_top = False, xlims = None, ylims = None, 
                            legend = True, ax = None, df_sim = None, fig_size = (10, 8)):
    if df_sim is None:
        df_sim = get_simu_results(variety = variety, nplants = nplants)
    df_obs = get_archi_tagged_data(variety = variety)
    plot_sim_obs_MS(df_sim, df_obs, variable = variable, xaxis = xaxis,
                    leaves = leaves, from_top = from_top,
                    title = 'Simu vs. treated archi tagged',
                    xlims = xlims, ylims = ylims, legend = legend, ax = ax, fig_size = fig_size)
                      
def compare_sim_symptom_tagged(variety = 'Tremie12',  nplants = 30, 
                               variable = 'necro', xaxis = 'degree_days', leaves = range(8, 14), 
                               from_top = False, xlims = None, ylims = None, 
                               legend = True, ax = None, df_sim = None, fig_size = (10, 8)):
    if df_sim is None:
        df_sim = get_simu_results(variety = variety, nplants = nplants)
    df_obs = get_symptom_tagged_data(variety = variety)
    plot_sim_obs_MS(df_sim, df_obs, variable = variable, xaxis = xaxis, 
                      leaves = leaves, from_top = from_top,
                      title = 'Simu vs. treated treated symptom tagged',
                      xlims = xlims, ylims = ylims, legend = legend, ax = ax, fig_size = fig_size)

def plot_sim_tagged(df_sim, df_obs_archi, df_obs_symptom,
                    variety = 'Tremie12', variable = 'necro', xaxis = 'degree_days', 
                    leaves = range(8, 14), from_top = False, error_bars = True,
                    xlims = None, ylims = None, 
                    legend = True, ax = None, fig_size = (10, 8)):
    if ax is None:
        if fig_size is not None:
            fig, ax = plt.subplots(figsize = fig_size)
        else:
            fig, ax = plt.subplots()
    plot_by_leaf(df_sim, variable = variable, xaxis = xaxis, leaves = leaves, 
                  plant_axis = ['MS'], from_top = from_top, error_bars = False, 
                  marker = '', linestyle = '-', title = None,
                  xlabel = xaxis, ylabel = variable, xlims = xlims, 
                  ylims = ylims, legend = legend, ax = ax)
    plot_by_leaf(df_obs_archi, variable = variable, xaxis = xaxis, leaves = leaves, 
                 plant_axis = ['MS'], from_top = from_top, error_bars = error_bars,
                 marker = 'o', empty_marker = True, linestyle = '', title = None,
                 xlabel = xaxis, ylabel = variable, xlims = xlims, ylims = ylims,
                 legend = False, ax = ax)
    plot_by_leaf(df_obs_symptom, variable = variable, xaxis = xaxis, leaves = leaves, 
                 plant_axis = ['MS'], from_top = from_top, error_bars = error_bars,
                 marker = 'o', linestyle = '--', title = "Simu vs. all tagged",
                 xlabel = xaxis, ylabel = variable, xlims = xlims, ylims = ylims,
                 legend = False, ax = ax)
    
    labels = ['L%d' %lf for lf in leaves]
    leg1 = ax.legend(ax._get_legend_handles(), labels, title = 'Leaf\nNumber', loc='center left', bbox_to_anchor=(0, 0.65))
    
    proxy = [plt.Line2D((0,1),(0,0), color = 'k', linestyle ='-'),
             plt.Line2D((0,1),(0,0), marker = 'o', markerfacecolor = 'none', markeredgecolor = 'k', linestyle =''),
             plt.Line2D((0,1),(0,0), marker = 'o', color = 'k', linestyle ='--')]
    labels = ['simulation', 'archi tagged', 'symptom tagged']
    leg2 = plt.legend(proxy, labels, loc='best')
    ax.add_artist(leg1)
    ax.add_artist(leg2)
                      
def compare_sim_tagged(variety = 'Tremie12', nplants=30, variable = 'necro', xaxis = 'degree_days', 
                        leaves = range(8, 14), from_top = False, error_bars = True, 
                        xlims = None, ylims = None, legend = True, ax = None, 
                        df_sim = None, fig_size = (10, 8)):
    if df_sim is None:
        df_sim = get_simu_results(variety = variety, nplants = nplants)
    df_obs_archi = get_archi_tagged_data(variety = variety)
    df_obs_symptom = get_symptom_tagged_data(variety = variety)
    plot_sim_tagged(df_sim, df_obs_archi, df_obs_symptom,
                    variety = variety, variable = variable, xaxis = xaxis, 
                    leaves = leaves, from_top = from_top, error_bars = error_bars,
                    xlims = xlims, ylims = ylims, 
                    legend = legend, ax = ax)
                      
def compare_sim_tagged_by_fnl(variety = 'Tremie12', nplants=30, variable = 'necro', xaxis = 'degree_days', 
                        leaves = range(8, 14), from_top = False, error_bars = True, 
                        xlims = None, ylims = None, 
                        legend = True, ax = None, df_sim = None, fig_size = (10, 17)):
    if df_sim is None:
        df_sim = get_simu_results(variety = variety, nplants = nplants)
    df_obs_archi = get_archi_tagged_data(variety = variety)
    df_obs_symptom = get_symptom_tagged_data(variety = variety)
    fnls = set(df_sim['fnl']) & set(df_obs_symptom['fnl']) & set(df_obs_archi['fnl'])
    fig, axs = plt.subplots(len(fnls), 1, figsize = fig_size)
    for i, fnl in enumerate(fnls):
        df_sim_fnl = df_sim[df_sim['fnl'] == fnl]
        df_obs_archi_fnl = df_obs_archi[df_obs_archi['fnl'] == fnl]
        df_obs_symptom_fnl = df_obs_symptom[df_obs_symptom['fnl'] == fnl]
        plot_sim_tagged(df_sim_fnl, df_obs_archi_fnl, df_obs_symptom_fnl,
                        variety = variety, variable = variable, xaxis = xaxis, 
                        leaves = leaves, from_top = from_top, error_bars = error_bars,
                        xlims = xlims, ylims = ylims, 
                        legend = legend, ax = axs[i], fig_size = fig_size)
        axs[i].annotate('FNL %d' % fnl, xy=(0.05, 0.4), xycoords='axes fraction', fontsize=18)
        
def compare_sim_scan_archi_tagged(variety = 'Tremie12', nplants=30, variable = 'length', xaxis = 'degree_days', 
                            leaves = range(8, 14), exclude_stat = None,
                            xlims = None, ylims = None, fig_size = (16, 10)):
    df_sim = get_simu_results(variety = variety, nplants = nplants)
    fig, axs = plt.subplots(1, 2, figsize = fig_size)
    compare_sim_scan(variety = variety, variable = variable, xaxis = xaxis, 
                     leaves = leaves, exclude_stat = exclude_stat, 
                     xlims = xlims, ylims = ylims, legend = False,
                     ax = axs[0], df_sim = df_sim)
    compare_sim_archi_tagged(variety = variety, variable = variable, xaxis = xaxis, 
                           leaves = leaves, xlims = xlims, ylims = ylims,
                           ax = axs[1], df_sim = df_sim)
    letters = iter(['a', 'b', 'c', 'd'])
    for ax in axs:
        ax.annotate(next(letters), xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
        
def compare_sim_scan_symptom_tagged(variety = 'Tremie12', nplants=30,
                                    variable = 'length', xaxis = 'degree_days', 
                                    leaves = range(8, 14), exclude_stat = None,
                                    xlims = None, ylims = None, fig_size = (15, 10)):
    df_sim = get_simu_results(variety = variety, nplants = nplants)
    fig, axs = plt.subplots(1, 2, figsize = fig_size)
    compare_sim_scan(variety = variety, variable = variable, xaxis = xaxis, 
                     leaves = leaves, exclude_stat = exclude_stat, 
                     xlims = xlims, ylims = ylims, legend = False,
                     ax = axs[0], df_sim = df_sim)
    compare_sim_symptom_tagged(variety = variety, variable = variable, xaxis = xaxis, 
                           leaves = leaves, xlims = xlims, ylims = ylims,
                           ax = axs[1], df_sim = df_sim)
    letters = iter(['a', 'b', 'c', 'd'])
    for ax in axs:
        ax.annotate(next(letters), xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)

def compare_splant(variety = 'Tremie12', nplants=30, variable = 'length', xaxis = 'degree_days',
                   axes = 'MS', group = 'sum', exclude_stat = None, xlims = None, ylims = None):
    """ axes = 'MS', 'tillers' or 'total'"""
    df_sim = get_simu_results(variety = variety, nplants = nplants, group = group)
    df_obs = get_scan_dimensions(variety = variety)
    if axes == 'MS':
        (df_sim, df_obs) = map(lambda x: x[x['axis'] == 'MS'], (df_sim, df_obs))
    elif axes == 'tillers':
        (df_sim, df_obs) = map(lambda x: x[x['axis'].isin(
                                [ax for ax in set(x['axis']) if ax.startswith('T')])], 
                                (df_sim, df_obs))
    if exclude_stat is not None:
        df_obs = df_obs[df_obs['stat'] < exclude_stat]
    
    if group == 'mean':
        ax = plot_mean(df_sim,  variable = variable, xaxis = xaxis,
                       error_bars = False, marker = '', linestyle = '-', color = 'r', 
                       title = 'Simu vs. scans', xlabel = xaxis, ylabel = variable,
                       xlims = xlims, ylims = ylims, return_ax = True)
        plot_mean(df_obs,  variable = variable, xaxis = xaxis, 
                   error_bars = True, marker = 'o', linestyle = '--', color = 'b', 
                   title = 'Simu vs. scans', xlabel = xaxis, ylabel = variable,
                   xlims = xlims, ylims = ylims, ax = ax)
    elif group == 'sum':
        ax = plot_sum(df_sim,  variable = variable, xaxis = xaxis,
                       marker = '', linestyle = '-', color = 'r', 
                       title = 'Simu vs. scans', xlabel = xaxis, ylabel = variable,
                       xlims = xlims, ylims = ylims, return_ax = True)
        plot_sum(df_obs,  variable = variable, xaxis = xaxis, 
                   marker = 'o', linestyle = '--', color = 'b', 
                   title = 'Simu vs. scans', xlabel = xaxis, ylabel = variable,
                   xlims = xlims, ylims = ylims, ax = ax)