""" Reconstruct a MTG wheat canopy at different ages and save LAI, cover fraction,
    intercepted radiation and compare it with field data """

import pandas
import numpy
import matplotlib.pyplot as plt
from math import sqrt

# Import for wheat reconstruction
from openalea.deploy.shared_data import shared_data
from alinea.echap.weather_data import *
import alinea.echap.architectural_data as archidb
from alinea.echap.architectural_reconstructions import (echap_reconstructions, 
                                                        fit_HS, density_fits)
from alinea.adel.postprocessing import axis_statistics, plot_statistics, ground_cover
from alinea.caribu.caribu_star import diffuse_source, run_caribu
from alinea.adel.postprocessing import ground_cover
from alinea.caribu.label import Label

# Run and save canopy properties ###################################################################
def get_lai_properties(g, adel, nplants):
    df_lai = adel.get_exposed_areas(g, convert=True)
    colnames = ['aire du plot', 'Nbr.plant.perplot', 'ThermalTime', 'LAI_tot', 'LAI_vert',
                    'PAI_tot', 'PAI_vert']
    if not df_lai.empty:
        df_axstat, _ = axis_statistics(df_lai, adel.domain_area, adel.convUnit)
        df_axstat = plot_statistics(df_axstat, nplants, adel.domain_area)
        return df_axstat.loc[:, colnames]
    else: 
        colnames = ['aire du plot', 'Nbr.plant.perplot', 'ThermalTime', 'LAI_tot', 'LAI_vert',
                    'PAI_tot', 'PAI_vert']
        return pandas.DataFrame([np.nan for i in colnames], columns = colnames)

def draft_TC(g, adel, domain, zenith, rep, scale = 1):
    echap_top_camera =  {'type':'perspective', 'distance':200., 
                         'fov':50., 'azimuth':0, 'zenith':zenith}
    gc, im, box = ground_cover(g, domain, camera=echap_top_camera, image_width = int(2144 * scale), 
                                image_height = int(1424*scale), getImages=True, replicate=rep)
    return gc
        
def get_cover_fraction_properties(g, adel, nplants, scale = 1):
    df_TC = pandas.DataFrame(columns = ['TCgreen', 'TCsen', 'TCtot',
                                         'TCgreen_57', 'TCsen_57', 'TCtot_57',])
    if nplants <= 30:
        reps = [1,2]
    else:
        reps = [1,1]
        
    for zenith, rep in zip([0,57],reps):
        dict_TC = draft_TC(g, adel, adel.domain, zenith, rep, scale)
        if zenith == 0:
            suffix = ''
        elif zenith == 57:
            suffix = '_57'
        df_TC.loc[0, 'TCgreen'+suffix] = dict_TC['green']
        df_TC.loc[0, 'TCsen'+suffix] = dict_TC['senescent']
        df_TC.loc[0, 'TCtot'+suffix] = sum(dict_TC.values())
        df_TC.loc[0, 'Gapgreen'+suffix] = 1 - df_TC.loc[0, 'TCgreen'+suffix]
        df_TC.loc[0, 'Gapsen'+suffix] = 1 - df_TC.loc[0, 'TCsen'+suffix]
        df_TC.loc[0, 'Gaptot'+suffix] = 1 - df_TC.loc[0, 'TCtot'+suffix]
    return df_TC

def draft_light(g, adel, domain, z_level):
    scene = adel.scene(g)
    sources = diffuse_source(46)
    out = run_caribu(sources, scene, domain=domain, zsoil = z_level)
    labs = {k:Label(v) for k,v in out['label'].iteritems() if not numpy.isnan(float(v))}
    ei_soil = [out['Ei'][k] for k in labs if labs[k].is_soil()]
    return float(numpy.mean(ei_soil))
    
def get_radiation_properties(g, adel, z_levels = [0, 5, 20, 25]):
    df_radiation = pandas.DataFrame(columns = ['LightPenetration_%d' %lev for lev in z_levels])
    for z_level in z_levels:
        df_radiation.loc[0, 'LightPenetration_%d' %z_level] = draft_light(g, adel, 
                                                                          adel.domain, z_level)
        df_radiation.loc[0, 'LightInterception_%d' %z_level] = 1 - df_radiation.loc[0, 
                                                                    'LightPenetration_%d' %z_level]
    return df_radiation
    
def run_one_simulation(variety = 'Tremie12', nplants = 30, variability_type = None,
                        age_range = [400., 2600.], time_steps = [20, 100], 
                        scale_povray = 1., reset = False, reset_data = False):
    reconst = echap_reconstructions(reset=reset, reset_data=reset_data)
    HSconv = reconst.HS_fit[variety]
    adel = reconst.get_reconstruction(name=variety, nplants = nplants)
    ages_1 = numpy.arange(age_range[0], age_range[1], time_steps[0])
    ages_2 = numpy.arange(age_range[0], age_range[1], time_steps[1])
    for age in numpy.unique(ages_1.tolist()+ages_2.tolist()):
        # Get canopy properties
        g = adel.setup_canopy(age=age)
        if age in ages_1:
            df_lai = get_lai_properties(g, adel, nplants)
        if age in ages_2:
            df_cover = get_cover_fraction_properties(g, adel, nplants, scale_povray)
            df_radiation = get_radiation_properties(g, adel)
        
        # Group in same dataframe
        if age == age_range[0]:
            df_prop = pandas.DataFrame(columns = ['variety', 'HS']+
                                                 [col for col in df_lai.columns]+
                                                 [col for col in df_cover.columns]+
                                                 [col for col in df_radiation.columns])
        if age in ages_1:
            df_prop.loc[age, df_lai.columns] = df_lai.loc[0, :]
        if age in ages_2:
            df_prop.loc[age, df_cover.columns] = df_cover.loc[0, :]
            df_prop.loc[age, df_radiation.columns] = df_radiation.loc[0, :]
        df_prop.loc[age, 'variety'] = variety
        df_prop.loc[age, 'HS'] = HSconv(age)
    df_prop.reset_index(drop = True, inplace = True)
    df_prop = df_prop.applymap(lambda x: numpy.float(x) if numpy.isreal(x) else x)
    return df_prop

def run_multi_simu(variety = 'Tremie12', nplants = 30, 
                    variability_type = None,
                    age_range = [400, 2600], time_steps = [20, 100],
                    filename = 'canopy_properties_tremie12_5rep.csv',
                    nrep = 5, scale_povray = 1., reset_reconst = False):
    df_props = []
    for rep in range(nrep):
        df_prop = run_one_simulation(variety = variety, nplants = nplants, 
                                       variability_type = variability_type, 
                                       age_range = age_range, time_steps = time_steps,
                                       scale_povray = scale_povray, reset_reconst = reset_reconst)
        df_prop.loc[:,'rep'] = rep
        df_props.append(df_prop)
    df_props = pandas.concat(df_props)
    return df_props.reset_index(drop = True)

def get_file_name(variety = 'Tremie12', nplants = 30, nrep = 5):
    return variety.lower() + '_' + str(nplants) + 'pl_' + str(nrep) + 'rep.csv'
    
def run_and_save_one_simu(variety = 'Tremie12', nplants = 30, variability_type = None,
                          age_range = [400, 2600], time_steps = [20, 100],
                          scale_povray = 1., reset_reconst = False):                         
    df_save = run_one_simulation(variety = variety, nplants = nplants, 
                                 variability_type = variability_type, age_range = age_range,
                                 time_steps = time_steps, scale_povray = scale_povray,
                                 reset_reconst = reset_reconst)
    filename = get_file_name(variety = variety, nplants = nplants, nrep = 1)
    file_path = str(shared_data(alinea.echap)/'architectural_simulations'/filename)
    df_save.to_csv(file_path)
    
def run_and_save_multi_simu(variety = 'Tremie12', nplants = 30, 
                            variability_type = None,
                            age_range = [400, 2600], time_steps = [20, 100],
                            nrep = 5, scale_povray = 1., reset_reconst = False):
    df_save = run_multi_simu(variety = variety, nplants = nplants, 
                             variability_type = variability_type, age_range = age_range,
                             time_steps = time_steps, nrep = nrep,
                             scale_povray = scale_povray, reset_reconst = reset_reconst)
    filename = get_file_name(variety = variety, nplants = nplants, nrep = nrep)
    file_path = str(shared_data(alinea.echap)/'architectural_simulations'/filename)
    df_save.to_csv(file_path)
    
# Get and plot canopy properties ###################################################################
def get_simu_results(variety = 'Tremie12', nplants = 30, nrep = 1):
    filename = get_file_name(variety = variety, nplants = nplants, nrep = nrep)
    file_path = str(shared_data(alinea.echap)/'architectural_simulations'/filename)
    return pandas.read_csv(file_path)

def variance(lst):
    """
    Uses standard variance formula (sum of each (data point - mean) squared)
    all divided by number of data points
    """
    mu = numpy.mean(lst)
    return 1.0/(len(lst)-1) * sum([(i-mu)**2 for i in lst])
        
def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor can be looked up from a table, for 95pcent confidence
    on a reasonable size sample (>=500) 1.96 is used.
    """
    from scipy.stats import t
    
    n, v = len(lst), variance(lst)
    c = t.interval(perc_conf * 1.0 / 100, n-1)[1]
    
    return sqrt(v/n) * c

def plot_mean(data, variable = 'LAI_vert', xaxis = 'ThermalTime', 
              error_bars = False, error_method = 'confidence_interval', 
              marker = 'd', empty_marker = False, linestyle = '-', color = 'b', 
              title = None, xlabel = None, ylabel = None,
              xlims = None, ylims = None, ax = None, return_ax = False):
    if variable in data.columns:
        if ax == None:
            fig, ax = plt.subplots()
        if empty_marker == False:
            markerfacecolor = color
        else:
            markerfacecolor = 'none'
            
        df = data[pandas.notnull(data.loc[:,variable])].loc[:, [xaxis, variable]]
        df_mean = df.groupby(xaxis).mean()
        df['nb_rep'] = map(lambda x: df[xaxis].value_counts()[x], df[xaxis])
        if error_bars == True and len(df['nb_rep'])>0 and min(df['nb_rep'])>1:
            if error_method == 'confidence_interval':
                df_err = df.groupby(xaxis).agg(conf_int)
            elif error_method == 'std_deviation':
                df_err = df.groupby(xaxis).std()
            else:
                raise ValueError("'error_method' unknown: 'try confidence_interval' or 'std_deviation'")
            ax.errorbar(df_mean.index, df_mean[variable], yerr = df_err[variable].values,
                        marker = marker, linestyle = linestyle, color = color,
                        markerfacecolor = markerfacecolor,  markeredgecolor = color)
        else:
            ax.plot(df_mean.index, df_mean[variable],
                    marker = marker, linestyle = linestyle, color = color,
                    markerfacecolor = markerfacecolor,  markeredgecolor = color)
        if title is not None:
            ax.set_title(title, fontsize = 18)
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize = 18)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize = 18)
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        if return_ax == True:
            return ax
            
def plot_sum(data, variable = 'LAI_vert', xaxis = 'ThermalTime', 
              marker = 'd', linestyle = '-', color = 'b', 
              title = None, xlabel = None, ylabel = None,
              xlims = None, ylims = None, ax = None, return_ax = False):
    
    def get_mean_data(variable):
        if variable in data.columns:
            df_sum.loc[:, variable] = df_mean.loc[:, variable]
    
    if variable in data.columns:
        if ax == None:
            fig, ax = plt.subplots()
        df = data[pandas.notnull(data.loc[:,variable])]
        df_mean = df.groupby(xaxis).mean()

        df_sum = df.groupby(xaxis).sum()
        for var in ['HS', 'num_leaf_bottom', 'num_leaf_top', 'fnl', 
                    'cur_max_leaf_top', 'cur_num_leaf_top']:
            if var != 'xaxis':
                get_mean_data(var)

        ax.plot(df_sum.index, df_sum[variable], color = color, 
                    marker = marker, linestyle = linestyle)
        if title is not None:
            ax.set_title(title, fontsize = 18)
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize = 18)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize = 18)
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        if return_ax == True:
            return ax

# Get observation data for canopy properties #######################################################        
def get_lai_obs(variety = 'Tremie12', origin = 'biomass'):
    
    def get_lai_from_single_origin(origin_):
        if origin_ == 'biomass':
            df_obs = archidb.LAI_biomasse_data()[variety]
            variable = 'LAI_vert_biomasse'
        elif origin_ == 'photo':
            df_obs = archidb.PAI_photo_data()[variety]
            variable = 'PAI_vert_photo'

        df_obs = df_obs.rename(columns = {variable:'LAI_vert', 'TT_date':'ThermalTime'})
        HSconv = fit_HS()[variety]
        df_obs['HS'] = HSconv(df_obs['ThermalTime'])
        df_obs['origin'] = origin_
        if 'date' in df_obs.columns:
            df_obs['date'] = df_obs['date'].apply(lambda x: datetime.datetime.strptime(x, '%Y-%m-%d'))
        return df_obs
        
    if variety in ['Mercia', 'Rht3'] and origin == 'biomass':
        # columns = ['LAI_vert', 'ThermalTime', 'HS', 'origin']
        return None
    if origin in ['biomass', 'photo']:
        return get_lai_from_single_origin(origin)
    elif origin == 'all':
        df_lai = pandas.concat([get_lai_from_single_origin('photo'), 
                                get_lai_from_single_origin('biomass')])
        df_lai = df_lai.sort('HS')
        return df_lai.reset_index(drop = True)
    else:
        raise ValueError("Choose origin between biomass, photo")
        
def get_TC_obs(variety = 'Tremie12'):
    df_obs_0 = archidb.TC_data()[variety+'_0']
    df_obs_0 = df_obs_0.rename(columns = {'TC':'TCgreen'})
    df_obs_57 = archidb.TC_data()[variety+'_57']
    df_obs_57 = df_obs_57.rename(columns = {'TC':'TCgreen_57'})
    df_obs = pandas.merge(df_obs_0, df_obs_57)
    df_obs['date'] = df_obs['date'].apply(lambda x: datetime.datetime.strptime(x, '%Y-%m-%d'))
    df_obs['Gapgreen'] = 1 - df_obs['TCgreen']
    df_obs['Gapgreen_57'] = 1 - df_obs['TCgreen_57']
    df_obs = df_obs.rename(columns = {'TT':'ThermalTime'})
    HSconv = fit_HS()[variety]
    df_obs['HS'] = HSconv(df_obs['ThermalTime'])
    return df_obs
    
def get_radiation_obs(variety = 'Tremie12'):
    dfa, tab0, tab20 = mat_ray_obs(variety)
    tab0 = tab0.rename(columns = {'%':'LightPenetration_0'})
    tab20 = tab20.rename(columns = {'%':'LightPenetration_20'})
    df_obs = pandas.merge(tab0, tab20)
    weather = read_weather_variety(variety)
    df_obs = df_obs.rename(columns = {'datetime':'date'})
    df_obs['date'] = df_obs['date'].apply(lambda x: datetime.datetime.strptime(x, '%d-%m-%Y'))
    df_obs['ThermalTime'] = weather.data.degree_days[df_obs['date']].values
    HSconv = fit_HS()[variety]
    df_obs['HS'] = HSconv(df_obs['ThermalTime'])
    df_obs['LightInterception_0'] = 1 - df_obs['LightPenetration_0']
    df_obs['LightInterception_20'] = 1 - df_obs['LightPenetration_20']
    return df_obs

def get_all_obs(variety = 'Tremie12', origin_lai_data = 'biomass'):
    df_lai = get_lai_obs(variety = variety, origin = origin_lai_data)
    df_TC = get_TC_obs(variety = variety)
    df_rad = get_radiation_obs(variety = variety)
    df_all = pandas.concat([df_lai, df_TC, df_rad])
    return df_all.reset_index(drop = True)
    
def plot_sim_obs(df_sim, df_obs, variable = 'LAI_vert', xaxis = 'HS', 
                 title = None, xlabel = 'Haun Stage', ylabel = 'LAI',
                 colors = ['b', 'r'], markers = ['d', 'o'], linestyles = ['-', '--'],
                 error_bars = [True, True], xlims = None, ylims = None, legend = True, ax = None):
    if ax == None:
        fig, ax = plt.subplots()
    plot_mean(df_sim, variable = variable, xaxis = xaxis, error_bars = error_bars[0], 
                color = colors[0], marker = markers[0], linestyle = linestyles[0],
                ax = ax)
        
    plot_mean(df_obs, variable = variable, xaxis = xaxis, error_bars = error_bars[1], 
                color = colors[1], marker = markers[1], linestyle = linestyles[1],
                title = title, xlabel = xlabel, ylabel = ylabel, 
                xlims = xlims, ylims = ylims, ax = ax)
    
    if legend == True:   
        ax.legend(['Simulated', 'Observed'], loc='center left', bbox_to_anchor=(1, 0.5))

def compare_sim_obs(variety = 'Tremie12', nplants = 30, nrep = 1,
                    origin = 'biomass', variable = 'LAI_vert', xaxis = 'HS',
                    title = None, xlabel = 'Haun Stage', ylabel = None,
                    colors = ['b', 'r'], markers = ['d', 'o'], linestyles = ['-', '--'],
                    error_bars = [True, True], xlims = None, ylims = None, ax = None):
    df_sim = get_simu_results(variety = variety)
    df_obs = get_all_obs(variety = variety, origin_lai_data = origin)
    plot_sim_obs(df_sim, df_obs, variable = variable, xaxis = xaxis,
                     title = title, xlabel = xlabel, ylabel = ylabel,
                     colors = colors, markers = markers, linestyles = linestyles,
                     error_bars = error_bars, xlims = xlims, ylims = ylims)
                     
def test_architecture_canopy_single(variety = 'Tremie12', nplants = 30, nrep = 1,
                                     fig_size = (10, 10), color = 'b', title = None, 
                                     axs = None, tight_layout = False):
    if axs == None:
        fig, axs = plt.subplots(2,2, figsize = fig_size)
    df_sim = get_simu_results(variety = variety, nplants = nplants, nrep = nrep)
    df_obs = get_all_obs(variety = variety, origin_lai_data = 'biomass')
    plot_sim_obs(df_sim, df_obs, variable = 'TCgreen_57', xaxis = 'HS',
                     xlabel = 'Haun stage', ylabel = 'Cover fraction (oblique view)',
                     colors = [color, color], markers = ['', 'o'], linestyles = ['-', '--'],
                     error_bars = [False, True], legend = False, ax = axs[0][0])
    plot_sim_obs(df_sim, df_obs, variable = 'TCgreen', xaxis = 'HS',
                     xlabel = 'Haun stage', ylabel = 'Cover fraction (vertical view)',
                     colors = [color, color], markers = ['', 'o'], linestyles = ['-', '--'],
                     error_bars = [False, True], legend = True, ax = axs[0][1])
    plot_sim_obs(df_sim, df_obs, variable = 'LAI_vert', xaxis = 'HS',
                     xlabel = 'Haun stage', ylabel = 'Green Leaf Area Index',
                     colors = [color, color], markers = ['', 'o'], linestyles = ['-', '--'],
                     error_bars = [False, True], legend = False, ax = axs[1][0])
    # Temp
    if not 'LightInterception_0' in df_sim.columns:
        df_sim['LightInterception_0'] = 1 - df_sim['LightPenetration_0']
    plot_sim_obs(df_sim, df_obs, variable = 'LightInterception_0', xaxis = 'HS',
                     xlabel = 'Haun stage', ylabel = 'Intercepted fraction of radiation',
                     colors = [color, color], markers = ['', ''], linestyles = ['-', '--'],
                     error_bars = [False, False], legend = False, ax = axs[1][1])
    letters = iter(['a', 'b', 'c', 'd'])
    for ax in axs.flat:
        ax.annotate(next(letters), xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    if title is not None:
        ax.set_title(title, fontsize = 18)
    if tight_layout == True:
        plt.tight_layout()
    
def test_architecture_canopy(varieties = ['Mercia', 'Rht3', 'Tremie12', 'Tremie13'],
                             nplants = 30, nrep = 1, fig_size = (15, 15), color = 'b', 
                             title = None, tight_layout = False):
    fig, axs = plt.subplots(2,2, figsize = fig_size)
    colors = {'Mercia':'r', 'Rht3':'g', 'Tremie12':'b', 'Tremie13':'m'}
    proxy = []
    for variety in varieties:
        color = colors[variety]
        test_architecture_canopy_single(variety = variety, nplants = nplants, nrep = nrep,
                                        fig_size = fig_size, color = color, title = title, 
                                        axs = axs, tight_layout = tight_layout)
        proxy += [plt.Line2D((0,1),(0,0), color=color, marker='', linestyle='-')]
    axs[1][1].legend(proxy, varieties, loc='center left', bbox_to_anchor=(1, 0.5))
    
# variety = 'Tremie12'
# nplants = 30
# nrep = 1
# df_sim = get_simu_results(variety = variety, nplants = nplants, nrep = nrep)
# df_dens = density_fits()['Tremie12']

