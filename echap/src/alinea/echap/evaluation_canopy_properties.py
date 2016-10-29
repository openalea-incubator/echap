""" Reconstruct a MTG wheat canopy at different ages and save LAI, cover fraction,
    intercepted radiation and compare it with field data """

import pandas
import numpy
import os
import glob
import matplotlib.pyplot as plt
from math import sqrt
try:
    import cPickle as pickle
except ImportError:
    import pickle

from alinea.adel.astk_interface import AdelWheat
# Import for wheat reconstruction
import alinea.echap
from openalea.deploy.shared_data import shared_data
from alinea.echap.hs_tt import tt_hs_tag
#from alinea.echap.weather_data import *
import alinea.echap.architectural_data as archidb
from alinea.echap.canopy_data import lai_pai_scan, pai57
from alinea.echap.architectural_reconstructions import (echap_reconstructions, 
                                                        EchapReconstructions,
                                                        HS_fit, density_fits,
                                                        pdict, reconstruction_parameters)
from alinea.adel.postprocessing import ground_cover
from alinea.caribu.caribu_star import run_caribu
from alinea.caribu.light import diffuse_source
from alinea.caribu.label import Label


def cache_reconstruction_path(tag):
    path = shared_data(alinea.echap) / 'cache' / 'reconstructions' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    return path


def cache_simulation_path(tag, rep):
    path = shared_data(alinea.echap) / 'cache' / 'simulations' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    path = shared_data(alinea.echap) / 'cache' / 'simulations' / tag / 'rep_' + str(rep)
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    return path

# Run and save canopy properties ###################################################################
def aggregate_lai(g, axstat):
    colnames = ['aire du plot', 'Nbr.plant.perplot', 'ThermalTime', 'LAI_tot', 'LAI_vert',
                 'PAI_tot', 'PAI_vert']
    pstat = AdelWheat.plot_statistics(g,axstat)
    if pstat is None:
        return pandas.DataFrame([np.nan for i in colnames], columns = colnames)
    else:
        return pstat.loc[:, colnames]

def get_lai_properties(g):
    df_axstat = AdelWheat.axis_statistics(g)
    df_lai_tot = aggregate_lai(g, df_axstat)
    df_lai_MS = aggregate_lai(g, df_axstat[df_axstat['axe_id']=='MS'])
    df_lai_MS.rename(columns = {col:col+'_MS' for col in ['LAI_tot', 'LAI_vert', 
                                                         'PAI_tot', 'PAI_vert']}, inplace=True)
    df_lai_ferti = aggregate_lai(g, df_axstat[df_axstat['has_ear']])
    df_lai_ferti.rename(columns = {col:col+'_ferti' for col in ['LAI_tot', 'LAI_vert', 
                                                            'PAI_tot', 'PAI_vert']}, inplace=True)
    df_lai = pandas.merge(df_lai_tot, df_lai_MS)
    df_lai = pandas.merge(df_lai, df_lai_ferti)
    return df_lai
    
def draft_TC(g, adel, zenith, rep, scale = 1):
    echap_top_camera =  {'type':'perspective', 'distance':200., 
                         'fov':50., 'azimuth':0, 'zenith':zenith}
    gc, im, box = ground_cover(g, adel.domain, camera=echap_top_camera, image_width = int(2144 * scale), 
                                image_height = int(1424*scale), getImages=True, replicate=rep)
    return gc
        
def get_cover_fraction_properties(g, adel, nplants, scale = 1):
    df_TC = pandas.DataFrame(columns = ['TCgreen', 'TCsen', 'TCtot',
                                         'TCgreen_57', 'TCsen_57', 'TCtot_57'])
    if nplants <= 30:
        reps = [1,2]
    else:
        reps = [1,1]
        
    for zenith, rep in zip([0,57],reps):
        dict_TC = draft_TC(g, adel, zenith, rep, scale)
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

def draft_light(g, adel, z_level):
    scene = adel.scene(g)
    sources = diffuse_source(46)
    out = run_caribu(sources, scene, domain=adel.domain, zsoil = z_level)
    labs = {k:Label(v) for k,v in out['label'].iteritems() if not numpy.isnan(float(v))}
    ei_soil = [out['Ei'][k] for k in labs if labs[k].is_soil()]
    return float(numpy.mean(ei_soil))
    
def get_radiation_properties(g, adel, z_levels = [0, 5, 20, 25]):
    df_radiation = pandas.DataFrame(columns = ['LightPenetration_%d' %lev for lev in z_levels])
    for z_level in z_levels:
        df_radiation.loc[0, 'LightPenetration_%d' %z_level] = draft_light(g, adel, z_level)
        df_radiation.loc[0, 'LightInterception_%d' %z_level] = 1 - df_radiation.loc[0, 
                                                                    'LightPenetration_%d' %z_level]
    return df_radiation


def get_reconstruction(variety='Tremie12', nplants=30, tag='reference', rep=1,
                       reset=False, reset_echap=False):
    """ Return a new adel wheat instance"""
    sim_path = cache_simulation_path(tag, rep)

    filename = sim_path / 'adel_' + variety.lower() + '_' + str(nplants) + 'pl.pckl'
    if not reset:
        try:
            with open(filename) as saved:
                return pickle.load(saved)
        except IOError:
            pass
    echap = echap_reconstructions(tag, reset = reset_echap)
    adel = echap.get_reconstruction(name=variety, nplants=nplants)
    # TO DO : make (part of ?) adel picklable !
    # with open(str(filename), 'w') as output:
    #     pickle.dump(adel, output)
    return adel


def as_daydate(daydate, tths):
    if daydate.startswith('T1'):
        return tths.set_index('tag_T1')['daydate'][daydate]
    elif daydate.startswith('T2'):
        return tths.set_index('tag_T2')['daydate'][daydate]
    else:
        return daydate


def daydate_range(variety, tag, start, stop, by=None, at=None):
    tths = tt_hs_tag(variety, tag)
    if at is None:
        if start is None:
            start = tths['daydate'][0]
        else:
            start = as_daydate(start, tths)
        if stop is None:
            stop = tths['daydate'][-1]
        else:
            stop = as_daydate(stop, tths)
        at = tths.set_index('daydate').ix[start:stop:by,].index.values.tolist()
    else:
        at = map(lambda x: as_daydate(x, tths), at)

    return at


def build_canopies(variety='Tremie12', nplants=30, tag='reference', rep=1,
                   start=None, stop=None, by=None, at=None, reset=False):

    sim_path = cache_simulation_path(tag, rep)
    head_path = sim_path / 'canopy' / variety.lower() + '_' + str(
            nplants) + 'pl_'
    if not os.path.exists(sim_path / 'canopy'):
        os.makedirs(sim_path / 'canopy')

    dd_range = daydate_range(variety, tag, start, stop, by, at)
    missing = dd_range
    if not reset:
        pattern = head_path + '*.pckl'
        done = glob.glob(pattern)
        done = map(lambda x: x.split('pl_')[1].split('.')[0], done)
        missing = [d for d in dd_range if d not in done]

    if len(missing) > 0:
        adel = get_reconstruction(variety=variety, nplants=nplants, tag=tag,
                                  rep=rep)
        tths = tt_hs_tag(variety, tag)
        for d in dd_range:
            print d
            basename = head_path + '_'.join(d.split('-'))
            age = tths.set_index('daydate')['TT'][d]
            g = adel.setup_canopy(age=age)
            adel.save(g, basename=str(basename))

    return missing


def get_canopy(variety='Tremie12', nplants=30, daydate='T1', tag='reference', rep=1, reset=False, load_geom=True):
    tths = tt_hs_tag(variety, tag)
    daydate = as_daydate(daydate, tths)

    sim_path = cache_simulation_path(tag, rep)

    if not os.path.exists(sim_path / 'canopy'):
        os.makedirs(sim_path / 'canopy')
    basename = sim_path / 'canopy' / variety.lower() + '_' + str(
        nplants) + 'pl_' + '_'.join(daydate.split('-'))
    if not reset:
        try:
            g = AdelWheat.load(basename=str(basename), load_geom=load_geom)
            return g
        except IOError:
            pass
    adel = get_reconstruction(variety=variety, nplants=nplants, tag=tag, rep=rep)
    age = tths.set_index('daydate')['TT'][daydate]
    g = adel.setup_canopy(age=age)
    adel.save(g, basename=str(basename))
    return g


def simulated_lai(variety='Tremie12', nplants=30, tag='reference', rep=1,
                  start=None, stop=None, by=None, at=None, reset=False):

    sim_path = cache_simulation_path(tag, rep)
    filename = sim_path / 'lai_' + variety.lower() + '_' + str(
        nplants) + 'pl.csv'
    dd_range = daydate_range(variety, tag, start, stop, by, at)
    done = None
    missing = dd_range
    if not reset:
        try:
            df = pandas.read_csv(filename)
            if all(map(lambda x: x in df['daydate'].values, dd_range)):
                return df.loc[df['daydate'].isin(dd_range),:]
            else:
                done = df.loc[df['daydate'].isin(dd_range),:]
                missing = [d for d in dd_range if d not in done['daydate'].values]
        except IOError:
            pass
    print('building missing canopies..')
    build_canopies(variety=variety, nplants=nplants, tag=tag, rep=rep, at=missing, reset=False)
    print('Compute missing plot statistics...')
    new = []
    for d in missing:
        print d
        g = get_canopy(daydate=d, variety=variety,
                                  nplants=nplants, tag=tag, rep=rep, load_geom=False)
        df_lai = get_lai_properties(g)
        df_lai['variety'] = variety
        df_lai['daydate'] = d
        new.append(df_lai)

    df = pandas.concat(new)
    tths = tt_hs_tag(variety, tag)
    df = df.merge(tths)
    df = pandas.concat((done, df))
    df.to_csv(filename, index=False)
    return df



def run_one_simulation(variety = 'Tremie12', nplants = 30, variability_type = None,
                        age_range = [400., 2600.], time_steps = [20, 100],
                        scale_povray = 1., z_levels = [0, 5, 20, 25], 
                        reset = False, reset_data = False, only_lai = False):
    # Temp
    if variety in ['Tremie12', 'Tremie13']:
        pars = reconstruction_parameters()
        # pars['density_tuning'] = pdict(None)
        # pars['density_tuning']['Tremie12'] = 0.85
        # pars['density_tuning']['Tremie13'] = 0.85
        reconst = EchapReconstructions(reset_data=reset_data, pars=pars)
    else:
        reconst = echap_reconstructions(reset=reset, reset_data=reset_data)
    HSconv = reconst.HS_fit[variety]
    adel = reconst.get_reconstruction(name=variety, nplants = nplants)
    ages_1 = numpy.arange(age_range[0], age_range[1], time_steps[0])
    ages_2 = numpy.arange(age_range[0], age_range[1], time_steps[1])
    for age in numpy.unique(ages_1.tolist()+ages_2.tolist()):
        # Get canopy properties
        g = adel.setup_canopy(age=age)
        if age in ages_1:
            df_lai = get_lai_properties(g, adel)
        if age in ages_2:
            if only_lai==False:
                df_cover = get_cover_fraction_properties(g, adel, nplants, scale_povray)
                df_radiation = get_radiation_properties(g, adel, z_levels)
            else:
                cols = ['TCgreen', 'TCsen', 'TCtot', 'TCgreen_57', 'TCsen_57', 'TCtot_57',]
                df_cover = pandas.DataFrame([[numpy.nan for col in cols]], columns = cols)
                cols = ['LightPenetration_%d' %lev for lev in z_levels]
                df_radiation = pandas.DataFrame([[numpy.nan for col in cols]], columns = cols)
                
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
                    nrep = 5, scale_povray = 1., reset = False, reset_data = False,
                    only_lai=False):
    df_props = []
    for rep in range(nrep):
        df_prop = run_one_simulation(variety = variety, nplants = nplants, 
                                       variability_type = variability_type, 
                                       age_range = age_range, time_steps = time_steps,
                                       scale_povray = scale_povray, reset = reset, 
                                       reset_data = reset_data, only_lai = only_lai)
        df_prop.loc[:,'rep'] = rep
        df_props.append(df_prop)
    df_props = pandas.concat(df_props)
    return df_props.reset_index(drop = True)

def get_file_name(variety = 'Tremie12', nplants = 30, nrep = 5, only_lai=False):
    if only_lai==False:
        return variety.lower() + '_' + str(nplants) + 'pl_' + str(nrep) + 'rep.csv'
    else:
        return variety.lower() + '_' + str(nplants) + 'pl_' + str(nrep) + 'rep_lai.csv'
    
def run_and_save_one_simu(variety = 'Tremie12', nplants = 30, variability_type = None,
                          age_range = [400, 2600], time_steps = [20, 100],
                          scale_povray = 1., reset = False, reset_data = False, only_lai=False):                         
    df_save = run_one_simulation(variety = variety, nplants = nplants, 
                                 variability_type = variability_type, age_range = age_range,
                                 time_steps = time_steps, scale_povray = scale_povray,
                                 reset = reset, reset_data = reset_data, only_lai=only_lai)
    filename = get_file_name(variety = variety, nplants = nplants, nrep = 1, only_lai=only_lai)
    file_path = str(shared_data(alinea.echap)/'architectural_simulations'/filename)
    df_save.to_csv(file_path)
    
def run_and_save_multi_simu(variety = 'Tremie12', nplants = 30, variability_type = None,
                            age_range = [400, 2600], time_steps = [20, 100],
                            nrep = 5, scale_povray = 1., reset = False, reset_data = False,
                            only_lai=False):
    df_save = run_multi_simu(variety = variety, nplants = nplants, 
                             variability_type = variability_type, age_range = age_range,
                             time_steps = time_steps, nrep = nrep,
                             scale_povray = scale_povray, reset = reset, reset_data = reset_data,
                             only_lai=only_lai)
    filename = get_file_name(variety = variety, nplants = nplants, nrep = nrep, only_lai=only_lai)
    file_path = str(shared_data(alinea.echap)/'architectural_simulations'/filename)
    df_save.to_csv(file_path)
    
def run_and_save_one_simu_all_varieties(nplants = 30,  variability_type = None,
                                        age_range = [400, 2600], time_steps = [20, 100],
                                        scale_povray = 1., reset = False, 
                                        reset_data = False, only_lai=False):
    varieties = ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']
    for variety in varieties:
        run_and_save_one_simu(variety=variety, nplants=nplants, variability_type=variability_type,
                              age_range=age_range, time_steps=time_steps,
                              scale_povray=scale_povray, reset=reset,
                              reset_data=reset_data, only_lai=only_lai)
                              
def run_and_save_multi_simu_all_varieties(nplants = 30,  variability_type = None,
                                        age_range = [400, 2600], time_steps = [20, 100],
                                        nrep = 5, scale_povray = 1., reset = False, 
                                        reset_data = False, only_lai=False):
    varieties = ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']
    for variety in varieties:
        run_and_save_multi_simu(variety=variety, nplants=nplants, variability_type=variability_type,
                              age_range=age_range, time_steps=time_steps, nrep=nrep,
                              scale_povray=scale_povray, reset=reset,
                              reset_data=reset_data, only_lai=only_lai)
    
# Get and plot canopy properties ###################################################################
def get_simu_results(variety = 'Tremie12', nplants = 30, nrep = 1, only_lai=False):
    filename = get_file_name(variety = variety, nplants = nplants, nrep = nrep, only_lai=only_lai)
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

    if len(lst) <= 1:
        return 0
    n, v = len(lst), variance(lst)
    c = t.interval(perc_conf * 1.0 / 100, n-1)[1]
    
    return sqrt(v/n) * c

def plot_mean(data, variable = 'LAI_vert', xaxis = 'ThermalTime',
              marker = 'd', empty_marker = False, linestyle = '-', color = 'b', 
              title = None, xlabel = None, ylabel = None,
              xlims = None, ylims = None, ax = None):

    var_mean = '_'.join((variable,'mean'))
    var_err = '_'.join((variable,'conf_int'))
    df = data.copy(deep=True)
    if var_mean not in df.columns:
        if variable not in df.columns:
            raise KeyError('variable ' + variable + ' not found in data!')
        else:
            df = df.groupby(xaxis).agg([numpy.mean, conf_int])
            df.columns = ['_'.join(c) for c in df.columns]
            df = df.reset_index()

    if ax is None:
        fig, ax = plt.subplots()
    if empty_marker:
        markerfacecolor = 'none'
    else:
        markerfacecolor = color

    df = df.loc[:, [xaxis, var_mean, var_err]].dropna()
    ax.errorbar(df[xaxis], df[var_mean].values, yerr = df[var_err].values,
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

    return ax
            
def plot_sum(data, variable = 'LAI_vert', xaxis = 'ThermalTime', 
              marker = 'd', linestyle = '-', color = 'b', 
              title = None, xlabel = None, ylabel = None,
              xlims = None, ylims = None, ax = None, return_ax = False):
    
    
    if variable in data.columns:
        if ax == None:
            fig, ax = plt.subplots()
        df = data[pandas.notnull(data.loc[:,variable])]

        df_sum = df.groupby([xaxis, 'num_plant']).sum()
        df_sum = df_sum.reset_index()
        df_sum = df_sum.groupby(xaxis).mean()
        
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
        HSconv = HS_fit()[variety]
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
    HSconv = HS_fit()[variety]
    df_obs['HS'] = HSconv(df_obs['ThermalTime'])
    return df_obs
    
def get_radiation_obs(variety = 'Tremie12'):
    weather = read_weather_variety(variety)
    dfa, tab0, tab20 = mat_ray_obs(variety)
    tab0 = tab0.rename(columns = {'%':'LightPenetration_0'})
    tab20 = tab20.rename(columns = {'%':'LightPenetration_20'})
    df_obs = pandas.merge(tab0, tab20)
    d = pandas.to_datetime((df_obs['datetime'].map(str) + ' ' + '12:00:00').values, dayfirst=True, utc=True)
    df_obs['date'] = df_obs['datetime'].apply(lambda x: datetime.datetime.strptime(x, '%d-%m-%Y'))
    df_obs['ThermalTime'] = weather.data.degree_days[d].values
    HSconv = HS_fit()[variety]
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
    
def plot_sim_obs(df_sim, df_obs=None, variable = 'LAI_vert', xaxis = 'HS',
                 title = None, xlabel = 'Haun Stage', ylabel = 'LAI',
                 colors = ['b', 'r'], markers = ['d', 'o'], linestyles = ['-', '--'],
                 error_bars = [True, True], xlims = None, ylims = None, legend = True, ax = None):
    if ax == None:
        fig, ax = plt.subplots()
    plot_mean(df_sim, variable = variable, xaxis = xaxis,
                color = colors[0], marker = markers[0], linestyle = linestyles[0],
                title = title, xlabel = xlabel, ylabel = ylabel, 
                xlims = xlims, ylims = ylims, ax = ax)
    if df_obs is not None:
        plot_mean(df_obs, variable = variable, xaxis = xaxis,
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

def compare_lai_by_axis(nplants = 30, nrep = 1,
                        fig_size = (10, 10), color = 'b', title = None, 
                        tight_layout = False, only_lai=False):
    fig, axs = plt.subplots(2,2, figsize = fig_size)
    linestyles = {'total':'-', 'MS':'--', 'ferti':'-.'}
    varieties = ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']
    colors = {'Mercia':'r', 'Rht3':'g', 'Tremie12':'b', 'Tremie13':'m'}
    axs_iter = iter(axs.flat)
    proxys = []
    labels = []
    for variety in varieties:
        df_sim = get_simu_results(variety = variety, nplants = nplants, nrep = nrep, only_lai=only_lai)
        df_obs = get_lai_obs(variety = variety, origin = 'biomass')
        ax = next(axs_iter)
        color = colors[variety]
        proxys += [plt.Line2D((0,1),(0,0), color=color, marker='', linestyle='-')]
        labels += [variety]
        plot_sim_obs(df_sim, df_obs, variable = 'LAI_vert', xaxis = 'HS',
                     xlabel = 'Haun stage', ylabel = 'Green Leaf Area Index',
                     colors = [color, color], markers = ['', 'o'], 
                     linestyles = [linestyles['total'], ''],
                     error_bars = [False, True], legend = False, ylims=[0,10], ax = ax)
        plot_sim_obs(df_sim, df_obs, variable = 'LAI_vert_MS', xaxis = 'HS',
                     xlabel = 'Haun stage', ylabel = 'Green Leaf Area Index',
                     colors = [color, color], markers = ['', 'o'], 
                     linestyles = [linestyles['MS'], ''],
                     error_bars = [False, True], legend = False, ylims=[0,10], ax = ax)
        plot_sim_obs(df_sim, df_obs, variable = 'LAI_vert_ferti', xaxis = 'HS',
                     xlabel = 'Haun stage', ylabel = 'Green Leaf Area Index',
                     colors = [color, color], markers = ['', 'o'], 
                     linestyles = [linestyles['ferti'], ''],
                     error_bars = [False, True], legend = False, ylims=[0,10], ax = ax)
    proxys += [plt.Line2D((0,1),(0,0), color='k', marker='', linestyle='-'),
               plt.Line2D((0,1),(0,0), color='k', marker='', linestyle='--'),
               plt.Line2D((0,1),(0,0), color='k', marker='', linestyle='-.')]
    labels += ['Total', 'MS', 'Ferti']
    axs[0][1].legend(proxys, labels, loc='center left', bbox_to_anchor=(1, 0.5))
