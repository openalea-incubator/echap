""" Reconstruction of Wheat for Boigneville data
"""

import pandas
import numpy
import math
import json
import os

from copy import deepcopy
try:
    import pickle as pickle
except ImportError:
    import pickle


from alinea.adel.astk_interface import AdelWheat
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.AdelR import devCsv
from openalea.deploy.shared_data import shared_data

import alinea.adel.plantgen_extensions as pgen_ext

import alinea.echap
import alinea.echap.architectural_data as archidb
import alinea.echap.plot_architectural_reconstructions as archi_plot
from alinea.echap.hs_tt import HS_fit



# run_plots = False # prevent ipython %run to make plots
# reset_data = False# control the ipython %run behavior concerning data + HSfit dependent data

share_dir = shared_data(alinea.echap,share_path='./share/data')

#---------- reconstructions

def cache_reconstruction_path(tag):
    path = share_dir / 'cache' / 'reconstructions' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    return path
#
# Fitting parameters
#
def pdict(value):
    """ create a parameter dict for all echap cultivar with value
    """
    return {k:deepcopy(value) for k in ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')}


def reconstruction_parameters(tag='reference', reset=False):
    """ Manual parameters used for fitting Echap experiments
    """

    path = cache_reconstruction_path(tag) / 'reconstruction_parameters.json'

    if not reset:
        try:
            with open(path, 'r') as input_file:
                cached = json.load(input_file)
            return cached
        except IOError:
            pass
    pars={}
    #
    # PlantGen parameters
    #--------------------
    pars['pgen_base'] = {'TT_hs_break':None} # Add {'inner_params':{'DELAIS_PHYLL_SEN_DISP':20}} # 3 par default param de get_reconstruction

    #TO CHECK : check that tiller reduction factor (new plantgen par should or not be set to its default value (now 0.3, but probably 1 at the time of our fitting)
    
    # Nff composition of the canopy
    # if None, it uses the nff composition of the tillering data
    # pars['nff_probabilities']['Tremie12'] = {'11': 0.3, '12': 0.7}
    pars['nff_probabilities'] = pdict(None)
    # TO do : take into account
    #
    # Adel parameters
    #----------------
    pars['adel_pars'] = {'senescence_leaf_shrink' : 0.5,'leafDuration' : 2, 'fracLeaf' : 0.2, 'stemDuration' : 2. / 1.2, 'dHS_col' : 0.2,'dHS_en': 0,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1, 'rate_inclination_tiller': 30, 'drop_empty':True}
    #
    # Green Leaves (ssi) = f(HS)
    #
    # in the original pgen, dHS_bolting is estimated as the decimal number of elongated internodes
    pars['GLpars'] = pdict({'GL_start_senescence':4.4, 'GL_bolting':2.2, 'GL_flag': 4.4, 'dHS_bolting': 3.7, 'curvature' : -0.002})
    pars['GLpars']['Tremie12'].update({'GL_start_senescence':4.8, 'GL_bolting':2.1, 'GL_flag': 5, 'dHS_bolting': 4, 'curvature' : -0.005})
    pars['GLpars']['Tremie13'].update({'GL_bolting':2.1})
    pars['GLpars']['Rht3'].update({'GL_bolting':2.3, 'GL_start_senescence':3.4})
    pars['GLpars']['Mercia'].update({'GL_start_senescence':3.6})
    #
    # Plant Density
    #--------------
    # plant damages. Comment out to fit without damages
    # when indicates the period (Haun stages) during which damage effectively occur (ie axes disappear). 
    # damage indicates the proportion of plant lost during the period
    pars['plant_damages'] = {'Mercia':{'when':[6,13], 'damage': 0.25}, #ammount estimated from plant density estmates 
                             'Rht3':{'when':[7,13], 'damage': 0.3},
                             'Tremie12': None, 'Tremie13': None}
    # density tuning factor(multiplicative effect on density curve)
    pars['density_tuning'] = {'Mercia':None, 'Rht3':None, 'Tremie12':None, 'Tremie13':None}
    # variability of plant emergence
    pars['std_emergence'] = pdict(None)
    # Tillering
    #----------
    # max order of tillering or None
    pars['max_order'] = None
    # delay between stop of growth and disparition of a tiller
    pars['delta_stop_del'] = pdict(2.)
    # in the original pgen, dHS_reg is estimated as the decimal number of elongated internode + 0.36
    pars['dHS_reg'] = {'Mercia':3.7, 'Rht3':3.7, 'Tremie12':4., 'Tremie13':3.7}
    # reduction of emmission
    pars['emission_reduction'] = {'Mercia':None, 'Rht3':None, 'Tremie12':None, 'Tremie13':{'T4':0.5}}
    # damages to tillers
    # when indicates the period (Haun stages) during which damage effectively occur (ie axes disappear).
    # when will be adapted for each tiller so that damaged tillers can emerge before dying    
    # damage indicates the proportion of tiller lost during the period (as a relative proportion of tillers emited)
    # natural regression participate to damage counts, and damaged tillers participate to natural regression (ie damage to a future regressing tiller is not inducing a new regression)
    #pars['tiller_damages'] = None
    pars['tiller_damages'] = {'Mercia':{'when':[6,13],'damage':{t:0.25 for t in ('T1','T2','T3')}}, 
                              'Rht3':{'when':[7,13],'damage':{t:0.3 for t in ('T1','T2','T3')}}, 
                              # 'Tremie12':{'when':[4.9,5.1],'damage':{t:0.2 for t in ['T3']}},
                              'Tremie12': None,
                              'Tremie13': None}
    # multiplicative tuning of the ears/plant data (allows to change MS/tiller ratio)
    pars['ears_per_plant_tuning'] = {'Mercia':None, 'Rht3':None, 'Tremie12':1.5, 'Tremie13':None}
    # Dimensions
    #-----------
    #
    #deformation of the standard profile at cardinal points
    #
    pars['card_scale'] = {'Mercia':{'L_blade':(1, 1.3, 1)},
                          'Rht3':{'L_blade':(1, 0.7, 1)}, 
                          'Tremie12':{'L_blade':(0.7, 1, 1)},
                          'Tremie13':{'L_blade':(0.7, 1, 1), 'L_sheath':(0.7,1,1)}
                          }
    #
    # Leaf geometry
    #--------------
    # flag for forcing median leaf shape = f(age) model instead of random shape sampling
    pars['median_leaf'] = pdict(True)
    # median leaf shape data to be used
    pars['xy_data'] = {'Mercia': 'MerciaRht_byleafclass',
             'Rht3': 'MerciaRht_byleafclass',
             'Tremie12':'Tremie_byleafclass',
             'Tremie13':'Tremie_byleafclass'}
    # number of leaves to be considered as leaf class 2 (top leaves)
    #pars['top_leaves'] = {'Mercia':3, 'Rht3':2, 'Tremie12':4, 'Tremie13':3} # values optimised for echap report december 2014
    pars['top_leaves'] = {'Mercia':4, 'Rht3':4, 'Tremie12':4, 'Tremie13':4}
    #
    pars['disc_level'] = pdict(7)
    #
    with open(path, 'w') as output_file:
        json.dump(pars, output_file, sort_keys=True, indent=4,
                  separators=(',', ': '))
    #
    return pars

#
# --------------------------------------------------------------- Fit HS = f(TT)
#
def check_hs_fits(tag='reference'):
    fits = HS_fit(tag)
    obs = archidb.haun_stage_aggregated()
    archi_plot.haun_stage_plot(obs, fits)

#
# --------------------------------------------------------------- Fit GL = f(HS since flag)
#
def GL_fits(tag='reference', reset=False):
    """ Regression for Green leaf data
    """
    parameters = reconstruction_parameters(tag, reset=reset)
    hsfit = HS_fit(tag)
    GLpars = parameters.get('GLpars')

    return {k:pgen_ext.GreenLeaves(hsfit[k], **GLpars[k]) for k in GLpars}

    
def check_green_leaves(tag='reference', reset=False):
    hsfit = HS_fit(tag)
    gl_obs_nff, gl_obs = archidb.green_leaves_aggregated(hsfit)
    gl_est_nff, gl_est = archidb.green_leaves_estimated(hsfit)
    # Compare varieties
    obs = (gl_obs, gl_est)
    fits = GL_fits(tag, reset=reset)
    archi_plot.green_leaves_plot_mean(obs, fits)
    # Compare nff
    obs = (gl_obs_nff, gl_est_nff, gl_obs, gl_est)
    archi_plot.green_leaves_plot(obs, fits)


#
# --------------------------------------------------------------- Fit Plant density = f(HS)
#
#
# use mean plant density and/or density at the end as fits
class deepdd(pandas.DataFrame):

    def deepcopy(self):
        return self.copy(deep=True)
#


def density_fits(tag='reference', reset=False):
    """
    Manual fit of plant density based on mean plant density and estimate of plant density at harvest
    """

    HS_converter = HS_fit(tag)
    parameters = reconstruction_parameters(tag, reset=reset)

    data = archidb.reconstruction_data()
    pdb = data.Plot_data
    tdb = data.Tillering_data

    damages = parameters.get('plant_damages')

    # compile density data
    density_data = {}
    for k in pdb:
        hs = []
        d = []
        if 'plant_density_at_emergence' in pdb[k]:
            hs.append(HS_converter[k](
                pdb[k]['TT_date'][pdb[k]['code_date']['emergence']]))
            d.append(pdb[k]['plant_density_at_emergence'])
        if 'plant_density' in pdb[k]:
            for date, densities in pdb[k]['plant_density'].items():
                hs.extend(
                    [HS_converter[k](pdb[k]['TT_date'][date])] * len(densities))
                d.extend(densities)
        if 'raw_ear_density_at_harvest' in pdb[k]:
            epp = tdb[k]['ears_per_plant']
            if 'ears_per_plant_tuning' in parameters:
                if parameters['ears_per_plant_tuning'][k] is not None:
                    epp *= parameters['ears_per_plant_tuning'][k]
            ear_d = numpy.array(pdb[k]['raw_ear_density_at_harvest']) / epp
            hs_harvest = HS_converter[k](
                pdb[k]['TT_date'][pdb[k]['code_date']['harvest']])
            hs.extend([hs_harvest] * len(ear_d))
            d.extend(ear_d.tolist())
        density_data[k] = pandas.DataFrame({'HS': hs, 'density': d})

    # estimate mean plant density or mean plant density after damages
    if damages is not None:
        for k in damages:
            if damages[k] is not None:
                density_data[k] = density_data[k].loc[
                                  density_data[k]['HS'] > max(
                                      damages[k]['when']), :]
    mean_plant_density = {k: round(density_data[k]['density'].mean()) for k in
                          density_data}
    #
    fits = {k: {'HS': [-20, 50], 'density': [mean_plant_density[k]] * 2} for k
            in mean_plant_density}
    #
    if damages is not None:
        for k in damages:
            if damages[k] is not None:
                density_at_emergence = round(
                    mean_plant_density[k] * 1. / (1 - damages[k]['damage']))
                density_at_harvest = mean_plant_density[k]
                fits[k] = {'HS': [-20] + damages[k]['when'] + [50],
                           'density': [density_at_emergence] * 2 + [
                                                                       density_at_harvest] * 2}

    fits = {k: {'sowing_density': pdb[k]['sowing_density'],
                'inter_row': pdb[k]['inter_row'],
                'density_table': pandas.DataFrame(fits[k])} for k in fits}

    adjust = parameters.get('density_tuning')
    if adjust is not None:
        for k in adjust:
            if adjust[k] is not None:
                df = fits[k]['density_table']
                df['density'] *= adjust[k]
                fits[k]['sowing_density'] *= adjust[k]
                fits[k]['inter_row'] /= math.sqrt(adjust[k])

    for k in fits:
        df = fits[k]['density_table']
        df['TT'] = HS_converter[k].TT(df['HS'])
        fits[k]['density_at_emergence'] = df['density'].iloc[0]
        fits[k]['density_at_harvest'] = df['density'].iloc[-1]

    return fits


#
def check_density_fits(tag='reference', reset=False):
    fits = density_fits(tag, reset=reset)
    obs = archidb.PlantDensity()
    HSfit = HS_fit(tag)
    archi_plot.density_plot(obs, fits, HSfit)
#
# --------------------------------------------------------------- Fit Tillering
#
# fitting functions
def _axepop_fit(tdb, delta_stop_del, dHS_reg, max_order, tiller_damages,
                std_em):
    ms_nff_probas = tdb['nff_probabilities']
    ears_per_plant = tdb['ears_per_plant']
    primary_emission = tdb['emission_probabilities']
    primary_emission = {k:v for k,v in primary_emission.items() if k not in ['TC','TP','TS','TT','FT','TT3F']}
    emission = pgen_ext.TillerEmission(primary_emission)
    regression = pgen_ext.TillerRegression(ears_per_plant, dHS_reg,
                                           delta_stop_del)
    return pgen_ext.AxePop(ms_nff_probas, emission, regression,
                           tiller_damages=tiller_damages, max_order=max_order,
                           std_em=std_em)


#
# Fits
def axepop_fits(tag='reference', reset=False):
    parameters = reconstruction_parameters(tag, reset=reset)
    # tillering model parameters
    delta_stop_del = parameters.get('delta_stop_del')
    dHS_reg = parameters.get('dHS_reg')
    max_order = parameters.get('max_order')
    tiller_damages = parameters.get('tiller_damages')
    # Emission probabilities tuning
    #
    data = archidb.reconstruction_data()
    tdb = deepcopy(data.Tillering_data)
    # apply manual reductions
    emf = parameters.get('emission_reduction')
    if emf is not None:
        for k in emf:
            if emf[k] is not None:
                for T in emf[k]:
                    tdb[k]['emission_probabilities'][T] *= emf[k][T]
    eppt = parameters.get('ears_per_plant_tuning')
    if eppt is not None:
        for k in eppt:
            if eppt[k] is not None:
                tdb[k]['ears_per_plant'] *= eppt[k]


    # variability of emergence
    std_em = parameters.get('std_emergence')
    for k in std_em:
        if std_em[k] == None:
            hs = HS_fit()
            std_em[k] = hs[k].std_TT_hs_0

    fits={k: _axepop_fit(tdb[k], delta_stop_del=delta_stop_del[k],dHS_reg=dHS_reg[k], max_order=max_order, tiller_damages=tiller_damages[k], std_em = std_em[k]) for k in tdb}
    
    return fits


def check_tillering_fits(tag='reference', reset=False):
    # populate with fits
    fits = axepop_fits(tag, reset=reset)
    hs_fits = HS_fit(tag)
    
    #1 graph tillering par var -> 4 graph sur une feuille
    archi_plot.multi_plot_tillering(archidb.tillers_per_plant(), fits, hs_fits)
    
    #tallage primary pour les 4 var   
    #archi_plot.tillering_primary(archidb.tillers_per_plant(), fits, hs_fits)
    
    #tallage total pour les 4 var    
    #archi_plot.tillering_tot(archidb.tillers_per_plant(), fits, hs_fits)
   
   # primary emissions
    #archi_plot.graph_primary_emission(archidb.emission_probabilities_table())
#
# --------------------------------------------------------------- Fit dimensions
#


def estimate_Wb_Tremie13(dim_fit, data, tag='reference', sep_up_down=4):
    def get_estimate(A, L, ff):
        return float(A)/(L*ff)
    data = data.reset_index()
    data = data[data['label']=='Tremie13']
    hs = HS_fit(tag)
    fits = leafshape_fits(tag)
    W_blades = []
    estimated = []
    for ind in data.index:
        if (numpy.isnan(data.loc[ind, 'W_blade']) and
            not numpy.isnan(data.loc[ind, 'A_blade'])):
            if numpy.isnan(data.loc[ind, 'nff']):
                nff = int(round(hs['Tremie13'].mean_nff))
            else:
                nff = data.loc[ind, 'nff']
            num_leaf_bottom = data.loc[ind, 'rank']
            num_leaf_top = data.loc[ind, 'rank']-nff+1
            if num_leaf_top>sep_up_down:
                ff = fits['Tremie13'].form_factor()[1]
            else:
                ff = fits['Tremie13'].form_factor()[2]
            A = data.loc[ind, 'A_blade']
            if numpy.isnan(data.loc[ind, 'L_blade']):
                L = dim_fit.predict('L_blade', num_leaf_bottom, data.loc[ind, 'nff'])[0]
            else:
                L = data.loc[ind, 'L_blade']
            W_blades.append(get_estimate(A, L, ff))
            estimated.append(True)
        else:
            W_blades.append(data.loc[ind, 'W_blade'])
            estimated.append(False)
    data['W_blade'] = W_blades
    #data['blade_estimated'] = estimated
    #data = data[data['blade_estimated']==True]
    #data = data.drop('blade_estimated', 1)
    return data


def dimension_fits(tag='reference', reset=False):
    """ semi-manual fit of dimension data
    """
    parameters = reconstruction_parameters(tag, reset=reset)
    hsfit = HS_fit(tag)
    data = archidb.reconstruction_data()
    obs = data.Dimension_data
    card_scale = parameters.get('card_scale')
    # add W estimates for Tremie 13
    fit=pgen_ext.WheatDimensions(hsfit['Tremie13'],card_scale=card_scale['Tremie13'])
    fit.fit_dimensions(obs[obs['label'] == 'Tremie13'])#fit L to estimate width from area
    w_est = estimate_Wb_Tremie13(fit, obs[obs['label']=='Tremie13'].copy(), tag)
    obs = pandas.concat((obs[obs['label']!='Tremie13'],w_est))
    
    fits = {k: pgen_ext.WheatDimensions(hsfit[k], card_scale=card_scale[k]) for k in hsfit}
    
    for k in fits:
        fits[k].fit_dimensions(obs[obs['label'] == k])
    # Tremie13: use Hcol scales for sheath and internode
    for lab in ['Tremie13']:
        sc = fits[lab].scale['H_col']
        fits[lab].set_scale({'L_internode':sc, 'L_sheath':sc})
    # Mercia / Rht3 : use Hcol and sheath
    for lab in ['Mercia','Rht3']:
        fits[lab].set_scale({'L_internode':1})
        stem = fits[lab].predict('H_col') - fits[lab].predict('L_sheath')
        Li = [stem[0]] + numpy.diff(stem).tolist()
        pred = fits[lab].predict('L_internode')
        r = Li / pred
        sc = numpy.mean(r[(r > 0) & (numpy.isfinite(r))])
        fits[lab].set_scale({'L_internode':sc})
    return fits


def check_dimension(tag='reference', reset=False):

    fits = dimension_fits(tag, reset=reset)
    obs = archidb.dimensions_aggregated()
    
    archi_plot.dimension_plot(obs, fit = fits, dimension = 'L_blade')
    
    t13 = archidb.Dim_data()[archidb.Dim_data()['label'] == 'Tremie13'].copy()
    estimates_Wb_Tremie13 = archidb.dimensions_aggregated(estimate_Wb_Tremie13(fits['Tremie13'], t13, tag))
    archi_plot.dimension_plot(obs, fit = fits, dimension = 'W_blade',
                                estimates_Wb_Tremie13 = estimates_Wb_Tremie13)  
                                
    archi_plot.dimension_plot(obs, fit = fits, dimension = 'L_sheath')
    archi_plot.dimension_plot(obs, fit = fits, dimension = 'L_internode')
    archi_plot.dimension_plot(obs, fit = fits, dimension = 'H_col')
    
    archi_plot.dimension_plot_varieties(obs, fit = fits)
    
    
# leaf geometry

# Processing of raw database to get data reabale or adel

# TO DO : automatic selection of tol_med = f(ntraj), based on histogram approx of normal distribution


def leaf_trajectories(dfxy, dfsr, bins = [-10, 0.5, 1, 2, 3, 4, 10], ntraj = 10, tol_med = 0.1):
    """
    Return a dynamic leaf database compatible with adel and appropriate dynamic key selecting functions ({Lindex:{num_traj:{age_class:xy_dict}}})
    
    - dfxy and srdb are dataframe containing the data
    - bins is None if leaf are static or define age clas otherwise
    - ntraj and tol_med control  the number of trajectory to sample among leaves of a given Lindex and of given age that have mean_angle +/- tol*med
    """
                   
    import random
          
    dfxy['age'] = dfxy['HS'] - dfxy['rank'] + 1
    #cut intervalle de 0 a 1, etc.
    dfxy_cut = pandas.cut(dfxy.age, bins,labels=False)
    dfxy['age_class'] = dfxy_cut
    
    # use mean_angle to filter / group leaves
    mean_angle = dfxy.groupby('inerv').apply(lambda x: numpy.mean(abs(numpy.arctan2(numpy.diff(x['y'].values),numpy.diff(x['x'].values)))))

    #filter leaves that are above/below med*tol
    def filter_leaves(x):
        angles = mean_angle[set(x['inerv'])]
        med = angles.median()
        valid_angles = angles[(angles >= (1 - tol_med) * med) & (angles <= (1 + tol_med) * med)]
        return x[x['inerv'].isin(set(valid_angles.index))]
    
    validxy = dfxy.groupby(('Lindex','age_class'), group_keys=False).apply(filter_leaves)
    grouped = validxy.groupby(('Lindex','age_class'))
    
    # build trajectories
    trajectories = {k:[] for k in set(validxy['Lindex'])}
    for i in range(ntraj):
        for k in set(validxy['Lindex']):
            trajectories[k].append({})
            for t in set(validxy['age_class']):
                x = grouped.get_group((k,t))
                trajectories[k][i][t] = x.loc[x['inerv'] == random.sample(set(x['inerv']),1),['x','y']].to_dict('list')
    
    srdb = {k:v.loc[:,['s','r']].to_dict('list') for k, v in dfsr.groupby('Lindex')}
    
    return trajectories, srdb, bins

# azimuth, lindex as function of rankclass (class 1 basal leaves, class 2 upper leaves)    


def geoLeaf(nlim=4,dazt=60,dazb=10, Lindex_base = 1, Lindex_top = 2):
    """ generate geoLeaf function for Adel """
    rcode = """
    geoLeaf <- list(
     Azim = function(a,n,nf) {{
            ntop = nf - n
            ifelse(ntop <= {ntoplim:d},
            180 + {dazTop:.2f} * (runif(1) - .5),
            180 + {dazBase:.2f} * (runif(1) - .5))
            }},
     Lindex = function(a,n,nf) {{
              ntop = nf - n
              ifelse(ntop <= {ntoplim:d}, {top_class:d}, {base_class:d})
              }}
              )
    """
    return rcode.format(ntoplim = nlim, dazTop = dazt, dazBase = dazb, top_class=Lindex_top, base_class= Lindex_base)


def _addLindex(dfxy, dfsr): 
    dfxy['Lindex'] = 1
    dfxy.loc[dfxy['ranktop'] <= 4,'Lindex'] = 2
    dfsr['Lindex'] = dfsr['rankclass']
    return dfxy, dfsr


def leaf_fits(bins=[-10, 0.5, 1, 2, 3, 4, 10], ntraj=10, tol_med=0.1, disc_level=7):
    d = {}
    gL = geoLeaf()
    for k in ['Mercia','Rht3', 'Tremie12', 'Tremie13']:
        dfxy, dfsr = _addLindex(*archidb.leaf_curvature_data(k))
        xy, sr, bins = leaf_trajectories(dfxy, dfsr , bins = bins, ntraj = ntraj, tol_med = tol_med)
        d[k] = Leaves(xy, sr, geoLeaf=gL, dynamic_bins = bins, discretisation_level = disc_level)
    return d


def median_leaf_fits(xydata, sr_data, disc_level=7, top_leaves=3):
    gL = geoLeaf(nlim=top_leaves)   
    trajs,bins = xydata
    sr_data['Lindex'] = sr_data['rankclass']
    srdb = {k:v.loc[:,['s','r']].to_dict('list') for k, v in sr_data.groupby('Lindex')}
    return Leaves(trajs, srdb, geoLeaf=gL, dynamic_bins = bins, discretisation_level = disc_level)


def leafshape_fits(tag='reference'):
    parameters = reconstruction_parameters(tag)
    data = archidb.reconstruction_data()
    median_leaf = parameters.get('median_leaf')
    if not all(median_leaf.values()):
        raise NotImplementedError('leaf_fits not operational anymore')
        #use leaf_fits in that case
    top_leaves = parameters.get('top_leaves')
    xydata = parameters.get('xy_data')
    disc_level = parameters.get('disc_level')
    return {k: median_leaf_fits(data.xy_data[xydata[k]], data.sr_data[k], disc_level=disc_level[k], top_leaves=top_leaves[k]) for k in median_leaf}

 
# Attention pour rht3 on veut geoleaf qui retourne 1 (= feuille du bas) quelque soit le rang !
    # creation de la colonne age et du selecteur d'age 


class EchapReconstructions(object):
    def __init__(self, tag='reference'):
        self.tag = tag
        self.pars = reconstruction_parameters(tag, reset=True)
        # if reset_data:
        #     d = archidb.reconstruction_data(reset=True)
        #     reset_data = False
        self.pgen_base = self.pars['pgen_base']
        self.HS_fit = HS_fit(tag, reset=True)
        # parameters already reset : can use cache
        self.density_fits = density_fits(tag)
        self.axepop_fits = axepop_fits(tag)
        self.dimension_fits = dimension_fits(tag)
        self.GL_fits = GL_fits(tag)
        self.leaves = leafshape_fits(tag)

    def get_reconstruction(self, name='Mercia', nplants=30, nsect=3, seed=1,
                           aborting_tiller_reduction=1, aspect='square',
                           ssipars={'r1': 0.07, 'ndelsen': 3}, **kwds):
        ''''''

        if name == 'Rht3':
            incT = 75
            dep = 10
        else:
            incT = 60
            dep = 7

        d = self.density_fits[name]
        stand = AgronomicStand(sowing_density=d['sowing_density'],
                               plant_density=d['density_at_emergence'],
                               inter_row=d['inter_row'], noise=0.04,
                               density_curve_data=d['density_table'])
        # n_emerged, domain, positions, area = stand.stand(nplants, aspect)
        n_emerged = nplants  # adel uses smart stand
        axp = self.axepop_fits[name]
        plants = axp.plant_list(n_emerged)
        pgen = pgen_ext.PlantGen(HSfit=self.HS_fit[name],
                                 GLfit=self.GL_fits[name],
                                 Dimfit=self.dimension_fits[name],
                                 base_config=self.pgen_base,
                                 adel_pars=self.pars['adel_pars'])
        axeT, dimT, phenT = pgen.adelT(plants)
        axeT = axeT.sort_values(['id_plt', 'id_cohort', 'N_phytomer'])
        devT = devCsv(axeT, dimT, phenT)

        leaves = self.leaves[name]

        run_adel_pars = self.pars['adel_pars']

        return AdelWheat(nplants=nplants, nsect=nsect, devT=devT, stand=stand,
                         seed=seed, sample='sequence', leaves=leaves,
                         aborting_tiller_reduction=aborting_tiller_reduction,
                         aspect=aspect, incT=incT, dep=dep,
                         run_adel_pars=run_adel_pars, **kwds)

    def save(self, filename):

        with open(str(filename), 'wb') as output:
            pickle.dump(self, output)
    
    @staticmethod
    def load(filename):
        with open(filename,'rb') as saved:
            return pickle.load(saved,encoding='bytes')


def echap_reconstructions(tag='reference', reset=False):
    filename = cache_reconstruction_path(tag) / 'EchapReconstructions.pckl'
    if not reset:
        try:
            with open(filename) as saved:
                return pickle.load(saved)
        except:
            pass
    echap = EchapReconstructions(tag)
    echap.save(filename)
    return echap


def soisson_reconstruction(nplants=30, sowing_density=250., plant_density=250.,
                           inter_row=0.15, nsect=3, seed=1):
    stand = AgronomicStand(sowing_density=sowing_density,
                           plant_density=plant_density, inter_row=inter_row,
                           noise=0.04, density_curve_data=None)
    n_emerged = nplants
    m = pgen_ext.TillerEmission(
        primary_tiller_probabilities={'T1': 1., 'T2': 0.5, 'T3': 0.5,
                                      'T4': 0.3})  # From Bertheloot Soissons N+
    axp = pgen_ext.AxePop(Emission=m)  # With 11 and 12 it's fine
    plants = axp.plant_list(n_emerged)
    hs_fit = HS_fit()
    GLfit = GL_fits()['Mercia']
    Dimfit = dimension_fits()['Mercia']
    Dimfit.scale = {k: v * 1.15 for k, v in
                    Dimfit.scale.items()}  # Seen on Soisson 2010 compared to Mercia 2010
    pgen = pgen_ext.PlantGen(HSfit=hs_fit['Mercia'], GLfit=GLfit, Dimfit=Dimfit)
    axeT, dimT, phenT = pgen.adelT(plants)
    axeT = axeT.sort_values(['id_plt', 'id_cohort', 'N_phytomer'])
    devT = devCsv(axeT, dimT, phenT)
    leaves = leafshape_fits()['Mercia']  # TODO Create and Take Soisson
    return AdelWheat(nplants=nplants, nsect=nsect, devT=devT, stand=stand,
                     seed=seed, sample='sequence', leaves=leaves)
    
# checks consistency adel/fits

def check_adel():
    e=echap_reconstructions()
    adel = e.get_reconstruction()
    df = adel.checkAxeDyn()
    fit = e.density_fits['Mercia']['density_table']
    ax=fit.plot('TT','density',style='-')
    df.plot('TT','nbplants',style='--',ax=ax)
