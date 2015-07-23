""" Reconstruction of Wheat for Boigneville data
"""
import pandas
import numpy
import math
import scipy.stats as stats
from scipy.interpolate import interp1d

#import warnings
#warnings.simplefilter('ignore', numpy.RankWarning)

from copy import deepcopy
try:
    import cPickle as pickle
except:
    import pickle


from alinea.adel.astk_interface import AdelWheat
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.AdelR import devCsv

import alinea.echap
from openalea.deploy.shared_data import shared_data

import alinea.echap.architectural_data as archidb
import alinea.echap.architectural_reconstructions_plot as archi_plot

import alinea.adel.plantgen_extensions as pgen_ext

run_plots = False # prevent ipython %run to make plots
reset_data = False# control the ipython %run behavior concerning data + HSfit dependent data

    

class NotImplementedError(Exception):
    pass

#generic function to be moved to adel
 
'''
def sen(pgen):
    """ Creates devT tables from plantgen dict
    """
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data, plantgen2adel
    from alinea.adel.AdelR import devCsv
    
    _, _, _, _, _, _, _, HS_GL_SSI_T, _, _, _ = gen_adel_input_data(**pgen)
    
    # verif de la sortie pour obtenir le graph de SSI
    # HS_GL_SSI_T.to_csv('C:/Users/Administrateur/openaleapkg/echap/test/HS_GL_SSI_T_test.csv')
    
    return HS_GL_SSI_T
    
def plantgen_as_dict(inputs, dynT, dimT):
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
    
'''
   
def setAdel(**kwds):    
     
    adel = AdelWheat(**kwds)    
    return adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants

#---------- reconstructions

#
# Fitting parameters
#
def pdict(value):
    """ create a parameter dict for all echap cultivar with value
    """
    return {k:deepcopy(value) for k in ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')}
#
def reconstruction_parameters():
    """ Manual parameters used for fitting Echap experiments
    """
    pars={}
    #
    # PlantGen parameters
    #--------------------
    pars['pgen_base'] = {'TT_hs_break':0.0,
                         # delay  hs->collar apprearance (pyllochronic units)
                        'inner_params': {'DELAIS_PHYLL_HS_COL_NTH' : 0.6 - 0.5 / 1.6}}
    #
    # Haun Stage = f(TT), convergence between axis
    #---------------------------------------------
    # delay between emergence of flag leaf on mainstem and flag leaf emergence on cohorts (60 is for Maxwell that has quite large desynchronisation, 30 may be more realistic)
    pars['dTT_cohort'] = pdict({'first': 30, 'increment': 10})
    #
    # Green Leaves (ssi) = f(HS)
    #
    # in the original pgen, dHS_bolting is estimated as the decimal number of elongated internodes
    pars['GLpars'] = pdict({'GL_start_senescence':4.4, 'GL_bolting':2.2, 'GL_flag': 4.4, 'dHS_bolting': 3.7, 'curvature' : -0.002})
    pars['GLpars']['Tremie12'].update({'GL_bolting':1.7, 'GL_flag': 5, 'dHS_bolting': 4, 'curvature' : -0.007})
    pars['GLpars']['Tremie13'].update({'GL_bolting':2.1})
    pars['GLpars']['Rht3'].update({'GL_bolting':2.3})
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
    pars['density_tuning'] = None
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
                              'Tremie12':{'when':[4.9,5.1],'damage':{t:0.2 for t in ['T3']}},
                              'Tremie13': None}
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
    # number of top leaves to be considered erectophyl
    #pars['top_leaves'] = {'Mercia':3, 'Rht3':2, 'Tremie12':4, 'Tremie13':3} # values optimised for echap report december 2014
    pars['top_leaves'] = {'Mercia':4, 'Rht3':4, 'Tremie12':4, 'Tremie13':4}
    #
    pars['disc_level'] = pdict(7)
    #
    #
    return pars

parameters = reconstruction_parameters()
#
# --------------------------------------------------------------- Fit HS = f(TT)
#

def linreg_df(x,y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return pandas.Series({'slope':slope, 'intercept':intercept, 'r2': r_value**2, 'std_err':std_err})
   
def fit_HS(reset_data=False, **parameters):
    """ Linear regression for phyllochron data
    """
    data = archidb.reconstruction_data(reset=reset_data)
    tagged = data.Pheno_data['archi_tagged']
    sampled = data.Pheno_data['archi_sampled']
    dTT_cohort = parameters.get('dTT_cohort')
    
    # HS of mean plant 
    g = tagged.groupby('label')
    def _fit(df):
        dat = df.loc[df['HS'] < df['nff'],('TT','HS')].dropna()
        nff = df.groupby('N').agg('mean')['nff'].mean()
        n = df.groupby('N').agg('mean')['nff'].count()
        res = linreg_df(dat['TT'],dat['HS'])
        res['nff'] = nff
        res['n'] = n
        return res
    hs_ms = g.apply(_fit)
    # complete fit if destructive sample available
    # for Tremie 12, keep only first sample
    # for Tremie 13, sampled plants are late / tagged plant : not used
    gd = sampled.groupby('label')
    for lab in gd.groups:
#        if lab == 'Tremie12':
        dat = g.get_group(lab)
        dat = dat.loc[dat['HS'] < dat['nff'],('TT','HS')].dropna()
        datd = gd.get_group(lab)
        datd = datd.loc[(datd['HS'] < datd['nff']) | (numpy.isnan(datd['nff'])),('TT','HS')].dropna()
        #datd = datd.loc[datd['HS'] < 10,:]
        datreg = pandas.concat((dat,datd))
        res = linreg_df(datreg['TT'],datreg['HS'])
        hs_ms.loc[lab,'intercept':'std_err'] = res
    # TTem per label
    TTem = - hs_ms['intercept'] / hs_ms['slope']
    # Flag leaf delay per nff
    TT_mean_flag = (hs_ms['nff'] - hs_ms['intercept']) / hs_ms['slope']
    g = tagged.groupby(('label','nff'))
    hs_nff = g.apply(_fit)
    TT_nff_flag = (hs_nff['nff'] - hs_nff['intercept']) / hs_nff['slope']
    TT_flag = TT_nff_flag.reset_index().rename(columns={0:'flag'}).merge(TT_mean_flag.reset_index().rename(columns={0:'mean_flag'}))
    TT_flag['dTT_flag'] = TT_flag['flag'] - TT_flag['mean_flag']
    TT_flag = TT_flag.merge(hs_nff['n'].reset_index().rename(columns=({0:'n'})))
    TT_flag = TT_flag.merge(hs_ms['nff'].reset_index().rename(columns=({'nff':'mean_nff'})))
    TT_flag = TT_flag.merge(TTem.reset_index().rename(columns=({0:'mean_TTem'})))
    dTTnff  = TT_flag.groupby('label').apply(lambda x: numpy.average(x['dTT_flag'] * 1. / (x['nff'] - x['mean_nff']), weights = x['n']))
    # residual variabilityy of emergenece between plants
    # mean slopes per nff: 
    dat = TT_flag.set_index(['label', 'nff'], drop=False)
    a_nff = dat['nff'] / (dat['flag'] - dat['mean_TTem']) 
    #estimate dTTEm per plant wyth modeled slope
    g = tagged.groupby(('label', 'N'))
    def _TTem(df):
        res = None
        dat = df.loc[df['HS'] < df['nff'],('TT','HS')].dropna()
        nff = df['nff'].values[0]
        label = df['label'].values[0]
        if not numpy.isnan(nff):
            a = a_nff[(label, nff)]
            b = numpy.mean(dat['HS'] - a * dat['TT'])
            res = -1. * b / a
        return res
    TTem_p = g.apply(_TTem)
    # standard deviation
    std_TTem = TTem_p.reset_index().groupby('label').agg('std').loc[:,0]
  
    # HS fits
    hs_fits = {k: pgen_ext.HaunStage(1. / hs_ms['slope'][k], TTem[k], std_TTem[k], hs_ms['nff'][k], dTTnff[k] / hs_ms['slope'][k], dTT_cohort[k]) for k in dTT_cohort}
    return hs_fits
    
    
def HS_fit(reset=False, reset_data=False):
    """ Handle fitted phyllochron persitent object
    """
    filename = str(shared_data(alinea.echap)/'HS_fit.pckl')
    if not reset:
        try:
            with open(filename) as input:
                return pickle.load(input)
        except:
            pass
    parameters = reconstruction_parameters()
    fit = fit_HS(reset_data=reset_data, **parameters)
    with open(filename,'w') as output:
        pickle.dump(fit, output)
    return fit
#
if reset_data:
    HSfit = HS_fit(reset=True, reset_data=True)
    vdata = archidb.validation_data(reset=True, HS_fit=HSfit)
    rdata = archidb.reconstruction_data(reset=True)
else:
    HSfit = HS_fit()
    vdata = archidb.validation_data()
    rdata = archidb.reconstruction_data()
    

if run_plots:
    fits = HSfit
    obs = vdata.haun_stage
    archi_plot.haun_stage_plot(obs, fits)


#
# --------------------------------------------------------------- Fit GL = f(HS since flag)
#
def GL_fits(hsfit, **parameters):
    """ Regression for Green leaf data
    """
    GLpars = parameters.get('GLpars')
    
    #df_GL_obs_nff,  df_GL_obs_global = data.green_leaves
    #df_GL_est_nff,  df_GL_est_global = data.green_leaves_estimated   
    #obs = pandas.concat([df_GL_obs_global, df_GL_est_global])
    ## compute curvature using Tremie12 data
    #obsT12 = obs[obs['label']=='Tremie12']
    #fitT12 = pgen_ext.GreenLeaves(**GLpars['Tremie12'])
    #fitT12.fit_a(obsT12['HS'].astype(float).values - obsT12['HSflag'].astype(float).values,
    #                  obsT12['GL_mean'].astype(float).values)
    #print fitT12.a
    #
    return {k:pgen_ext.GreenLeaves(hsfit[k], **GLpars[k]) for k in GLpars}

    
if run_plots:
    gl_obs_nff, gl_obs = vdata.green_leaves
    gl_est_nff, gl_est = vdata.green_leaves_estimated
    # Compare varieties
    obs = (gl_obs, gl_est)
    fits = GL_fits(HSfit,**parameters)
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

def density_fits(HS_converter=None, reset_data=False, **parameters):
    """
    Manual fit of plant density based on mean plant density and estimate of plant density at harvest
    """

    if HS_converter is None:
        HS_converter = HS_fit()
    
    data = archidb.reconstruction_data(reset=reset_data)
    pdb = data.Plot_data
    tdb = data.Tillering_data

    damages = parameters.get('plant_damages')
    
    #compile density data
    density_data = {}
    for k in pdb:
        hs =[]
        d = []
        if 'plant_density_at_emergence' in pdb[k]:
            hs.append(HS_converter[k](pdb[k]['TT_date'][pdb[k]['code_date']['emergence']]))
            d.append(pdb[k]['plant_density_at_emergence'])
        if 'plant_density' in pdb[k]:
            for date,densities in pdb[k]['plant_density'].iteritems():
                hs.extend([HS_converter[k](pdb[k]['TT_date'][date])] * len(densities))
                d.extend(densities)
        if 'raw_ear_density_at_harvest' in pdb[k]:
            ear_d = numpy.array(pdb[k]['raw_ear_density_at_harvest']) / tdb[k]['ears_per_plant']
            hs_harvest = HS_converter[k](pdb[k]['TT_date'][pdb[k]['code_date']['harvest']])
            hs.extend([hs_harvest] * len(ear_d))
            d.extend(ear_d.tolist())
        density_data[k] = pandas.DataFrame({'HS':hs, 'density':d})

    #estimate mean plant density or mean plant density after damages
    if damages is not None:
        for k in damages:
            if damages[k] is not None:
                density_data[k] = density_data[k].loc[density_data[k]['HS'] > max(damages[k]['when']),:]
    mean_plant_density = {k: round(density_data[k]['density'].mean()) for k in density_data}
    #
    fits = {k: {'HS':[-20,50],'density':[mean_plant_density[k]] * 2} for k in mean_plant_density}
    #
    if damages is not None:
        for k in damages:
            if damages[k] is not None:
                density_at_emergence = round(mean_plant_density[k] * 1. / (1 - damages[k]['damage']))
                density_at_harvest = mean_plant_density[k]
                fits[k] = {'HS':[-20] + damages[k]['when'] + [50], 'density': [density_at_emergence] * 2 + [density_at_harvest] * 2}
    
    fits = {k:{'sowing_density': pdb[k]['sowing_density'], 
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
if run_plots:
    fits = density_fits(**parameters)
    obs = vdata.PlantDensity
    archi_plot.density_plot(obs, fits, HSfit)
#
# --------------------------------------------------------------- Fit Tillering
#
#
# fitting functions

                          
def _axepop_fit(tdb, delta_stop_del, dHS_reg, max_order, tiller_damages, std_em):
    ms_nff_probas = tdb['nff_probabilities']
    ears_per_plant = tdb['ears_per_plant']
    primary_emission = tdb['emission_probabilities']
    emission = pgen_ext.TillerEmission(primary_emission)
    regression = pgen_ext.TillerRegression(ears_per_plant, dHS_reg, delta_stop_del)
    return pgen_ext.AxePop(ms_nff_probas, emission, regression, tiller_damages = tiller_damages, max_order=max_order, std_em=std_em)
#
# Fits
def axepop_fits(reset_data=False, **parameters):

    # tillering model parameters
    delta_stop_del = parameters.get('delta_stop_del')
    dHS_reg = parameters.get('dHS_reg')
    max_order = parameters.get('max_order')
    tiller_damages = parameters.get('tiller_damages')
    # Emission probabilities tuning
    #  
    data = archidb.reconstruction_data(reset=reset_data)
    tdb = deepcopy(data.Tillering_data)    
    # apply manual reductions
    emf = parameters.get('emission_reduction')
    if emf is not None:
        for k in emf:
            if emf[k] is not None:
                for T in emf[k]:
                    tdb[k]['emission_probabilities'][T] *= emf[k][T]

    # variability of emergence
    hs = HS_fit()
    std_em = {k:hs[k].std_TT_hs_0 for k in hs}
                    
    fits={k: _axepop_fit(tdb[k], delta_stop_del=delta_stop_del[k],dHS_reg=dHS_reg[k], max_order=max_order, tiller_damages=tiller_damages[k], std_em = std_em[k]) for k in tdb}
    
    return fits


if run_plots:
    # populate with fits
    fits = axepop_fits(**parameters)
    obs = vdata
    
    #1 graph tillering par var -> 4 graph sur une feuille
    archi_plot.multi_plot_tillering(obs.tillers_per_plant, fits, HS_fit())  
    
    #tallage primary pour les 4 var   
    archi_plot.tillering_primary(obs.tillers_per_plant, fits, HS_fit())
    
    #tallage total pour les 4 var    
    archi_plot.tillering_tot(obs.tillers_per_plant, fits, HS_fit())
   
   # primary emissions
    archi_plot.graph_primary_emission(obs.emission_probabilities)
#
# --------------------------------------------------------------- Fit dimensions
#

def estimate_Wb_Tremie13(dim_fit, data, sep_up_down=4):
    def get_estimate(A, L, ff):
        return float(A)/(L*ff)
    data = data.reset_index()
    data = data[data['label']=='Tremie13']
    hs = HS_fit()
    fits = leafshape_fits(reset_data=False,**parameters)
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

def dimension_fits(hsfit,reset_data=False, **parameters):
    """ semi-manual fit of dimension data
    """
    
    data = archidb.reconstruction_data(reset=reset_data)
    obs = data.Dimension_data
    card_scale = parameters.get('card_scale')
    # add W estimates for Tremie 13
    fit=pgen_ext.WheatDimensions(hsfit['Tremie13'],card_scale=card_scale['Tremie13'])
    fit.fit_dimensions(obs[obs['label'] == 'Tremie13'])#fit L to estimate width from area
    w_est = estimate_Wb_Tremie13(fit, obs[obs['label']=='Tremie13'].copy())
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


    
if run_plots:

    fits = dimension_fits(HSfit, **parameters)
    obs = vdata.dimensions
    
    archi_plot.dimension_plot(obs, fit = fits, dimension = 'L_blade')
    
    t13 = rdata.Dimension_data[rdata.Dimension_data['label'] == 'Tremie13'].copy()
    estimates_Wb_Tremie13 = archidb.dimensions_aggregated(estimate_Wb_Tremie13(fits['Tremie13'], t13))
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
                trajectories[k][i][t] = x.ix[x['inerv'] == random.sample(set(x['inerv']),1),['x','y']].to_dict('list')
    
    srdb = {k:v.ix[:,['s','r']].to_dict('list') for k, v in dfsr.groupby('Lindex')}
    
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
    srdb = {k:v.ix[:,['s','r']].to_dict('list') for k, v in sr_data.groupby('Lindex')}
    return Leaves(trajs, srdb, geoLeaf=gL, dynamic_bins = bins, discretisation_level = disc_level)
 
def leafshape_fits(reset_data=False, **parameters): 
    data = archidb.reconstruction_data(reset=reset_data)
    median_leaf = parameters.get('median_leaf')
    if not all(median_leaf.values()):
        raise NotImplementedError('leaf_fits not operational anymore')
        #use leaf_fits in that case
    top_leaves = parameters.get('top_leaves')
    disc_level = parameters.get('disc_level')
    return {k: median_leaf_fits(data.xy_data[k], data.sr_data[k], disc_level=disc_level[k], top_leaves=top_leaves[k]) for k in median_leaf}

 
# Attention pour rht3 on veut geoleaf qui retourne 1 (= feuille du bas) quelque soit le rang !
    # creation de la colonne age et du selecteur d'age 
    
            


class EchapReconstructions(object):
    
    def __init__(self, reset_data=False, pars = reconstruction_parameters()):
        self.pars = pars
        if reset_data:
            d = archidb.reconstruction_data(reset=True)
            reset_data = False
        self.pgen_base = self.pars['pgen_base']
        self.HS_fit = HS_fit(reset=True, reset_data=reset_data)
        self.density_fits = density_fits(HS_converter=self.HS_fit, reset_data=reset_data, **self.pars)
        self.axepop_fits = axepop_fits(reset_data=reset_data, **self.pars)
        self.dimension_fits = dimension_fits(self.HS_fit, reset_data=reset_data,**self.pars)
        self.GL_fits = GL_fits(self.HS_fit, **self.pars)
        self.leaves = leafshape_fits(reset_data=reset_data, **self.pars)

   
    def get_reconstruction(self, name='Mercia', nplants=30, nsect=3, seed=1, aborting_tiller_reduction=1, aspect = 'square', stand_density_factor = {'Mercia':1, 'Rht3':1, 'Tremie12':1, 'Tremie13':1}, dimension=1, ssipars={'r1':0.07,'ndelsen':3},**kwds):
        '''stand_density_factor = {'Mercia':0.9, 'Rht3':1, 'Tremie12':0.8, 'Tremie13':0.8}, **kwds)'''
                
        run_adel_pars = {'rate_inclination_tiller': 15, 'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1}
        if name == 'Rht3':
            incT=75
            dep=10
        else:
            incT=60
            dep=7
            
        d = self.density_fits[name]
        stand = AgronomicStand(sowing_density=d['sowing_density'], plant_density=d['density_at_emergence'], inter_row=d['inter_row'], noise=0.04, density_curve_data = d['density_table'])       
        #n_emerged, domain, positions, area = stand.stand(nplants, aspect)
        n_emerged = nplants#adel uses smart stand       
        axp = self.axepop_fits[name]
        plants = axp.plant_list(n_emerged)
        pgen = pgen_ext.PlantGen(HSfit = self.HS_fit[name], GLfit = self.GL_fits[name], Dimfit=self.dimension_fits[name])
        axeT, dimT, phenT = pgen.adelT(plants)
        axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])
        devT = devCsv(axeT, dimT, phenT)
                    
        leaves = self.leaves[name] 
       
        
        return AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample='sequence', leaves = leaves, aborting_tiller_reduction = aborting_tiller_reduction, aspect = aspect,incT=incT, dep=dep, run_adel_pars = run_adel_pars, **kwds)
    
    def save(self, filename):
        with open(filename, 'w') as output:
            pickle.dump(self, output)
    
def echap_reconstructions(reset=False, reset_data=False):
    filename = str(shared_data(alinea.echap)/'EchapReconstructions.pckl')
    if not reset:
        try:
            with open(filename) as input:
                return pickle.load(input)
        except:
            pass
    Echap = EchapReconstructions(reset_data=reset_data)
    Echap.save(filename)
    return Echap   

def soisson_reconstruction(nplants=30, sowing_density=250., plant_density=250., inter_row=0.15,
                            nsect=3, seed=1):
    stand = AgronomicStand(sowing_density=sowing_density,
                            plant_density=plant_density, 
                            inter_row=inter_row, noise=0.04, density_curve_data = None)       
    n_emerged = nplants
    axp = pgen_ext.AxePop() # With 11 and 12 it's fine
    plants = axp.plant_list(n_emerged)
    hs_fit = HS_fit()
    GLfit = GL_fits(hs_fit, **parameters)['Mercia']
    Dimfit = dimension_fits(hs_fit, **parameters)['Mercia']
    Dimfit.scale = {k:v*1.15 for k,v in Dimfit.scale.iteritems()} # Seen on Soisson 2010 compared to Mercia 2010
    pgen = pgen_ext.PlantGen(HSfit=hs_fit['Mercia'], GLfit=GLfit, Dimfit=Dimfit)
    axeT, dimT, phenT = pgen.adelT(plants)
    axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])
    devT = devCsv(axeT, dimT, phenT)
    leaves = leafshape_fits(**parameters)['Mercia'] # TODO Create and Take Soisson
    return AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , 
                    seed=seed, sample='sequence', leaves = leaves)
    
# checks consistency adel/fits

if run_plots:
    e=echap_reconstructions()
    adel = e.get_reconstruction()
    df = adel.checkAxeDyn()
    fit = e.density_fits['Mercia']['density_table']
    ax=fit.plot('TT','density',style='-')
    df.plot('TT','nbplants',style='--',ax=ax)
    
    
    # future deprecated    
def get_EchapReconstructions():
 
    filename = str(shared_data(alinea.echap)/'EchapReconstructions.pckl')
    f = open(filename)
    reconst = pickle.load(f)
    f.close()
    return reconst
    
    

    