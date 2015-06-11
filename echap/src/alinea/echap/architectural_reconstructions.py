""" Reconstruction of Wheat for Boigneville data
"""
import pandas
import numpy
import math
import scipy.stats as stats
from scipy.interpolate import interp1d

from copy import deepcopy
try:
    import cPickle as pickle
except:
    import pickle


from alinea.adel.astk_interface import AdelWheat
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.WheatTillering import WheatTillering
from alinea.adel.AdelR import devCsv

import alinea.echap
from openalea.deploy.shared_data import shared_data

import alinea.echap.architectural_data as archidb
import alinea.echap.architectural_reconstructions_plot as archi_plot

import alinea.adel.plantgen_extensions as pgen_ext
run_plots = False # prevent ipython %run to make plots


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
    return {k:value for k in ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')}
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
    # delay between emergence of flag leaf on mainstem and flag leaf emergence on cohorts
    pars['dTT_cohort'] = pdict({'first': 60, 'increment': 10})
    #
    # Green Leaves (ssi) = f(HS)
    #
    pars['GLpars'] = pdict({'GL_start_senescence':4.4, 'GL_bolting':2.5, 'GL_flag': 4.6, 'n_elongated_internode': 3.7})
    pars['GLpars']['Tremie12']={'GL_start_senescence':4.4, 'GL_bolting':2.5, 'GL_flag': 5, 'n_elongated_internode': 4}
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
    # n_elongated internode is used to control HS of tiller regression start for the mean plant (start = nff - n_elongated_internode + delta_reg)
    # if original plantgen model behaviour is desired (true decimal number of elongated internodes computed from dimension data), this value should be set to None
    pars['n_elongated_internode'] = {'Mercia':3.2, 'Rht3':3.1, 'Tremie12':4.2, 'Tremie13':4}
    # reduction of emmission
    # For Rht3, factor 2 of T5 compensates for absence of T6
    pars['emission_reduction'] = {'Mercia':None, 'Rht3':{'T5':2.}, 'Tremie12':None, 'Tremie13':{'T3':0.5, 'T4':0.5}}
    # damages to tillers
    # when indicates the period (Haun stages) during which damage effectively occur (ie axes disappear).
    # when will be adapted for each tiller so that damaged tillers can emerge before dying    
    # damage indicates the proportion of tiller lost during the period (as a relative proportion of tillers emited)
    # natural regression participate to damage counts, and damaged tillers participate to natural regression (ie damage to a future regressing tiller is not inducing a new regression)
    #pars['tiller_damages'] = None
    pars['tiller_damages'] = {'Mercia':{'when':[6,13],'damage':{t:0.25 for t in ('T1','T2','T3')}}, 
                              'Rht3':{'when':[7,13],'damage':{t:0.3 for t in ('T1','T2','T3')}}, 
                              'Tremie12':{'when':[4.9,5.1],'damage':{t:0.5 for t in ['T3']}},
                              'Tremie13': None}
    #
    # Leaf geometry
    #--------------
    # flag for forcing median leaf shape = f(age) model instead of random shape sampling
    pars['median_leaf'] = pdict(True)
    # number of top leaves to be considered erectophyl
    #pars['top_leaves'] = {'Mercia':3, 'Rht3':2, 'Tremie12':4, 'Tremie13':3} # values optimised for echap report december 2014
    pars['top_leaves'] = {'Mercia':4, 'Rht3':4, 'Tremie12':4, 'Tremie13':4}
    #
    #
    return pars

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
        if lab == 'Tremie12':
            dat = g.get_group(lab)
            dat = dat.loc[dat['HS'] < dat['nff'],('TT','HS')].dropna()
            datd = gd.get_group(lab)
            datd = datd.loc[(datd['HS'] < datd['nff']) | (numpy.isnan(datd['nff'])),('TT','HS')].dropna()
            datd = datd.loc[datd['HS'] < 10,:]
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
    hs_fits = {k: pgen_ext.HaunStage(hs_ms['slope'][k], TTem[k], std_TTem[k], hs_ms['nff'][k], dTTnff[k], dTT_cohort[k]) for k in dTT_cohort}
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


if run_plots:
    parameters = reconstruction_parameters()
    fits = HS_fit(reset=True) 
    obs = archidb.validation_data()
    archi_plot.haun_stage_plot(obs.haunstage, fits)


#
# --------------------------------------------------------------- Fit GL = f(HS since flag)
#
def GL_fit(reset_data=False, **parameters):
    """ Regression for Green leaf data
    """
    GLpars = parameters.get('GLpars')
    fits = {k:pgen_ext.GreenLeaves(**GLpars[k]) for k in GLpars}
    #  TO DO : fit a (GreenLeaves method) with GL data
    df_GL_obs_nff, df_GL_est_nff, df_GL_obs_global, df_GL_est_global = archidb.green_leaves_aggregated(HS_fit())
    obs = pandas.concat([df_GL_obs_global, df_GL_est_global])
    for k in fits:
        if k=='Tremie12':
            obsk=obs[obs['label']=='Tremie12']
        else:
            obsk=obs[obs['label']!='Tremie12']
            x = obsk['HS'].astype(float).values-obsk['HSflag'].astype(float).values
            obsk=obsk[x<10]
        fits[k].fit_a(obsk['HS'].astype(float).values-obsk['HSflag'].astype(float).values,
                      obsk['GL_mean'].astype(float).values)
    return fits
    

if run_plots:
    parameters = reconstruction_parameters()
    fits = GL_fit(**parameters) 
    obs = archidb.green_leaves_aggregated(HS_fit())
    # Compare varieties
    archi_plot.green_leaves_plot_mean(obs, fits, HS_fit())
    # Compare nff
    archi_plot.green_leaves_plot(obs, fits, HS_fit())

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
    fits = {k: {'HS':[0,20],'density':[mean_plant_density[k]] * 2} for k in mean_plant_density}
    #
    if damages is not None:
        for k in damages:
            if damages[k] is not None:
                density_at_emergence = round(mean_plant_density[k] * 1. / (1 - damages[k]['damage']))
                density_at_harvest = mean_plant_density[k]
                fits[k] = {'HS':[0] + damages[k]['when'] + [30], 'density': [density_at_emergence] * 2 + [density_at_harvest] * 2}
    
    density_fits = {k:{'sowing_density': pdb[k]['sowing_density'], 
                       'inter_row': pdb[k]['inter_row'],
                       'density_table': pandas.DataFrame(fits[k])} for k in fits}

                       
    adjust = parameters.get('density_tuning')
    if adjust is not None: 
        for k in adjust:
            if adjust[k] is not None:
                df = density_fits[k]['density_table']
                df['density'] *= adjust[k]
                density_fits[k]['sowing_density'] *= adjust[k]
                density_fits[k]['inter_row'] /= math.sqrt(adjust[k])

    for k in density_fits:
        df = density_fits[k]['density_table']
        df['TT'] = HS_converter[k].TT(df['HS'])
        density_fits[k]['density_at_emergence'] = df['density'].iloc[0]
        density_fits[k]['density_at_harvest'] = df['density'].iloc[-1]
        density_fits[k]['hs_curve'] = interp1d(df['HS'], df['density'])
        density_fits[k]['TT_curve'] = interp1d(df['TT'], df['density'])
                
    return density_fits


#
if run_plots:
    parameters = reconstruction_parameters()
    fits = density_fits(**parameters)
    obs = archidb.validation_data()
    archi_plot.density_plot(obs.PlantDensity, fits, HS_fit())
#
# --------------------------------------------------------------- Fit Tillering
#
#
# fitting functions

                          
def _axepop_fit(tdb, delta_stop_del, n_elongated_internode, max_order, tiller_damages, std_em):
    ms_nff_probas = tdb['nff_probabilities']
    ears_per_plant = tdb['ears_per_plant']
    primary_emission = tdb['emission_probabilities']
    emission = pgen_ext.TillerEmission(primary_emission)
    regression = pgen_ext.TillerRegression(ears_per_plant, n_elongated_internode, delta_stop_del)
    return pgen_ext.AxePop(ms_nff_probas, emission, regression, tiller_damages = tiller_damages, max_order=max_order, std_em=std_em)
#
# Fits
def axepop_fits(reset_data=False, **parameters):

    # tillering model parameters
    delta_stop_del = parameters.get('delta_stop_del')
    n_elongated_internode = parameters.get('n_elongated_internode')
    if n_elongated_internode is None: # to do :use wheat Tillering decimal internode number and dimension fits to compute it as did plantgen
        raise NotImplementedError("Not yet implemented")
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
                    
    fits={k: _axepop_fit(tdb[k], delta_stop_del=delta_stop_del[k],n_elongated_internode=n_elongated_internode[k], max_order=max_order, tiller_damages=tiller_damages[k], std_em = std_em[k]) for k in tdb}
    
    return fits


if run_plots:
    # populate with fits
    parameters = reconstruction_parameters()
    fits = axepop_fits(**parameters)
    obs = archidb.validation_data()
    
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

def all_scan():
    df_obs_all = pandas.DataFrame()
    for name in ['Tremie12','Tremie13']:
        df_obs = archidb.treatment_scan(name)
        df_obs['var'] = name
        df_obs_all = df_obs_all.append(df_obs)
    df_obs_all = df_obs_all[df_obs_all['moyenne id_Feuille'] <= df_obs_all['HS']]
    return df_obs_all

if run_plots:
    archi_plot.dimension_plot(archidb.dimensions_data(), archidb.dimension_fits(), leaf_fits(), all_scan(), archidb.blade_dimensions_MerciaRht3_2009_2010())
    
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
    dfxy['Lindex'][dfxy['ranktop'] <= 4] = 2
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

def median_leaf_fits(disc_level=7, top_leaves={'Mercia':3, 'Rht3':2, 'Tremie12':4, 'Tremie13':3}):

    gL = {k:geoLeaf(nlim=3) for k in ['Mercia','Tremie13']}
    gL['Mercia'] = geoLeaf(nlim=top_leaves['Mercia']) 
    gL['Rht3'] = geoLeaf(nlim=top_leaves['Rht3']) 
    gL['Tremie12'] = geoLeaf(nlim=top_leaves['Tremie12'])    
    gL['Tremie13'] = geoLeaf(nlim=top_leaves['Tremie13']) 
    
    trajs,bins = archidb.median_leaf_trajectories()
    
    d={}
    for k in ['Mercia','Rht3', 'Tremie12', 'Tremie13']:
        _, dfsr = _addLindex(*archidb.leaf_curvature_data(k))
        srdb = {k:v.ix[:,['s','r']].to_dict('list') for k, v in dfsr.groupby('Lindex')}
        d[k] = Leaves(trajs, srdb, geoLeaf=gL[k], dynamic_bins = bins, discretisation_level = disc_level)
    return d
 
# Attention pour rht3 on veut geoleaf qui retourne 1 (= feuille du bas) quelque soit le rang !
    # creation de la colonne age et du selecteur d'age 
    
            
# Standard echap reconstruction protocol for plant development and axis dynamic       
#def echap_development(nplants, pgen, primary_proba, tdata, pdata, dTT_stop)
    #ici include mortality



class EchapReconstructions(object):
    
    def __init__(self, reset_data=False):
        self.pars = reconstruction_parameters()
        
        self.pgen_base = self.pars['pgen_base']
        self.HS_fit = HS_fit(reset=True, reset_data=reset_data)
        converter = HS_fit(reset=True, reset_data=reset_data)
        self.density_fits = density_fits(HS_converter=converter, reset_data=reset_data, **self.pars)
        self.axepop_fits = axepop_fits(reset_data=reset_data, **self.pars)
        self.dimension_fits = archidb.dimension_fits()
        self.GL_fits = GL_fit(reset_data=reset_data, **self.pars)
        median_leaf = self.pars['median_leaf']
        top_leaves = self.pars['top_leaves']
        if median_leaf:            
            self.leaf_fits = median_leaf_fits(top_leaves=top_leaves)
        else:
            self.leaf_fits = leaf_fits()

   
    def get_reconstruction(self, name='Mercia', nplants=30, nsect=3, seed=1, disc_level=7, aborting_tiller_reduction=1, aspect = 'square', stand_density_factor = {'Mercia':1, 'Rht3':1, 'Tremie12':1, 'Tremie13':1}, dimension=1, ssipars={'r1':0.07,'ndelsen':3},**kwds):
        '''stand_density_factor = {'Mercia':0.9, 'Rht3':1, 'Tremie12':0.8, 'Tremie13':0.8}, **kwds)'''
                
        run_adel_pars = {'rate_inclination_tiller': 15, 'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1}
        if name == 'Rht3':
            incT=75
            dep=10
        else:
            incT=60
            dep=7
            
        d = self.density_fits[name]
        stand = AgronomicStand(sowing_density=d['sowing_density'], plant_density=d['density_at_emergence'], inter_row=d['inter_row'], noise=0.04, density_curve = d['TT_curve'])       
        #n_emerged, domain, positions, area = stand.stand(nplants, aspect)
        n_emerged = nplants#adel uses smart stand       
        axp = self.axepop_fits[name]
        plants = axp.plant_list(n_emerged)
        pgen = pgen_ext.PlantGen(HSfit = self.HS_fit[name], GLfit = self.GL_fits[name])
        axeT, dimT, phenT = pgen.adelT(plants)
        axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])
        devT = devCsv(axeT, dimT, phenT)
                    
        leaves = self.leaf_fits[name] 
        
        std_emergence = self.HS_fit[name].std_TT_hs_0
        
        return AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , std_em = std_emergence, seed=seed, sample='sequence', leaves = leaves, aborting_tiller_reduction = aborting_tiller_reduction, aspect = aspect,incT=incT, dep=dep, run_adel_pars = run_adel_pars, **kwds)
    
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
    
    

    