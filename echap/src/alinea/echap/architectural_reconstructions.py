""" Reconstruction of Wheat for Boigneville data
"""
import pandas
import numpy
import math
import scipy.stats as stats
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
    # Plant Density
    #--------------
    #thining period (comment out if constant density is desired). Make density decrease linearly between density at emergence and density at harvest
    pars['thining_period'] = {'Mercia':[6,13], 'Rht3':[6,13], 'Tremie12': None, 'Tremie13': None}
    # density adjust factor(multiplicative effect)
    pars['adjust_density'] = None
    # Tillering
    #----------
    # max order of tillering or None
    pars['max_order'] = None
    # delay between stop of growth and disparition of a tiller
    pars['delta_stop_del'] = pdict(2.)
    # n_elongated internode is used to control HS of tiller regression start (start = nff - n_elongated_internode + delta_reg)
    # if original plantgen model behaviour is desired (true decimal number of elongated internodes computed from dimension data), this value should be set to None
    pars['n_elongated_internode'] = {'Mercia':3.5, 'Rht3':3., 'Tremie12':5, 'Tremie13':4}
    # reduction of emmission
    pars['emission_reduction'] = {'Mercia':None, 'Rht3':None, 'Tremie12':None, 'Tremie13':{'T3':0.5, 'T4':0.5}}
    # damages to tillers (HS values define the interval at which axes are deleted)
    #pars['tiller_damages'] = None
    pars['tiller_damages'] = {'Mercia':{'HS':[4,9],'density':[1,0.7],'tillers':('T1','T2','T3')}, 
                              'Rht3':{'HS':[4,9],'density':[1,0.6],'tillers':('T1','T2','T3')}, 
                              'Tremie12':{'HS':[4.9,5.1],'density':[1,0.2],'tillers':['T3']},
                              'Tremie13': None}
    #
    # Leaf geometry
    #--------------
    # number of top leaves to be considered erectophyl
    #pars['top_leaves'] = {'Mercia':3, 'Rht3':2, 'Tremie12':4, 'Tremie13':3} # values optimised for echap report december 2014
    pars['top_leaves'] = {'Mercia':4, 'Rht3':4, 'Tremie12':4, 'Tremie13':4}
    #
    #
    return pars


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


#
# Parametres communs / Pre-traitement des donnees
#
# Fit dynamique of median plant
#
# ---------------------------------------------------- Fitted HS = f(TT)
#
class HaunStage(object):
    """ Handle HaunStage = f (ThermalTime) fits
    """
    
    def __init__(self, a_cohort = 1. / 110., TT_hs_0 = 0):
        self.a_cohort = a_cohort
        self.TT_hs_0 = TT_hs_0
        
    def __call__(self, TT):# HS
        return (numpy.array(TT) - self.TT_hs_0) * self.a_cohort
        
    def TT(self, HS):
        return self.TT_hs_0 + numpy.array(HS) / self.a_cohort
        
    def TTem(self, TT):
        return numpy.array(TT) - self.TT_hs_0
        
    def phyllochron(self):
        return 1. / self.a_cohort
'''
if run_plots:
    archi_plot.dynamique_plot_nff(archidb.HS_GL_SSI_data(), dynamique_fits())   
'''    
   
def fit_HS():    
    # linear reg on HS data : on cherche fit meme pente
    hsd = {k:v['HS'] for k,v in archidb.HS_GL_SSI_data().iteritems()}
    nff = archidb.mean_nff()
    for k in nff:
        df = hsd[k]
        hsd[k] = df[df['mean_pond'] < int(nff[k])]
    #compute mean slope 
    sl = numpy.mean([stats.linregress(df['TT'], df['mean_pond'])[0] for df in hsd.values()])
    #dec at origin considering common slope
    dec = {k:numpy.mean(df['mean_pond'] - sl * df['TT']) / sl for k,df in hsd.iteritems()}
    hsdec = hsd
    for k in hsd:
        hsdec[k]['TTdec'] = hsd[k]['TT'] + dec[k]
    hsdec = reduce(lambda x,y: pandas.concat([x,y]), hsdec.values())
    # common fit
    a_cohort = stats.linregress(hsdec['TTdec'], hsdec['mean_pond'])[0]
    TTcol0 = {k:-numpy.mean(df['mean_pond'] - a_cohort * df['TT']) / a_cohort for k,df in hsd.iteritems()}
    return {k: HaunStage(a_cohort, TTcol0[k]) for k in TTcol0}
    

    
def HS_fit(reset=False):
    filename = str(shared_data(alinea.echap)/'HS_fit.pckl')
    if not reset:
        try:
            with open(filename) as input:
                return pickle.load(input)
        except:
            pass
    fit = fit_HS()
    with open(filename,'w') as output:
        pickle.dump(fit, output)
    return fit
#


if run_plots:
    HS_converter = fit_HS()
    archi_plot.dynamique_plot(archidb.HS_GL_SSI_data(), converter = HS_converter) # weighted (frequency of nff modalities) mean 

#
# Styrategie : pour premieres feuilles
#

def HS_GL_fits(HS_converter=None):
    if HS_converter is None:
        HS_converter = fit_HS()
    # common fit mercia, rht3, Tremie12
    conv = HS_converter['Tremie12']
    HS_ref = conv(archidb.GL_number()['Tremie12']['TT'])
    GL_ref = archidb.GL_number()['Tremie12']['mediane']
    conv = HS_converter['Tremie13']
    HS_T13 = conv(archidb.GL_number()['Tremie13']['TT'])
    GL_T13 = archidb.GL_number()['Tremie13']['mediane']
    HS = {'Mercia':HS_ref, 'Rht3': HS_ref, 'Tremie12':HS_ref, 'Tremie13':HS_T13}
    GL = {'Mercia':GL_ref, 'Rht3': GL_ref, 'Tremie12':GL_ref, 'Tremie13':GL_T13}
    nff = archidb.mean_nff()
    #coefs = {'n0': 4.4, 'hs_t1': 8.6}#n0=4.4 : Maxwell value
    coefs = {'n0': 4.7, 'hs_t1': 8.6} # newparameter to avoid too short grean area lifetime for disease
    #n1 = {'Mercia':1.9, 'Rht3': 1.9, 'Tremie12':1.9, 'Tremie13':2.8}
    n1 = {'Mercia':3, 'Rht3': 3, 'Tremie12':3, 'Tremie13':3}# newparameter to avoid too short grean area lifetime for disease (should be 3.5 but pb in fitting)
    n2 = {'Mercia':5, 'Rht3': 4.9, 'Tremie12':5, 'Tremie13':4.3}
    fits = {k:pgen_ext.GL_model(HS[k], GL[k], nff=nff[k], n2=n2[k],n1=n1[k],**coefs) for k in nff}
    #
    hsgl_fits = {k:{'HS': HS_converter[k], 'GL': fits[k]} for k in nff}
    return hsgl_fits
#
if run_plots:
    #
    data = archidb.HS_GL_SSI_data()
    fits = HS_GL_fits()
    #
    #archi_plot.dynamique_plot_GL_fits(data, fits , abs='HS')
    #
    #archi_plot.dynamique_plot_GL_fits(data, fits, abs='TT', obs=False)
    archi_plot.dynamique_plot_nff(data)
    # conc : GL dynamique identique whatever nff => on change plutot acohort par nff, tq TTem_t2 et TTem_t1 restent les memes.

def pars():
    sim_all={}
    e = EchapReconstructions()
    for name in ('Mercia','Rht3','Tremie12','Tremie13'):
        pars = e.get_pars(name=name, nplants=30)
        nff_lst = pars.keys()
        nff_proba = archidb.Tillering_data()[name]['nff_probabilities']

        for nff in nff_lst:
            select = pandas.Series([111,112,113,1111,1112,1113,10111,10112,10113])
            pars[nff]['HS_GL_SSI_T'] = pars[nff]['HS_GL_SSI_T'][pars[nff]['HS_GL_SSI_T'].id_phen.isin(select)]
            TT = pars[nff]['HS_GL_SSI_T']['TT']
        
        GL = sum([pars[k]['HS_GL_SSI_T']['GL'] * nff_proba[k] for k in nff_lst])
        HS = sum([pars[k]['HS_GL_SSI_T']['HS'] * nff_proba[k] for k in nff_lst])
        SSI = sum([pars[k]['HS_GL_SSI_T']['SSI'] * nff_proba[k] for k in nff_lst])
        df = pandas.DataFrame({'TT':TT, 'GL':GL, 'HS':HS, 'SSI':SSI})
        sim_all.update({name:df})
    return sim_all
    
if run_plots:
    HS_converter = fit_HS()
    archi_plot.dynamique_plot_sim(archidb.HS_GL_SSI_data(), pars(), converter = HS_converter)
 
#
# Fit plant density
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
    
    # in case one density only is used, choose mean density or density at harvest if missing
    plant_density = {k: round(pdb[k].get('mean_plant_density', 1. * pdb[k]['ear_density_at_harvest'] / tdb[k]['ears_per_plant'])) for k in pdb}
    #
    fits = {k: {'HS':[0,20],'density':[plant_density[k]] * 2} for k in plant_density}
    #
    thining = parameters.get('thining_period')
    if thining is not None:
        for k in thining:
            if thining[k] is not None:
                density_at_emergence = round(pdb[k].get('plant_density_at_emergence',plant_density[k]))
                density_at_harvest = round(1. * pdb[k]['ear_density_at_harvest'] / tdb[k]['ears_per_plant'])
                fits[k] = {'HS':[0] + thining[k] + [20], 'density': [density_at_emergence] * 2 + [density_at_harvest] * 2}
    
    density_fits = {k:{'sowing_density': pdb[k]['sowing_density'], 
                       'inter_row': pdb[k]['inter_row'],
                       'density_table': pandas.DataFrame(fits[k])} for k in fits}

                       
    adjust = parameters.get('adjust_density')
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

                
    return density_fits


#
if run_plots:
    parameters = reconstruction_parameters()
    fits = density_fits(**parameters)
    obs = archidb.validation_data()
    archi_plot.density_plot(obs.PlantDensity, fits, HS_fit())
                    
                

#
# Fit tillering
#****************
#
# fitting functions
def _tsurvival(tiller_damages=None, HS_converter=None):
    survivals = {k:None for k in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']}
    if tiller_damages is not None:
        if HS_converter is None:
            HS_converter = fit_HS()
        conv = HS_converter
        for k in survivals:
            damages = tiller_damages[k]
            if damages is not None:
                tillers = damages.pop('tillers')
                df = pandas.DataFrame(damages)
                df['TT'] = conv[k].TT(df['HS'])
                survivals[k] = {T:df for T in tillers} 
    return survivals
#
def _tfit(tdb, delta_stop_del, n_elongated_internode, max_order, tiller_survival):
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    options={}
    return WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff,delta_stop_del=delta_stop_del,n_elongated_internode = n_elongated_internode, max_order=max_order, tiller_survival=tiller_survival)

#
# Fits
def tillering_fits(**parameters):


    # Emission probabilities
    #
    tdb = archidb.Tillering_data()
    
    # apply manual reductions
    emf = parameters.get('emission_reduction')
    if emf is not None:
        for k in emf:
            if emf[k] is not None:
                for T in emf[k]:
                    tdb[k]['emission_probabilities'][T] *= emf[k][T]
    
    # Tiller damages
    HS_converter = parameters.get('HS_converter')
    tiller_damages = parameters.get('tiller_damages')
    tiller_survival = _tsurvival(tiller_damages=tiller_damages, HS_converter=HS_converter)
    
    # Wheat tillering model parameters
    delta_stop_del = parameters.get('delta_stop_del')
    n_elongated_internode = parameters.get('n_elongated_internode') 
    if n_elongated_internode is None: # to do :use wheat Tillering decimal internode number and dimension fits to compute it as did plantgen
        raise NotImplementedError("Not yet implemented")
    #
    max_order = parameters.get('max_order')
    
    t_fits={k: _tfit(tdb[k], delta_stop_del=delta_stop_del[k],n_elongated_internode=n_elongated_internode[k], max_order=max_order, tiller_survival=tiller_survival[k]) for k in tdb}

    return t_fits


if run_plots:
    # populate with fits
    pars = reconstruction_parameters()
    HS_converter = fit_HS()
    fits = tillering_fits(HS_converter=HS_converter, **pars)
    obs = archidb.tillers_per_plant()
    
    #1 graph tillering par var -> 4 graph sur une feuille
    archi_plot.multi_plot_tillering(obs, fits, HS_converter, pars['delta_stop_del'])  
    
    #tallage primary pour les 4 var   
    archi_plot.tillering_primary(obs, fits, HS_converter, pars['delta_stop_del'])
    
    #tallage total pour les 4 var    
    archi_plot.tillering_tot(obs, fits, HS_converter, delta_stop_del)
   
   # primary emissions
    archi_plot.graph_primary_emission(archidb)
   
  
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
    

class EchapReconstructions(object):
    
    def __init__(self, reset_data=False,  median_leaf=True, top_leaves={'Mercia':3, 'Rht3':2, 'Tremie12':4, 'Tremie13':3}):
        self.pars = reconstruction_parameters()
        
        self.pgen_base = self.pars['pgen_base']
        converter = HS_fit(reset=True)
        self.density_fits = density_fits(converter, reset_data=reset_data, **self.pars)
        self.tillering_fits = tillering_fits(HS_converter=converter, **self.pars)
        self.dimension_fits = archidb.dimension_fits()
        self.HS_GL_fits = HS_GL_fits(converter)
        if median_leaf:
            self.leaf_fits = median_leaf_fits(top_leaves=top_leaves)
        else:
            self.leaf_fits = leaf_fits()

    def get_pars(self, name='Mercia', nplants=1, dimension=1):
        """ 
        Construct devT tables from models, considering one reconstruction per nff (ie composite)
        """
        
        conv = self.HS_GL_fits[name]['HS']
        density = self.density_fits[name]['density_at_emergence']
                
        Tillering = self.tillering_fits[name]
        pgens = Tillering.to_pgen(nplants, density = density, phyllochron = conv.phyllochron(), TTem = conv.TT_hs_0, pgen_base = self.pgen_base)
        pars = {}
        
        for k in pgens:
            dimT = deepdd(self.dimension_fits[name][k]).copy()
            dimT['W_blade'] = dimT['W_blade']*(math.sqrt(dimension))
            dimT['L_blade'] = dimT['L_blade']*(math.sqrt(dimension))
            
            nff = k
            glfit = self.HS_GL_fits[name]['GL']
            
            TT_t1 = conv.TT(glfit.hs_t1)
            TT_t2 = conv.TT(glfit.hs_t2)
            a_nff = nff * 1.0 / (TT_t2 - conv.TT_hs_0)
            dynpars = {'a_cohort': a_nff, 
                       'TT_hs_0': conv.TT_hs_0,
                       'TT_hs_N_phytomer_potential': TT_t2,
                       'n0': glfit.n0,
                       'n1': glfit.n1,
                       'n2': glfit.n2}
                                   
            dynT = pgen_ext.dynT_user(dynpars, Tillering.primary_tiller_probabilities.keys())
            
            GL = glfit.GL_number()
            GL['TT'] = conv.TT(GL['HS'])
            GL = dict(zip(GL['TT'],GL['GL']))
            
            pgens[k].update({'dimT_user':dimT, 'dynT_user':dynT, 'GL_number':GL, 'TT_t1_user':TT_t1})
            
                
            pars[k] = pgen_ext.pgen_tables(pgens[k])
            axeT, dimT, phenT = pars[k]['adelT']
            tx = k * 10000
            axeT['id_dim'] = axeT['id_dim']+tx; dimT['id_dim'] = dimT['id_dim']+tx
            axeT['id_phen'] = axeT['id_phen']+tx; phenT['id_phen'] = phenT['id_phen']+tx
        
        return pars
   
    def get_reconstruction(self, name='Mercia', nplants=30, nsect=3, seed=1, sample='sequence', disc_level=7, aborting_tiller_reduction=1, aspect = 'square', stand_density_factor = {'Mercia':1, 'Rht3':1, 'Tremie12':1, 'Tremie13':1}, dimension=1, ssipars={'r1':0.07,'ndelsen':3},**kwds):
        '''stand_density_factor = {'Mercia':0.9, 'Rht3':1, 'Tremie12':0.8, 'Tremie13':0.8}, **kwds)'''
        
        
        run_adel_pars = {'rate_inclination_tiller': 15, 'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1}
        if name == 'Rht3':
            incT=75
            dep=10
        else:
            incT=60
            dep=7
            
        d = self.density_fits[name]
        stand = AgronomicStand(sowing_density=d['sowing_density'], plant_density=d['density_at_emergence'], inter_row=d['inter_row'], noise=0.04)       
        n_emerged, domain, positions, area = stand.stand(nplants, aspect)
               
        pars = self.get_pars(name=name, nplants=n_emerged, dimension = dimension)
        axeT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][0] for k in pars])
        dimT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][1] for k in pars])
        phenT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][2] for k in pars])
        axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])
        devT = devCsv(axeT, dimT, phenT)
        
        #adjust density according to density fit
        if d['density_at_harvest']  < d['density_at_emergence'] and nplants > 1:
            TT_stop_del = self.tillering_fits[name].delta_stop_del / self.HS_GL_fits[name]['HS'].a_cohort
            #to do test for 4 lines table only
            dt = d['density_table'].iloc[1:-1]# remove flat part to avoid ambiguity in time_of_death
            devT = pgen_ext.adjust_density(devT, dt, TT_stop_del = TT_stop_del)
            
        # adjust tiller survival if early death occured
        if self.tillering_fits[name].tiller_survival is not None:
            TT_stop_del = self.tillering_fits[name].delta_stop_del / self.HS_GL_fits[name]['HS'].a_cohort
            cohort_survival = self.tillering_fits[name].cohort_survival()
            devT = pgen_ext.adjust_tiller_survival(devT, cohort_survival, TT_stop_del = TT_stop_del)
            
        leaves = self.leaf_fits[name] 
        
        return AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, aborting_tiller_reduction = aborting_tiller_reduction, aspect = aspect,incT=incT, dep=dep, run_adel_pars = run_adel_pars, **kwds)

def save_EchapReconstructions():

    filename = str(shared_data(alinea.echap)/'EchapReconstructions.pckl')
    f = open(filename, 'w')
    pickle.dump(EchapReconstructions(), f)
    f.close()
    
def get_EchapReconstructions():
 
    filename = str(shared_data(alinea.echap)/'EchapReconstructions.pckl')
    f = open(filename)
    reconst = pickle.load(f)
    f.close()
    return reconst
    
    
# checks consistency adel/fits

if run_plots:
    e=EchapReconstructions()
    adel = e.get_reconstruction()
    df = adel.checkAxeDyn()
    fit = e.density_fits['Mercia']['density_table']
    ax=fit.plot('TT','density',style='-')
    df.plot('TT','nbplants',style='--',ax=ax)
    