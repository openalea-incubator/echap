""" Reconstruction of Wheat for Boigneville data

Current reconstructions use fit from Mariem (pgen data) + new modification consigned here
"""
import pandas
import numpy
import math
import scipy.stats as stats
from copy import deepcopy

from alinea.adel.astk_interface import AdelWheat
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.WheatTillering import WheatTillering
from alinea.adel.AdelR import devCsv

import alinea.echap.architectural_data as archidb
import alinea.echap.architectural_reconstructions_plot as archi_plot

import alinea.adel.plantgen_extensions as pgen_ext

run_plots = False # prevent ipython %run to make plots

#generic function to be moved to adel
  
def plantgen_to_devT_comp(pgen_pars):
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen_pars)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    return axeT, dimT, phenT
    
    
def check_nff(devT):
    """ Count the probability of occurence of MS with given nff
    """
    nffs = devT['axeT']['N_phytomer'][devT['axeT']['id_axis']=='MS']
    counts = {n:nffs.tolist().count(n) for n in set(nffs)}
    probas = {n:nffs.tolist().count(n) * 1.0 / len(nffs) for n in set(nffs)}
    return counts, probas
    
    
def check_primary_tillers(devT):
    """ Count/estimate probabilitie of occurence of primary tillers
    """
    import re 
    axis = devT['axeT']['id_axis']
    ms = [e for e in axis if re.match('MS',e)]
    tillers = [e for e in axis if re.match('T.$',e)]
    counts = {n:tillers.count(n) for n in set(tillers)}
    probas = {n:tillers.count(n) * 1.0 / len(ms) for n in set(tillers)}
    return counts, probas
    
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

def median_leaf_fits(disc_level=7):
    d={}
    gL = geoLeaf()
    trajs,bins = archidb.median_leaf_trajectories()
    for k in ['Mercia','Rht3', 'Tremie12', 'Tremie13']:
        _, dfsr = _addLindex(*archidb.leaf_curvature_data(k))
        srdb = {k:v.ix[:,['s','r']].to_dict('list') for k, v in dfsr.groupby('Lindex')}
        d[k] = Leaves(trajs, srdb, geoLeaf=gL, dynamic_bins = bins, discretisation_level = disc_level)
    return d
 
# Attention pour rht3 on veut geoleaf qui retourne 1 (= feuile du bas) quelque soit le rang !
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
#
HS_converter = fit_HS()

if run_plots:
    archi_plot.dynamique_plot(archidb.HS_GL_SSI_data(), converter = HS_converter) # weighted (frequency of nff modalities) mean 

#
def HS_GL_fits():
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
    coefs = {'n0': 4.4, 'hs_t1': 8.1}#n0=4.4 : Maxwell value
    n1 = {'Mercia':1.9, 'Rht3': 1.9, 'Tremie12':1.9, 'Tremie13':2.8}
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
    for name in ('Mercia','Rht3','Tremie12','Tremie13'):
        e = EchapReconstructions()
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
    archi_plot.dynamique_plot_sim(archidb.HS_GL_SSI_data(), pars(), converter = HS_converter)
 
#
# Fit plant density
#
# use mean plant density and/or density at the end as fits
class deepdd(pandas.DataFrame):

    def deepcopy(self):
        return self.copy(deep=True)
#

def density_fits():
    """
    Manual fit of plant density based on mean plant density and estimate of plant density at harvest
    """
    density_fits = {'Mercia':deepdd({'HS':[0,6,13,20],'density':[203,203,153,153]}),
                'Rht3': deepdd({'HS':[0,6,13,20],'density':[211,211,146,146]}),
                'Tremie12': deepdd({'HS':[0,6,13,20],'density':[281,281,281,281]}),
                'Tremie13': deepdd({'HS':[0,20],'density':[251,251]})} #ne prend pas en compte la densite releve a epis1cm
                
    conv = HS_converter
    
    for k in density_fits:
        df = density_fits[k]
        df['TT'] = conv[k].TT(df['HS'])
    return density_fits


#densite non ajustee 
'''def density_fits():
    """
    Manual fit of plant density based on mean plant density and estimate of plant density at harvest
    """
    density_fits = {'Mercia':deepdd({'HS':[0,6,13,20],'density':[153,153,153,153]}),
                'Rht3': deepdd({'HS':[0,6,13,20],'density':[146,146,146,146]}),
                'Tremie12': deepdd({'HS':[0,6,13,20],'density':[251,251,251,251]}),
                'Tremie13': deepdd({'HS':[0,20],'density':[251,251]})} #ne prend pas en compte la densite releve a epis1cm!
                
    conv = HS_converter
    
    for k in density_fits:
        df = density_fits[k]
        df['TT'] = conv[k].TT(df['HS'])
    return density_fits'''
#
if run_plots:
    archi_plot.density_plot(archidb.PlantDensity(), density_fits(), HS_converter)
#
# Tiller survival to handl fly effects
#    
def tiller_survival():
    conv = HS_converter
    sfly = pandas.DataFrame({'HS':[4,9], 'density':[1,0.7]})# 0.7 is mean survival among living plant
    sfly_Mercia = sfly.copy(deep=True)
    sfly_Mercia['TT'] = conv['Mercia'].TT(sfly['HS'])
    sfly_Rht3 = sfly.copy(deep=True)
    sfly_Rht3['TT'] = conv['Rht3'].TT(sfly['HS'])
    sfly_Rht3['density'][1] = 0.6
    sfreeze_Tremie12 = pandas.DataFrame({'HS':[4.9,5.1], 'density':[1,0.2]})
    sfreeze_Tremie12['TT'] = conv['Tremie12'].TT(sfreeze_Tremie12['HS'])
    survivals = {'Mercia': {T:sfly_Mercia for T in ('T1','T2','T3')},
                  'Rht3':  {T:sfly_Rht3 for T in ('T1','T2','T3')},
                  'Tremie12': {'T3':sfreeze_Tremie12},
                  'Tremie13': None}

    return survivals
                    
fly_damage_tillering = tiller_survival()                

#
# Fit tillering
#
def _tfit(tdb, delta_stop_del, n_elongated_internode, max_order, tiller_survival):
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    return WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff,delta_stop_del=delta_stop_del,n_elongated_internode = n_elongated_internode, max_order=max_order, tiller_survival=tiller_survival)
    
def tillering_fits(delta_stop_del=2.8, n_elongated_internode={'Mercia':3.5, 'Rht3':3.5, 'Tremie12': 5, 'Tremie13':4} , max_order=None, tiller_survival=fly_damage_tillering):
        
    tdb = archidb.Tillering_data()
    pdata = archidb.Plot_data_Tremie_2012_2013()
    #
    tdb['Tremie13']['ears_per_plant'] = pdata['ear_density_at_harvest'] / pdata['mean_plant_density']
    for T in ('T3','T4'):
        tdb['Tremie13']['emission_probabilities'][T] *= 0.5
    
    if not isinstance(delta_stop_del,dict):
        delta_stop_del = {k:delta_stop_del for k in tdb}
        
    if not isinstance(n_elongated_internode,dict):
        n_elongated_internode = {k:n_elongated_internode for k in tdb}   
    
    t_fits={k: _tfit(tdb[k], delta_stop_del=delta_stop_del[k],n_elongated_internode=n_elongated_internode[k], max_order=max_order, tiller_survival=tiller_survival[k]) for k in tdb}

    return t_fits

if run_plots:
    delta_stop_del = 2.
    #delta_stop_del={'Mercia':3.5, 'Rht3':2.5, 'Tremie12': 0.5, 'Tremie13':2.5}#handle freezing effect on Tremie12
    fits = tillering_fits(delta_stop_del=delta_stop_del, n_elongated_internode={'Mercia':3.5, 'Rht3':3, 'Tremie12': 5, 'Tremie13':4}, max_order=None, tiller_survival=fly_damage_tillering)
    obs = archidb.tillers_per_plant()
    archi_plot.multi_plot_tillering(obs, fits, HS_converter, delta_stop_del)
   
if run_plots:
    archi_plot.graph_primary_emission(archidb)

    
reconst_db={}

    
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
    
    def __init__(self, median_leaf=True):
        self.density_fits = density_fits()
        self.tillering_fits = tillering_fits()
        self.dimension_fits = archidb.dimension_fits()
        self.HS_GL_fits = HS_GL_fits()
        if median_leaf:
            self.leaf_fits = median_leaf_fits()
        else:
            self.leaf_fits = leaf_fits()
        self.plot_data = {k:{kk:v[kk] for kk in ['inter_row', 'sowing_density']} for k,v in archidb.Plot_data().iteritems()}
    
    def get_pars(self, name='Mercia', nplants=1, density = 1, force_start_reg = True, dimension=1):
        """ 
        Construct devT tables from models, considering one reconstruction per nff (ie composite)
        """
        Tillering = self.tillering_fits[name]
        pgens = Tillering.to_pgen(nplants, density, force_start_reg=force_start_reg)
        pars = {}
        
        for k in pgens:
            dimT = deepdd(self.dimension_fits[name][k]).copy()
            dimT['W_blade'] = dimT['W_blade']*(math.sqrt(dimension))
            dimT['L_blade'] = dimT['L_blade']*(math.sqrt(dimension))
            
            nff = k
            glfit = self.HS_GL_fits[name]['GL']
            conv = self.HS_GL_fits[name]['HS']
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
            
            if 'hs_deb_reg' in pgens[k]:
                hs_deb_reg = pgens[k].pop('hs_deb_reg')
                TT_start_reg = conv.TT(hs_deb_reg)
                pgens[k].update({'TT_regression_start_user': TT_start_reg})
                
            pars[k] = pgen_ext.pgen_tables(pgens[k])
            axeT, dimT, phenT = pars[k]['adelT']
            tx = k * 10000
            axeT['id_dim'] = axeT['id_dim']+tx; dimT['id_dim'] = dimT['id_dim']+tx
            axeT['id_phen'] = axeT['id_phen']+tx; phenT['id_phen'] = phenT['id_phen']+tx
        
        return pars
   
    def get_reconstruction(self, name='Mercia', nplants=30, nsect=3, seed=1, sample='sequence', disc_level=7, aborting_tiller_reduction=1, aspect = 'square', adjust_density = {'Mercia':0.7, 'Rht3':0.7, 'Tremie12': 0.7, 'Tremie13':None}, dec_density={'Mercia':0, 'Rht3':0, 'Tremie12': 0, 'Tremie13':None}, freeze_damage ={'Mercia':{'T4':0.01,'T5':0.01,'T6':0.01}, 'Rht3':{'T4':0.01,'T5':0.01}, 'Tremie12': None, 'Tremie13':{'T4':0.25}}, stand_density_factor = {'Mercia':1, 'Rht3':1, 'Tremie12':1, 'Tremie13':1}, dimension=1, **kwds):
        '''stand_density_factor = {'Mercia':0.9, 'Rht3':1, 'Tremie12':0.8, 'Tremie13':0.8}, **kwds)'''
    
        density = self.density_fits[name].deepcopy()
        
        density_at_emergence = density['density'][density['HS'] == 0].iloc[0] * stand_density_factor[name]
        density_at_harvest = density['density'][density['HS'] == max(density['HS'])].iloc[0] * stand_density_factor[name]
        
        pdata = self.plot_data[name]
        sowing_density = pdata['sowing_density'] * stand_density_factor[name]
        inter_row = pdata['inter_row']/math.sqrt(stand_density_factor[name])
        stand = AgronomicStand(sowing_density=sowing_density, plant_density=density_at_emergence, inter_row=inter_row)       
        n_emerged, domain, positions, area = stand.stand(nplants, aspect)
        
        if freeze_damage[name] is not None:
            self.tillering_fits = tillering_fits()
            for k in freeze_damage[name]:
                self.tillering_fits[name].primary_tiller_probabilities[k] = freeze_damage[name][k]
        
        pars = self.get_pars(name=name, nplants=n_emerged, density = density_at_emergence, dimension = dimension)
        axeT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][0] for k in pars])
        dimT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][1] for k in pars])
        phenT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][2] for k in pars])
        axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])
        devT = devCsv(axeT, dimT, phenT)
        
        #adjust density according to density fit
        if density_at_harvest  < density_at_emergence and nplants > 1:
            #to do test for 4 lines table only
            d = density.iloc[1:-1]# remove flat part to avoid ambiguity in time_of_death
            devT = pgen_ext.adjust_density(devT, d)
            if adjust_density[name] is not None:
                d['density'].iloc[0] = 1
                d['density'].iloc[1] = adjust_density[name]
                conv = self.HS_GL_fits[name]['HS']
                d['HS'] -= dec_density[name]
                d['TT'] = conv.TT(d['HS'])
                devT = pgen_ext.adjust_tiller_survival(devT, d)
            
        leaves = self.leaf_fits[name] 
        
        return AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, aborting_tiller_reduction = aborting_tiller_reduction, aspect = aspect, **kwds)
                             