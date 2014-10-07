""" Reconstruction of Wheat for Boigneville data

Current reconstructions use fit from Mariem (pgen data) + new modification consigned here
"""
import pandas
import numpy
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
    dfxy_cut = pandas.cut(dfxy.age, bins)
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
    
    
def leaf_shapes():
    d = {k: {'shapes' : _addLindex(*archidb.leaf_curvature_data(k)),
            'geoLeaf' : geoLeaf()} for k in ['Mercia','Rht3', 'Tremie12', 'Tremie13']}
    return d
    
# Attention pour rht3 on veut geoleaf qui retourne 1 (= feuile du bas) quelque soit le rang !
    # creation de la colonne age et du selecteur d'age 
    
            
# Standard echap reconstruction protocol for plant development and axis dynamic       
#def echap_development(nplants, pgen, primary_proba, tdata, pdata, dTT_stop)
    #ici include mortality


#
# Parametres communs / Pre-traitement des donnees

#
# Fit plant density
#
# use mean plant density and/or density at the end as fits
def density_fits():
    """
    Manual fit of plant density based on mean plant density and estimate of plant density at harvest
    """
    density_fits = {'Mercia':pandas.DataFrame({'HS':[0,6,13,20],                                           'density':[203,203,153,153]}),
                'Rht3': pandas.DataFrame({'HS':[0,6,13,20],'density':[211,211,146,146]}),
                'Tremie12': pandas.DataFrame({'HS':[0,15,20],'density':[281,281,251]}),
                'Tremie13': pandas.DataFrame({'HS':[0,20],'density':[233,233]})}
                
    conv = HS_converter
    
    for k in density_fits:
        df = density_fits[k]
        df['TT'] = conv[k].TT(df['HS'])
    return density_fits
#
if run_plots:
    archi_plot.density_plot(archidb.PlantDensity(), density_fits())
 
#
# Fit tillering
#
def _tfit(tdb, delta_stop_del, n_elongated_internode, max_order):
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    return WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff,delta_stop_del=delta_stop_del,n_elongated_internode = n_elongated_internode, max_order=max_order)
    
def tillering_fits(delta_stop_del=2.8, n_elongated_internode={'Mercia':4, 'Rht3':4, 'Tremie12': 7, 'Tremie13':4} , max_order=None):
        
    tdb = archidb.Tillering_data()
    pdata = archidb.Plot_data_Tremie_2012_2013()
    tdb['Tremie13']['ears_per_plant'] = pdata['ear_density_at_harvest'] / pdata['mean_plant_density']
    
    if not isinstance(delta_stop_del,dict):
        delta_stop_del = {k:delta_stop_del for k in tdb}
        
    if not isinstance(n_elongated_internode,dict):
        n_elongated_internode = {k:n_elongated_internode for k in tdb}   
    
    t_fits={k: _tfit(tdb[k], delta_stop_del=delta_stop_del[k],n_elongated_internode=n_elongated_internode[k], max_order=max_order) for k in tdb}

    return t_fits

if run_plots:
    delta_stop_del = 2.8
    #delta_stop_del={'Mercia':3.5, 'Rht3':2.5, 'Tremie12': 0.5, 'Tremie13':2.5}#handle freezing effect on Tremie12
    fits = tillering_fits(delta_stop_del=delta_stop_del, max_order=None)
    obs = archidb.tillers_per_plant()
    archi_plot.multi_plot_tillering(obs, fits, delta_stop_del)
   
if run_plots:
    archi_plot.graph_primary_emission(archidb)

    
reconst_db={}

#
# Fit dynamique
#
# ---------------------------------------------------- Fitted HS = f(TT)
#
class HaunStage(object):
    """ Handle HaunStage = f (ThermalTime) fits
    """
    
    def __init__(self, a_cohort = 1. / 110., TT_col_0 = 0):
        self.a_cohort = a_cohort
        self.TT_col_0 = TT_col_0
        
    def __call__(self, TT):
        return (numpy.array(TT) - self.TT_col_0) * self.a_cohort
        
    def TT(self, HS):
        return self.TT_col_0 + numpy.array(HS) / self.a_cohort

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
#
HS_converter = {k: HaunStage(a_cohort, TTcol0[k]) for k in TTcol0}
#
if run_plots:
    archi_plot.dynamique_plot(archidb.HS_GL_SSI_data(), converter = HS_converter)
#
def GL_fits():
    # Maxwell values (read on mariem paper) : n0 : 4..4, n1: 3.34, n2: 5.83
    # common fit mercia, rht3, Tremie12
    conv = HS_converter['Tremie12']
    xref = conv(archidb.GL_number()['Tremie12']['TT'])
    yref = archidb.GL_number()['Tremie12']['mediane']
    nff = archidb.mean_nff()
    coefs  = {k:{'n0': 4.4,'n1':1.9,'n2':numpy.interp(nff[k],xref, yref), 'hs_t1': 8, 'hs_t2':nff[k]} for k in nff}
    coefs['Tremie13']=coefs['Tremie13'].update({'n1':2.66, 'n2':4})
    #
    fits = {k:{'HS_fit': HS_converter[k], 'coefs': coefs[k], 'GL': v} for k,v in archidb.GL_number().iteritems()}
    # Mercai, Rht3, use Tremie 12 data instead ?
    return fits
          
def dynamique_fits():

    dynT_MS = {'Others': {'a_cohort':0.009380186,
                'TT_col_0':101.4740799,
                'TT_col_N_phytomer_potential':1380.766379, # Will be recomputed as a function of Nff
                'n0':4.764532295,'n1':2.523833441,'n2':5.083739846}, #Mercia dynT_user
        'Tremie13': {'a_cohort':0.009380186,
                'TT_col_0':101.4740799,
                'TT_col_N_phytomer_potential':1254, # Will be recomputed as a function of Nff
                'n0':4.2,'n1':2.6,'n2':3.6}} # a_cohort et TT_col_0 de Mercia           
    return dynT_MS

def HS_GL_SSI_data_HS():
    df = archidb.HS_GL_SSI_data()
    for g in ('Mercia','Rht3', 'Tremie12', 'Tremie13'):
        HS_GL_SSI_data = archidb._add_ghs(df, g)
    return HS_GL_SSI_data_HS
#
if run_plots:
    archi_plot.dynamique_plot_nff(archidb.HS_GL_SSI_data(), dynamique_fits())   
     
if run_plots:
    archi_plot.dynamique_plot(HS_GL_SSI_data())



    

class EchapReconstructions(object):
    
    def __init__(self):
        self.density_fits = density_fits()
        self.tillering_fits = tillering_fits(delta_stop_del=2.8, n_elongated_internode={'Mercia':4, 'Rht3':4, 'Tremie12':7, 'Tremie13':4} , max_order=None)
        self.dimension_fits = archidb.dimension_fits()
        self.GL_fits = GL_fits()
        self.leaf_shapes = leaf_shapes()
        self.plot_data = {k:{kk:v[kk] for kk in ['inter_row', 'sowing_density']} for k,v in archidb.Plot_data().iteritems()}
    
    def get_pars(self, name='Mercia', nplants=1, density = 1, force_start_reg = False):
        """ 
        Construct devT tables from models, considering one reconstruction per nff (ie composite)
        """
        Tillering = self.tillering_fits[name]
        pgens = Tillering.to_pgen(nplants, density)
        pars = {}
        
        for k in pgens:
            dimT = self.dimension_fits[name][k]
            
            glfit = self.GL_fits[name]
            conv = glfit['HS_fit']
            dynpars = {'a_cohort': conv.a_cohort, 
                       'TT_col_0': conv.TT_col_0,
                       'TT_col_N_phytomer_potential': conv.TT(k)}
            dynpars.update(glfit['coefs'])
            dynT = pgen_ext.dynT_user(dynpars, Tillering.primary_tiller_probabilities.keys())
            
            GL = glfit['GL'].ix[:,['TT',str(k)]]
            GL = GL.ix[GL['TT'] > dynpars['TT_col_N_phytomer_potential'],:]
            GL = dict(zip(GL['TT'],GL[str(k)]))
            
            pgens[k].update({'dimT_user':dimT, 'dynT_user':dynT, 'GL_number':GL})
            
            if 'hs_deb_reg' in pgens[k]:
                hs_deb_reg = pgens[k].pop('hs_deb_reg')
                TT_start_reg = conv.TT(hs_deb_reg)
                pgens[k].update({'TT_start_reg_user': TT_start_reg})

            pars[k] = pgen_ext.pgen_tables(pgens[k])
            axeT, dimT, phenT = pars[k]['adelT']
            tx = k * 10000
            axeT['id_dim'] = axeT['id_dim']+tx; dimT['id_dim'] = dimT['id_dim']+tx
            axeT['id_phen'] = axeT['id_phen']+tx; phenT['id_phen'] = phenT['id_phen']+tx
        
        return pars
   
    def get_reconstruction(self, name='Mercia', nplants=30, nsect=3, seed=1, sample='sequence', disc_level=7, aborting_tiller_reduction=1, aspect = 'square', **kwds):
    
        density = self.density_fits[name]
        density_at_emergence = density['density'][density['HS'] == 0].iloc[0]
        density_at_harvest = density['density'][density['HS'] == max(density['HS'])].iloc[0]
        
        pdata = self.plot_data[name]
        stand = AgronomicStand(sowing_density=pdata['sowing_density'], plant_density=density_at_emergence, inter_row=pdata['inter_row'])        
        n_emerged, domain, positions, area = stand.stand(nplants, aspect)
        
        pars = self.get_pars(name=name, nplants=n_emerged, density = density_at_emergence, force_start_reg=False)#force_start_reg remains to be tested , but setting to True is needed for tremie12
        axeT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][0] for k in pars])
        dimT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][1] for k in pars])
        phenT = reduce(lambda x,y : pandas.concat([x,y]),[pars[k]['adelT'][2] for k in pars])
        axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])
        devT = devCsv(axeT, dimT, phenT)
        #adjust density according to density fit
        if density_at_harvest  < density_at_emergence and nplants > 1:
            devT = pgen_ext.adjust_density(devT, density)
        
        dfxy, dfsr = self.leaf_shapes[name]['shapes']
        xy, sr, bins = leaf_trajectories(dfxy, dfsr , bins = [-10, 0.5, 1, 2, 3, 4, 10], ntraj = 10, tol_med = 0.1)
        geoLeaf = self.leaf_shapes[name]['geoLeaf']
        leaves = Leaves(xy, sr, geoLeaf=geoLeaf, dynamic_bins = bins, discretisation_level = disc_level)
        
        return AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, aborting_tiller_reduction = aborting_tiller_reduction, aspect = aspect, **kwds)
                             