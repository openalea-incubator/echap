""" Reconstruction of Wheat for Boigneville data

Current reconstructions use fit from Mariem (pgen data) + new modification consigned here
"""
import pandas
import numpy
from copy import deepcopy

from alinea.adel.astk_interface import AdelWheat
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.WheatTillering import WheatTillering

import alinea.echap.architectural_data as archidb

import alinea.adel.plantgen_extensions as pgen_ext

import matplotlib.pyplot as plt
plt.ion()

#generic function to be moved to adel
  
def plantgen_to_devT_comp(pgen_pars):
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen_pars)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    return axeT, dimT, phenT
    
    
def check_nff(devT):
    """ Count the probailitiy of occurence of MS with given nff
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
    return density_fits
        
def density_plot():
    
    density = archidb.PlantDensity()  
    fits = density_fits()
    grouped = density.groupby('Var')
    
    dens_mercia = grouped.get_group('Mercia'); dens_rht3 = grouped.get_group('Rht3')
    dens_tremie12 = grouped.get_group('Tremie12'); dens_tremie13 = grouped.get_group('Tremie13')

    plt.errorbar(dens_mercia['HS'], dens_mercia['density'], yerr=dens_mercia['SD'], fmt='or', label = 'Mercia density')
    plt.errorbar(dens_rht3['HS'], dens_rht3['density'], yerr=dens_rht3['SD'], fmt='og', label = 'Rht3 density')
    plt.errorbar(dens_tremie12['HS'], dens_tremie12['density'], yerr=dens_tremie12['SD'], fmt='ob', label = 'Tremie12 density')
    plt.errorbar(dens_tremie13['HS'], dens_tremie13['density'], yerr=dens_tremie13['SD'], fmt='om', label = 'Tremie13 density')
   
    for g in fits:       
        if g=='Mercia':
            color='r'
        elif g=='Rht3':
            color='g'
        elif g=='Tremie12':
            color='b'
        else:
            color='m'
        fits[g].plot('HS', 'density', style='-'+color, label=g+' density fits')
      
    plt.title("Plant density"); plt.xlabel("HS"); plt.ylim(ymin=0); plt.ylim(ymax=350)
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})

#
# Fit tillering
#
def _tfit(tdb, delta_stop_del):
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    return WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff,delta_stop_del=delta_stop_del)
    
def tillering_fits(delta_stop_del=2.5):
    t_fits={}
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()
    t_fits['Mercia'] = _tfit(tdb['Mercia'], delta_stop_del=delta_stop_del)
    t_fits['Rht3'] = _tfit(tdb['Rht3'], delta_stop_del=delta_stop_del)
    tdb = archidb.Tillering_data_Tremie12_2011_2012()
    t_fits['Tremie12'] = _tfit(tdb, delta_stop_del=delta_stop_del/5.)#effect of gel
    tdb = archidb.Tillering_data_Tremie13_2012_2013()
    pdata = archidb.Plot_data_Tremie_2012_2013()
    tdb['ears_per_plant'] = pdata['ear_density_at_harvest'] / pdata['mean_plant_density']
    t_fits['Tremie13'] = _tfit(tdb, delta_stop_del=delta_stop_del)
    return t_fits
'''
def plot_tillering(name='Mercia', delta_stop_del=2.5):
    fits = tillering_fits(delta_stop_del)
    fit = fits[name].axis_dynamics(include_MS = False)
    fit.plot('HS', 'total', style='--r', label='Total '+name)
    fit.plot('HS', 'primary', style='-b', label='Primary '+name)
    fit.plot('HS', 'others', style=':g', label= 'Others '+name)

    obs = archidb.tillers_per_plant()
    grouped = obs.groupby('Var')
    obs = grouped.get_group(name)
    obs.plot('HS', ['TP', 'TS', 'TPS','TT3F','FT'],style=['pb','pg','pr','py','pr'])
'''    
def multi_plot_tillering(delta_stop_del=2.5):
    fig, axes = plt.subplots(nrows=2, ncols=2)
    ax0, ax1, ax2, ax3 = axes.flat
    
    varieties = [['Mercia',ax0],['Rht3',ax1],['Tremie12',ax2],['Tremie13',ax3]]
    
    for name,ax in varieties :

        fits = tillering_fits(delta_stop_del)
        fit = fits[name].axis_dynamics(include_MS = False)
        ax.plot(fit['HS'], fit['total'], '--r', label='Total')
        ax.plot(fit['HS'], fit['primary'], '-b', label='Primary')
        ax.plot(fit['HS'], fit['others'], ':g', label= 'Others')

        obs = archidb.tillers_per_plant()
        grouped = obs.groupby('Var')
        obs = grouped.get_group(name)
        ax.plot(obs['HS'], obs['TP'], 'pb', label='TP')
        ax.plot(obs['HS'], obs['TS'], 'pg', label='TS')
        ax.plot(obs['HS'], obs['TPS'], 'pr', label='TPS')
        ax.plot(obs['HS'], obs['TT3F'], 'py', label='TT3F')
        ax.plot(obs['HS'], obs['FT'], 'pm', label='FT')
        
        if name is 'Mercia':
            ax.set_xlim([0, 18]); ax.set_ylim([0, 10]) 
        elif name is 'Rht3':
            ax.set_xlim([0, 18]); ax.set_ylim([0, 10]) 
        else :
            ax.set_xlim([0, 18]); ax.set_ylim([0, 7]) 
            
        ax.set_title(name, fontsize=10)
    
    ax1.legend(numpoints=1, bbox_to_anchor=(1.2, 1.2), prop={'size': 9})
    fig.suptitle("Tillering")
   
reconst_db={}

dynT_MS = {'a_cohort':0.009380186,
            'TT_col_0':101.4740799,
            'TT_col_N_phytomer_potential':1380.766379, # Will be recomputed as a function of Nff
            'n0':4.7,'n1':2.5,'n2':5}

def echap_reconstruction(nplants, density, Tillering, dimensions, dynamic, green_leaves, dTT_stop=0):
    """ Construct devT tables from models, considering one reconstruction per nff (ie composite)
    """
    from alinea.adel.AdelR import devCsv
    
    pgens = Tillering.to_pgen(nplants, density)
    adelT = {}
    pars = {}
    
    for k in pgens:
        dimT = dimensions[k]
        dynamic['TT_col_N_phytomer_potential'] = dynamic['TT_col_0'] + k * 1. / dynamic['a_cohort']
        dynT = pgen_ext.dynT_user(dynamic, Tillering.primary_tiller_probabilities.keys())
        GL = green_leaves.ix[:,['TT',str(k)]]
        GL = GL.ix[GL['TT'] > dynamic['TT_col_N_phytomer_potential'],:]
        GL = dict(zip(GL['TT'],GL[str(k)]))
        pgens[k].update({'dimT_user':dimT, 'dynT_user':dynT, 'GL_number':GL})
        # complete plantgen default
        dtt = 600
        dtt -= dTT_stop
        pgens[k].update({'delais_TT_stop_del_axis':dtt,'TT_col_break': 0.0,})
        adelT[k], pars[k] = pgen_ext.build_tables(pgens[k])
        axeT, dimT, phenT = adelT[k]
        tx = k * 10000
        axeT['id_dim'] = axeT['id_dim']+tx; dimT['id_dim'] = dimT['id_dim']+tx
        axeT['id_phen'] = axeT['id_phen']+tx; phenT['id_phen'] = phenT['id_phen']+tx
    
    axeT = reduce(lambda x,y : pandas.concat([x,y]),[adelT[k][0] for k in adelT])
    dimT = reduce(lambda x,y : pandas.concat([x,y]),[adelT[k][1] for k in adelT])
    phenT = reduce(lambda x,y : pandas.concat([x,y]),[adelT[k][2] for k in adelT])
    
    return devCsv(axeT, dimT, phenT), pars
                
            
# NORMAL : composite with raw observations
def Mercia_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7, **kwds):

    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Mercia']
    
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    plant_density_at_harvest = float(pdata['ear_density_at_harvest']) / ears_per_plant
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    
    wfit = WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff)
    
    dimT = archidb.Mercia_2010_fitted_dimensions()
    GL = archidb.GL_number()['Mercia']
    
    devT, pars = echap_reconstruction(nplants,plant_density_at_harvest,wfit, dimT, dynT_MS, GL, dTT_stop = dTT_stop)

    pars.update({'tillering':wfit})
        
    xy, sr, bins = archidb.leaf_curvature_data('Mercia')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
        
    stand = AgronomicStand(sowing_density=pdata['sowing_density'], plant_density=pdata['plant_density_at_emergence'],inter_row=pdata['inter_row'])

    adel = AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, **kwds)
    
    return pars, adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants
    
reconst_db['Mercia'] = Mercia_2010

def Rht3_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7, face_up=False, classic=False):

    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Rht3']
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Rht3']
    
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    plant_density_at_harvest = float(pdata['ear_density_at_harvest']) / ears_per_plant
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    
    wfit = WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff)
    dimT = archidb.Rht3_2010_fitted_dimensions()
    GL = archidb.GL_number()['Rht3']
    
    devT, pars = echap_reconstruction(nplants,plant_density_at_harvest,wfit, dimT, dynT_MS, GL, dTT_stop = dTT_stop)
    
    pars.update({'tillering':wfit})
       
    xy, sr, bins = archidb.leaf_curvature_data('Rht3')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)

    #stand = AgronomicStand(sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row'])
    stand = AgronomicStand(sowing_density=pdata['sowing_density'], plant_density=pdata['plant_density_at_emergence'],inter_row=pdata['inter_row'])

    adel = AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, incT=22, dinT=22, face_up=face_up)
   
    return pars, adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants
    
reconst_db['Rht3'] = Rht3_2010

def Tremie12(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7):
#BUG: nplants not generated as expected for this genotype !
    pdata = archidb.Plot_data_Tremie_2011_2012()
    tdb = archidb.Tillering_data_Tremie12_2011_2012()
    
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    plant_density_at_harvest = float(pdata['ear_density_at_harvest']) / ears_per_plant
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    
    wfit = WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff)
    
    dimT = archidb.Tremie12_fitted_dimensions()
    GL = archidb.GL_number()['Tremie12']

    devT, pars = echap_reconstruction(nplants,plant_density_at_harvest,wfit, dimT, dynT_MS, GL, dTT_stop = dTT_stop)

    pars.update({'tillering':wfit})    
      
    xy, sr, bins = archidb.leaf_curvature_data('Tremie')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    
    stand = AgronomicStand(sowing_density=pdata['sowing_density'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row'])

    adel = AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand=stand, seed=seed, sample=sample, leaves = leaves)
    
    return pars, adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants
    
reconst_db['Tremie12'] = Tremie12

def Tremie13(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7):

    pdata = archidb.Plot_data_Tremie_2012_2013()
    tdb = archidb.Tillering_data_Tremie13_2012_2013()
    
    ms_nff_probas = tdb['nff_probabilities']
    nff = sum([int(k)*v for k,v in ms_nff_probas.iteritems()]) 
    ears_per_plant = tdb['ears_per_plant']
    plant_density_at_harvest = float(pdata['ear_density_at_harvest']) / ears_per_plant
    primary_emission = {k:v for k,v in tdb['emission_probabilities'].iteritems() if v > 0}
    
    wfit = WheatTillering(primary_tiller_probabilities=primary_emission, ears_per_plant = ears_per_plant, nff=nff)
    
    dimT = archidb.Tremie13_fitted_dimensions()
    GL = archidb.GL_number()['Tremie13']

    devT, pars = echap_reconstruction(nplants,plant_density_at_harvest,wfit, dimT, dynT_MS, GL, dTT_stop = dTT_stop)

    pars.update({'tillering':wfit})  
    
    xy, sr, bins = archidb.leaf_curvature_data('Tremie')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    
    stand = AgronomicStand(sowing_density=pdata['sowing_density'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row'])
    
    adel = AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand=stand , seed=seed, sample=sample, leaves = leaves)
    
    return pars, adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants
    
reconst_db['Tremie13'] = Tremie13

 
age_at_application = {'Mercia_2010': {'2011-04-19': 1166,
                                      '2011-05-11' : 1500}}
                                
                                   