import pandas
import numpy
import os
import matplotlib.pyplot as plt
import math
import csv

from alinea.adel.WheatTillering import WheatTillering
import alinea.echap.architectural_data as archidb

from alinea.echap.architectural_reconstructions import reconst_db

from multiprocessing import Pool,freeze_support
import itertools

import openalea.plantgl.all as pgl
from alinea.adel.mtg_interpreter import *
from alinea.echap.weather_data import *
from alinea.caribu.caribu_star import diffuse_source, run_caribu

plt.ion()

def get_reconstruction(name='Mercia', **args):
    fun = reconst_db[name]
    pars, adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants = fun(**args)
    return pars, adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants
'''
def get_pgen(name='Mercia', original = False, dTT_stop = 0):
    fun = reconst_db[name]
    wfit, _, _, _, _, _, _, _ = fun(nplants=1,nsect=1,as_pgen=original, dTT_stop=dTT_stop)
    return wfit
'''
#-------------------------------------------------------------------------------------   
# test new tillering parametrisation against original one by Mariem and data
#
# usage (ipython):
# %pylab
# %run TestCalibration.py
# test_axis_dynamics('Mercia')
#
# Very strange simulation for Tremie
#    
    
def test_axis_dynamics(name='Mercia', name_obs='Mercia', color='r'):

    pars, _, _, _, _, _ = get_reconstruction(name)
    
    wfit = pars['tillering_model']; pgen = pars['config']
    
    if name_obs is "Mercia":
        obs = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']
    elif name_obs is 'Rht3':
        obs = archidb.Plot_data_Mercia_Rht3_2010_2011()['Rht3']
    elif name_obs is 'Tremie12':
        obs = archidb.Plot_data_Tremie_2011_2012()
    else:
        obs = archidb.Plot_data_Tremie_2012_2013()

    '''
    def _getfit(pgen):
        primary_proba = pgen['decide_child_axis_probabilities']
        plant_density = pgen['plants_density']
        ears_density = pgen['ears_density']
        ears_per_plant = float(ears_density) / plant_density
        nff = sum([int(k)*v for k,v in pgen['MS_leaves_number_probabilities'].iteritems()]) #moyenne ponderee
        m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
        fit = m.axis_dynamics(plant_density = plant_density)
        return fit
    
    #fit by Mariem
    #fit = _getfit(pgen)
    #newfit
    new_fit = _getfit(newpgen)
    '''
    plant_density = pgen['plants_density']
    fit = wfit.axis_dynamics(plant_density = plant_density)
    fit.plot('HS', 'total', style='--'+color, label='Total '+name)
    fit.plot('HS', 'primary', style='-'+color, label='Primary '+name)
    fit.plot('HS', 'others', style=':'+color, label= 'Others '+name)

    if 'plant_density_at_emergence' in obs:
        plt.plot([1,13],[obs['plant_density_at_emergence'],obs['ear_density_at_harvest']], 'o'+color) 
    else :
        plt.plot([1,13],[obs['sowing_density'],obs['ear_density_at_harvest']], 'o'+color) 
        
    plt.xlabel("HS")
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})

    return fit

#------------------------------------------------------------------------------------- 
# LAI
    
def simLAI(adel, domain_area, convUnit, nplants):
    from alinea.adel.postprocessing import axis_statistics, plot_statistics 
     
    dd = range(500,2500,300)
    
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res

def compare_LAI(name='Mercia', name_obs='Mercia', dTT_stop=0, original=False, n=30, color='r', **kwds): #name_obs = Mercia/Rht3/Tremie12/Tremie13

    from math import *
    import numpy as np

    pars, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original, **kwds)
    pgen = pars[12]['config']
    sim = simLAI(adel, domain_area, convUnit, nplants)
    sim['HS'] = (sim.ThermalTime - pgen['dynT_user'].TT_col_0[0]) * pgen['dynT_user'].a_cohort[0]
    label = 'LAI vert simule '+name_obs
    sim.plot('HS','LAI_vert',color=color, label=label)
    
    #donnees obs
    obs = archidb.PAI_data()[name_obs]
    obs['HS'] = (obs.TT_date - pgen['dynT_user'].TT_col_0[0]) * pgen['dynT_user'].a_cohort[0]
    #photos
    x=obs['HS']; y=obs['PAI_vert_average_photo']; err=obs['PAI_vert_SD_photo']
    label = 'PAI vert SD photo '+name_obs
    plt.errorbar(x,y,yerr=err, fmt='--o'+color, label=label)
    #calcul intervalle de confiance
    s  = obs['PAI_vert_plt_photo']
    obs['IC'] = (1.96*obs['PAI_vert_SD_photo'])/np.sqrt(s) ; erric = obs['IC']
    #plt.errorbar(x,y,yerr=erric, fmt='--om', label='PAI vert IC photo')
    #biomasse
    if 'LAI_vert_SD_biomasse' in obs:
        x=obs['HS']; y=obs['LAI_vert_average_biomasse']; err=obs['LAI_vert_SD_biomasse']
        label = 'LAI vert SD biomasse '+name_obs
        plt.errorbar(x,y,yerr=err, fmt='--oc', label=label)
    
    plt.xlabel("HS")
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
    return sim
    
#------------------------------------------------------------------------------------- 
# Hauteur de couvert simule

def height(name='Mercia', dTT_stop=0, original=False, n=30):
    from alinea.astk.plantgl_utils import get_height

    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, dTT_stop=dTT_stop, as_pgen=original, nplants=n)
    dd = range(400,3000,100); max_h = []
    sim = [adel.setup_canopy(age) for age in dd]
    
    for g in sim:
        #scene = plot3d(g)
        #pgl.Viewer.display(scene)
        scene_geom = g.property('geometry')
        heights = get_height(scene_geom)
        max_height = max(map(numpy.min,heights.values()))
        max_h.append(max_height)

    h = pandas.DataFrame({'tt':dd, 'height':max_h})
    h['HS'] = (h.tt - pgen['dynT_user'].TT_col_0[0]) * pgen['dynT_user'].a_cohort[0]
    
    return h
        
#------------------------------------------------------------------------------------- 
# Taux de couverture

def draft_TC(g, adel, domain, zenith, rep):
    from alinea.adel.postprocessing import ground_cover
    
    #modelisation afin de voir si erreur
    #scene=plot3d(g)
    #pgl.Viewer.display(scene)
    
    echap_top_camera =  {'type':'perspective', 'distance':200., 'fov':50., 'azimuth':0, 'zenith':zenith}
    #high resolution
    #gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = 4288, image_height = 2848, getImages=True, replicate=rep)
    #resolution /4
    gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = 1072, image_height = 712, getImages=True, replicate=rep)
    
    return gc
    
def comp_TC(name='Mercia_maq1', name_obs='Mercia', original=False, n=30, zenith=0, dTT_stop=0): #zenith = 0 or 57

    if zenith==0:
        zen='0'; rep=1
    else:
        zen='57';rep=2

    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, dTT_stop=dTT_stop, as_pgen=original, nplants=n)  
    dd = range(400,2600,100)
    sim = [adel.setup_canopy(age) for age in dd]
    TC_sim = [draft_TC(g, adel, domain, zenith, rep) for g in sim]
    obs = archidb.TC_data()[name_obs+'_'+zen]

    n=0; tc=[]; tc_sen=[]
    while n<len(TC_sim):
        tc.append(TC_sim[n]['green'])
        tc_sen.append(TC_sim[n]['senescent'])
        n = n + 1

    sim_green = pandas.DataFrame({'tt':dd, 'TCgreen':tc, 'TCsen':tc_sen})
    sim_green['HS'] = (sim_green.tt - pgen['dynT_user'].TT_col_0[0]) * pgen['dynT_user'].a_cohort[0]
    obs['HS'] = (obs.TT - pgen['dynT_user'].TT_col_0[0]) * pgen['dynT_user'].a_cohort[0]
    
    # Tracer TC_tot, TC_green et senescence sur le meme graph
    sim_green['TCtot'] = sim_green['TCgreen'] + sim_green['TCsen']
    sim_green.plot('HS','TCtot',color='b'); sim_green.plot('HS','TCgreen',color='g'); sim_green.plot('HS','TCsen',color='y')
    obs.plot('HS','TC',style='or')
    
    # Tracer les courbes 1-TC et TC_tot pr comparer avec graph rayonnement obs/sim
    #sim_green['1-TC']=1-sim_green['TCgreen']
    #sim_green['TCtot']=sim_green['TCgreen']+sim_green['TC_sen']
    #sim_green['1-TCtot']=1-sim_green['TCtot']
    #sim_green.plot('tt','1-TC',color='y',linestyle='--',linewidth=2)
    #sim_green.plot('tt','1-TCtot',color='y',linewidth=2)
    #obs.plot('TT','1-TC',style='oy')
    # Graph k pr Carole ?
    #sim_green['k_TCgreen'] = numpy.log(1/(1-sim_green['TCgreen']))
    #sim_green['k_TCtot'] = numpy.log(1/(1-sim_green['TCtot']))
    
    return sim_green

#-------------------------------------------------------------------------------------
# GRAPH METEO
    
def mat_ray_obs(data_file):
    import re
    
    header_row = ['DATE','H','RECORD','QR_PAR_Avg','SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4','SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8']
    #data = pandas.read_csv(data_file, parse_dates={'datetime':[0,1]}, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')
    df = pandas.read_csv(data_file, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')
     
    #filtre seulement heure entre 9h et 19h pour chaque journee
    df = df[df['H']>9]; df = df[df['H']<19]
    
    #tableau non traité
    dfa=df; dfa = dfa.reset_index()
    dfa.columns = ['datetime','H','RECORD','QR_PAR_Avg','SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4','SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8']
    dfa = dfa.groupby(dfa['datetime'], axis=0).mean(); dfa = dfa.reset_index()
    da=0
    while da<len(dfa):
        dfa['datetime'][da] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", dfa['datetime'][da])
        da = da+1
    dfa['%SE1'] = dfa['SE_PAR_Avg_1']/dfa['QR_PAR_Avg']; dfa['%SE2'] = dfa['SE_PAR_Avg_2']/dfa['QR_PAR_Avg']
    dfa['%SE3'] = dfa['SE_PAR_Avg_3']/dfa['QR_PAR_Avg']; dfa['%SE4'] = dfa['SE_PAR_Avg_4']/dfa['QR_PAR_Avg']
    dfa['%SE5'] = dfa['SE_PAR_Avg_5']/dfa['QR_PAR_Avg']; dfa['%SE6'] = dfa['SE_PAR_Avg_6']/dfa['QR_PAR_Avg']
    dfa['%SE7'] = dfa['SE_PAR_Avg_7']/dfa['QR_PAR_Avg']; dfa['%SE8'] = dfa['SE_PAR_Avg_8']/dfa['QR_PAR_Avg']
    
    # traitement tableau de données obs
    #avg PAR niveau du sol ('SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4')
    dat0 = df.ix[:,'SE_PAR_Avg_1':'SE_PAR_Avg_4']; dat0 = dat0.apply(numpy.mean,1)
    #avg PAR env 20 cm du sol ('SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8')
    dat20 = df.ix[:,'SE_PAR_Avg_5':'SE_PAR_Avg_8']; dat20 = dat20.apply(numpy.mean,1)
    # tableau % dat0/dat200
    tab_prc0 = dat0/df.QR_PAR_Avg
    # tableau % dat20/dat200
    tab_prc20 = dat20/df.QR_PAR_Avg
    
    # chaine de traitement pour hauteur au niveau du sol
    tab0 = pandas.DataFrame(tab_prc0)
    tab0 = tab0.reset_index(); tab0.columns = ['datetime','%']
    tab0 = tab0.groupby(tab0['datetime'], axis=0).mean()
    tab0 = tab0.reset_index()
    # changement format date JJ/MM/AAAA en JJ-MM-AAAA
    dt=0
    while dt<len(tab0):
        tab0['datetime'][dt] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", tab0['datetime'][dt])
        dt = dt+1

    # chaine de traitement pour hauteur a 20 cm du sol
    tab20 = pandas.DataFrame(tab_prc20)
    tab20 = tab20.reset_index(); tab20.columns = ['datetime','%']
    tab20 = tab20.groupby(tab20['datetime'], axis=0).mean()
    tab20 = tab20.reset_index()
    # changement format date JJ/MM/AAAA en JJ-MM-AAAA
    dte=0
    while dte<len(tab20):
        tab20['datetime'][dte] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", tab20['datetime'][dte])
        dte = dte+1
    
    return dfa, tab0, tab20
    
def draft_light(g, adel, domain, z_level):
    from alinea.caribu.label import Label

    scene = adel.scene(g)
    #modelisation afin de voir si erreur
    #scene=plot3d(g)
    #pgl.Viewer.display(scene)

    sources = diffuse_source(46)
    out = run_caribu(sources, scene, domain=domain, zsoil = z_level)
    labs = {k:Label(v) for k,v in out['label'].iteritems() if not numpy.isnan(float(v))}
    ei_soil = [out['Ei'][k] for k in labs if labs[k].is_soil()]
    return float(numpy.mean(ei_soil))
    
#fonction a lancer pour obtenir graph 'rayonnement obs contre rayonnement sim'
def graph_meteo(name='Mercia', dTT_stop=0, original=False, n=30):
    from alinea.astk.TimeControl import thermal_time # Attention la fonction marche seulement avec freq='H' !!! (pas en jour)
    
    # PARTIE OBS ------------------------------------------------------------------------------------------------------------
    if name is 'Mercia' or 'Rht3': #Mercia & Rht3 2010/2011
        data_file = 'METEO_stationINRA_20102011.csv'
        # correspondance entre date et TT
        met = Boigneville_2010_2011()
        seq = pandas.date_range(start="2010-10-15", end="2011-07-18", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2010-10-15", end="2011-07-18")
        bid = bid[seq]
    if name is 'Tremie': #Tremie 2011/2012
        data_file = 'METEO_stationINRA_20112012.csv'
        # correspondance entre date et TT
        met = Boigneville_2011_2012()
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18")
        bid = bid[seq]
        
    dfa, tab0, tab20 = mat_ray_obs(data_file)
        
    # merger bid avec tab0 puis tab20 
        # mise en forme de bid
    bid = pandas.DataFrame(list(bid.values), index=bid.index)
    bid = bid.reset_index()
    bid.columns = ['datetime','TT']
    bid['datetime'] = bid['datetime'].map(lambda x: x.strftime('%d-%m-%Y'))
        # merge bid avec tab0 puis tab20, classement par ordre croissant
    tab0 = tab0.merge(bid); tab20 = tab20.merge(bid)
    tab0 = tab0.sort(['TT']); tab20 = tab20.sort(['TT'])
    
    # PARTIE SIM ------------------------------------------------------------------------------------------------------------
    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, dTT_stop=dTT_stop, as_pgen=original, nplants=n) 
    tx = pgen['dynT_user'].a_cohort[0]    
    #dd = range(900,2200,100)
    dd = range(400,2600,100)
    sim = [adel.setup_canopy(age) for age in dd]
    #sim=adel.setup_canopy(1380)

    light_sim_0 = [draft_light(g, adel, domain, z_level=0) for g in sim]
    #light_sim_5 = [draft_light(g, adel, domain, z_level=5) for g in sim]
    #light_sim_20 = [draft_light(g, adel, domain, z_level=20) for g in sim]
    #light_sim_25 = [draft_light(g, adel, domain, z_level=25) for g in sim]

    sim0 = pandas.DataFrame({'TT':dd, 'light0':light_sim_0})
    sim0['HS'] = sim0.TT*tx
    #sim5 = pandas.DataFrame({'TT':dd, 'light5':light_sim_5})
    #sim20 = pandas.DataFrame({'TT':dd, 'light20':light_sim_20})
    #sim25 = pandas.DataFrame({'TT':dd, 'light25':light_sim_25})

    # GRAPHES ---------------------------------------------------------------------------------------------------------------
    #obs tous les points
    dfa = dfa.merge(bid); dfa = dfa.sort(['TT']); print dfa.head()
        #niveau du sol (point rouge)
    #dfa.plot('TT','%SE1',style='om',markersize=4); dfa.plot('TT','%SE2',style='om',markersize=4); dfa.plot('TT','%SE3',style='om',markersize=4); dfa.plot('TT','%SE4',style='om',markersize=4)
        #20cm du sol (point bleu)
    #dfa.plot('TT','%SE5',style='oc',markersize=4); dfa.plot('TT','%SE6',style='oc',markersize=4); dfa.plot('TT','%SE7',style='oc',markersize=4); dfa.plot('TT','%SE8',style='oc',markersize=4)
    #obs moyennes
    tab0['HS'] = tab0.TT*tx; tab0.plot('TT','%',color='m')
    #tab20['HS'] = tab20.TT*tx; tab20.plot('HS','%',color='c')
    #sim
    #sim0.plot('HS','light0',color='b',linewidth=2)  
    #sim5.plot('TT','light5',color='r',linestyle='--',linewidth=2)     
    #sim20.plot('TT','light20',color='b',linewidth=2)
    #sim25.plot('TT','light25',color='b',linestyle='--',linewidth=2)
    
    return sim0
 
#-------------------------------------------------------------------------------------

# test parallelisation
'''
def func_star(a_b):
    return test_adel(*a_b)
def main():
    p = Pool()
    TT_stop=[0,600,400,200]
    #var = "Mercia"
    var = "Rht3"
    #var = "Tremie"
    p.map(func_star, itertools.izip(TT_stop, itertools.repeat(var)))
# if __name__=="__main__":
    # freeze_support()
    # main()
'''

