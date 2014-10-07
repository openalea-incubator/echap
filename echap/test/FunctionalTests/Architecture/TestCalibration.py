import pandas
import numpy
import os
import matplotlib.pyplot as plt
import math
import csv

from alinea.adel.WheatTillering import WheatTillering
import alinea.echap.architectural_data as archidb
import alinea.echap.architectural_reconstructions as rec

from alinea.echap.architectural_reconstructions import EchapReconstructions

from multiprocessing import Pool,freeze_support
import itertools

import openalea.plantgl.all as pgl
from alinea.adel.mtg_interpreter import *
from alinea.echap.weather_data import *
from alinea.caribu.caribu_star import diffuse_source, run_caribu
from alinea.adel.astk_interface import AdelWheat

plt.ion()

Reconst = EchapReconstructions()
HSconv = rec.HS_converter

def get_reconstruction(name='Mercia', **args):
    adel = Reconst.get_reconstruction(name, **args)
    return adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants

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
    
def test_axis_dynamics(name='Mercia', color='r'):

    dates = range(0,2500,100)
    adel = AdelWheat()
    res = AdelWheat.checkAxeDyn(adel, dates, density=1)

    conv = HSconv[name]
    res['HS'] = conv(res['TT'])
    
    res.plot('HS', 'nbaxes_present', style='--o'+color, label='nbaxes_present '+name)
    plt.ylim(ymin=0); plt.xlabel("HS")
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
#------------------------------------------------------------------------------------- 
# Intervalle de confiance
def mean(lst):
    return sum(lst) / float(len(lst))

def variance(lst):
    """
    Uses standard variance formula (sum of each (data point - mean) squared)
    all divided by number of data points
    """
    mu = mean(lst)
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
    
    return math.sqrt(v/n) * c

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

def compare_LAI(name='Mercia', n=30, aborting_tiller_reduction=1, **kwds): 

    if name is 'Mercia':
        color='r'
    elif name is 'Rht3':
        color='g'
    elif name is 'Tremie12':
        color='b'
    else:
        color = 'm'

    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, **kwds)
    conv = HSconv[name]
    sim = simLAI(adel, domain_area, convUnit, nplants)
    sim['HS'] = conv(sim.ThermalTime)
    sim.plot('HS','LAI_vert',color=color, label='LAI vert simule '+name)
    
    #sim['nbr_axe_tot'] = (( sim['Nbr.axe.tot.m2'] * sim['aire du plot']) / sim['Nbr.plant.perplot']) - 1
    #sim.plot('HS','nbr_axe_tot', style='--oy', label='Nbr axe tot '+name) 
    
    #donnees photo
    obs = archidb.PAI_photo_data()[name]
    obs['HS'] = conv(obs.TT_date)
        #mean
    obs_new = obs.groupby('HS').mean(); obs_new = obs_new.reset_index()
        #IC 
    IC = []
    lst_HS = obs['HS'].unique()
    grouped = obs.groupby('HS')
    for HS in lst_HS :
        obs_HS = grouped.get_group(HS)
        print 'HS = ', obs_HS
        IC.append(conf_int(obs_HS['PAI_vert_photo'], perc_conf=95))
    obs_new['IC'] = IC
        #plot mean with IC
    plt.errorbar(obs_new['HS'], obs_new['PAI_vert_photo'], yerr=obs_new['IC'], fmt='--o'+color, label = 'PAI vert photo mean+IC '+name)
    
    #donnees biomasse
    var = ['Tremie12', 'Tremie13']
    if name in var :
        obs = archidb.LAI_biomasse_data()[name]
        obs['HS'] = conv(obs.TT_date)
            #mean
        obs_new = obs.groupby('HS').mean(); obs_new = obs_new.reset_index()
                #IC 
        IC = []
        lst_HS = obs['HS'].unique()
        grouped = obs.groupby('HS')
        for HS in lst_HS :
            obs_HS = grouped.get_group(HS)
            IC.append(conf_int(obs_HS['LAI_vert_biomasse'], perc_conf=95))
        obs_new['IC'] = IC
            #plot mean with IC
        if name is 'Tremie12':
            color = 'c'
        else:
            color = '#F5A9F2'
        plt.errorbar(obs_new['HS'], obs_new['LAI_vert_biomasse'], yerr=obs_new['IC'], fmt='--o', color=color, label = 'LAI vert biomasse mean+IC '+name)
    
    plt.xlabel("HS")
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
    return sim
    
#------------------------------------------------------------------------------------- 
# Hauteur de couvert simule

def height(name='Mercia', n=30, aborting_tiller_reduction=1, **kwds):
    from alinea.astk.plantgl_utils import get_height
    
    if name is 'Mercia':
        color='r'
    elif name is 'Rht3':
        color='g'
    elif name is 'Tremie12':
        color='b'
    else:
        color='m'

    dd = range(500,2500,100)
    
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction=aborting_tiller_reduction, **kwds)
    sim = [adel.setup_canopy(age) for age in dd]
    max_h = []
    for g in sim:
        #scene = plot3d(g)
        #pgl.Viewer.display(scene)
        scene_geom = g.property('geometry')
        heights = get_height(scene_geom)
        max_height = max(map(numpy.min,heights.values()))
        max_h.append(max_height)

    h = pandas.DataFrame({'tt':dd, 'height':max_h})
    conv = HSconv[name]
    h['HS'] = conv(h['tt'])
    
    h.plot('HS', 'height', style='--o'+color, label = 'Height '+name)
    plt.xlabel("HS")
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
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
    gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = 4288, image_height = 2848, getImages=True, replicate=rep)

    return gc
    
def comp_TC(name='Mercia', n=30, zenith=0, aborting_tiller_reduction=1, **kwds): #zenith = 0 or 57
    conv = HSconv[name]

    if zenith==0:
        zen='0'; rep=1
    else:
        zen='57'; rep=2

    dd = range(400,2600,100)
    
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, **kwds)
    sim = [adel.setup_canopy(age) for age in dd]
    TC_sim = [draft_TC(g, adel, domain, zenith, rep) for g in sim]

    n=0; tc=[]; tc_sen=[]
    while n<len(TC_sim):
        tc.append(TC_sim[n]['green'])
        tc_sen.append(TC_sim[n]['senescent'])
        n = n + 1

    sim_green = pandas.DataFrame({'tt':dd, 'TCgreen':tc, 'TCsen':tc_sen})
    sim_green['HS'] = conv(sim_green.tt)

    # plot sim TC_tot, TC_green et senescence sur le meme graph
    sim_green['TCtot'] = sim_green['TCgreen'] + sim_green['TCsen']
    sim_green.plot('HS','TCtot',style='-ob', label='TC total sim '+name); sim_green.plot('HS','TCgreen',style='-og', label = 'TC green sim '+name); sim_green.plot('HS','TCsen',style='-oy', label = 'TC senescent sim '+name)
    
    #obs    
    obs = archidb.TC_data()[name+'_'+zen]
    obs['HS'] = conv(obs.TT)
        #mean
    obs_new = obs.groupby('HS').mean(); obs_new = obs_new.reset_index()
        #IC 
    IC = []
    lst_HS = obs['HS'].unique()
    grouped = obs.groupby('HS')
    for HS in lst_HS :
        obs_HS = grouped.get_group(HS)
        IC.append(conf_int(obs_HS['TC'], perc_conf=95))
    obs_new['IC'] = IC
        #plot mean with IC
    plt.errorbar(obs_new['HS'], obs_new['TC'], yerr=obs_new['IC'], fmt='or', label = 'TC mean+IC '+name)
    
    plt.xlabel("HS")
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
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
def graph_meteo(name='Mercia', n=30, aborting_tiller_reduction=1, **kwds):
    from alinea.astk.TimeControl import thermal_time # Attention la fonction marche seulement avec freq='H' !!! (pas en jour)
    conv = HSconv[name]
    
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
    if name is 'Tremie12': #Tremie 2011/2012
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
    
    dd = range(500,2500,1000) #dd = range(400,2600,300)
    
    # PARTIE SIM --------------------------------------------------------------------------------
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, **kwds)
    sim = [adel.setup_canopy(age) for age in dd]

    #list_level = [[0,'g'],[5,'y'],[20,'m'],[25,'r']]
    list_level = [[0,'g']]
    for level,c in list_level : 
        light_sim = 'light_sim_'+str(level); res_sim = 'sim_'+str(level)
        light_sim = [draft_light(g, adel, domain, z_level=level) for g in sim]
        res_sim = pandas.DataFrame({'TT':dd, 'light':light_sim})
        res_sim['HS'] = conv(res_sim.TT)
        #plot
        res_sim.plot('HS', 'light', style = '-'+c, linewidth=2, label = 'Light sim level = '+str(level)) #linewidth=2


    # PARTIE DONNEES ----------------------------------------------------------------------------
    #obs tous les points
    dfa = dfa.merge(bid); dfa['HS'] = conv(dfa.TT); dfa = dfa.sort(['HS'])
    if name is 'Mercia':
        for col in [1,2,3,4]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+g', label = 'Capteurs PAR sol') #markersize=4
        tab0['HS'] = conv(tab0.TT); tab0.plot('HS','%',color='g', label = 'Moyenne capteurs PAR sol')
        print tab0
    elif name is 'Rht3':
        for col in [5,6,7,8]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+g', label = 'Capteurs PAR sol') #markersize=4
        tab20['HS'] = conv(tab20.TT); tab20.plot('HS','%',color='g', label = 'Moyenne capteurs PAR sol')
    elif name is 'Tremie12':
        for col in [1,2,3,4]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+b',label=None) #markersize=4
        for col in [5,6,7,8]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+c',label=None) #markersize=4
        tab0['HS'] = conv(tab0.TT); tab0.plot('HS','%',color='b', label = 'Moyenne capteurs PAR sol')
        tab20['HS'] = conv(tab20.TT); tab20.plot('HS','%',color='c', label = 'Moyenne capteurs PAR 20cm')
    
    plt.xlabel("HS"); plt.ylabel("%")
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
    # BOITE DE PETRI ----------------------------------------------------------------------------
    if name is 'Tremie12':
        data_file1 = 'T1_20112012.csv'; data_file2 = 'T2_20112012.csv'
        header_row=['PETRI','Volume','Niveau','Bloc','ABSORBANCE','DILUTION','concentration(mg/l)','ConcentrationArrondie(mg/l)','quantiteRetenue(mg)','quantite(g/ha)','rapportPoucentage(sol/emis)']
        df1 = pandas.read_csv(data_file1, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal=',')
        df2 = pandas.read_csv(data_file2, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal=',')
        df1['TT'] = 1278; df2['TT'] = 1589
        print df1
        #df1['HS'] = conv(df1.TT); df2['HS'] = df2(res_sim.TT)
    
    
    #return sim