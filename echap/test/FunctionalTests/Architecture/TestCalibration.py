import pandas
import numpy
import os
import matplotlib.pyplot as plt
import math
import csv

from alinea.adel.WheatTillering import WheatTillering
import alinea.echap.architectural_data as archidb
import alinea.echap.interception_data as interceptdb
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
    
    return res
    
    #res.plot('HS', 'nbaxes_present', style='--o'+color, label='nbaxes_present '+name)
    #plt.ylim(ymin=0); plt.xlabel("HS")
    #plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
def multi_plot_test_axis_dynamics():
    fig, axes = plt.subplots(nrows=2, ncols=2)
    ax0, ax1, ax2, ax3 = axes.flat
    
    varieties = [['Mercia',ax0],['Rht3',ax1],['Tremie12',ax2],['Tremie13',ax3]]
    for name,ax in varieties :
        res = test_axis_dynamics(name=name)
        ax.plot(res['HS'], res['nbaxes_present'], '-ob', label = 'nbaxes_present '+name)
        ax.set_xlim([0, 18]); ax.set_ylim([0, 10]) 
        ax.legend(numpoints=1, bbox_to_anchor=(1.1,1.1), prop={'size': 9})
        
    fig.suptitle("Axis dynamics")
        
    
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
     
    dd = range(500,700,200)
    
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
# Image plante
def silhouette(name='Mercia', n=1, aborting_tiller_reduction=1, **kwds):   
    if name is 'Mercia':
        HS = [10.9, 15.73]
    elif name is 'Rht3':
        HS = [10.57, 15.4]
    elif name is 'Tremie12':
        HS = [7.42, 10.52, 13.06, 17.85]
    else :
        HS = [8.36, 9.1, 9.46]
    # conversion HS en TT
    conv = HSconv[name]
    TT = conv.TT(HS)

    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, **kwds)
    
    for age in TT :
        g = adel.setup_canopy(age)
        HS = conv(age)
        #scene = plot3d(g)
        gms = adel.get_axis(g, 'plant1')
        scene = adel.plot(gms)
        fname = 'plot_'+name+'_'+str(HS)+'.png'
        pgl.Viewer.display(scene)
        #pgl.Viewer.camera.lookAt((0,150,0),(0,0,0))
        pgl.Viewer.frameGL.saveImage(fname)
        
def plot_scan_obs_sim_surface(name='Tremie12', n=30): # Pour Tremie12 et Tremie13
    #fichier obs
    df_obs = archidb.treatment_scan(name)
    
    #simulation
    #if name is 'Mercia':
    #    HS = [9.74, 12.8]
    #elif name is 'Rht3':
    #    HS = [9.15, 12.48]
    if name is 'Tremie12':
        HS = [7.55, 10.15, 10.98, 12.63]
    else:
        HS = [8.36, 8.7, 9.7, 11.04]
    #conversion HS en TT
    conv = HSconv[name]
    dd = conv.TT(HS)  
    #sim pour chacun des TT
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n)    
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    out = out[out['axe_id']=='MS']; out = out[out['Slv']>0]
    def _fun(sub):
        sub['nmax'] = sub['numphy'].max()
        return sub
    grouped = out.groupby(['TT','plant'], as_index=False)
    df = grouped.apply(_fun)
    df['ntop_cur'] = df['nmax'] - df['numphy'] + 1
    lst_TT = df['TT'].unique()
    date = []; moy_numphy = []; ntop = []; Slv = []; Slvgreen = []
    for TT in lst_TT:
        dt = df[df['TT']==TT]
        data = dt.groupby(['ntop_cur'], as_index=False).mean()   
        nbr = data['ntop_cur'].count(); cpt = 0
        hs = conv(TT)
        while cpt<nbr:
            date.append(hs)
            moy_numphy.append(data['numphy'][cpt])
            ntop.append(data['ntop_cur'][cpt])
            Slv.append(data['Slv'][cpt])
            Slvgreen.append(data['Slvgreen'][cpt])
            cpt += 1
    surfSet = zip(date, moy_numphy, ntop, Slv, Slvgreen)
    df_fin = pandas.DataFrame(data = surfSet, columns=['HS', 'moyenne numphy', 'ntop_cur', 'Slv', 'Slvgreen'])

    #plot obs/sim sur le meme diagramme
    df_fin.HS=map(str,df_fin.HS)
    df_obs.HS=map(str,df_obs.HS)
    df_all = df_fin.merge(df_obs.ix[:,['HS','ntop_cur','Area A_bl']], how='outer')
    df_all = df_all[df_all['ntop_cur']<=5]

    bar_width = 0.4; opacity = 0.4
    fig, axes = plt.subplots(nrows=1, ncols=len(df_all['HS'].unique()))
    val = df_all['HS'].unique()
        
    for x, HS in enumerate(val):
        print HS
        df_fin = df_all[df_all['HS']==HS]
        n_groups = len(df_fin['ntop_cur'].unique())
        index = numpy.arange(n_groups)      
        rects1 = axes[x].bar(index, df_fin['Area A_bl'], bar_width,
                     alpha=opacity,
                     color='b')
        rects2 = axes[x].bar(index + bar_width, df_fin['Slv'], bar_width,
                     alpha=opacity,
                     color='r')   
        # Si on veut ajouter Slvgreen au graph, attention de modifier bar_width par 0.33
        #rects3 = axes[x].bar(index + bar_width*2, df_fin['Slvgreen'], bar_width,alpha=opacity,color='y')
        
        # Mise en forme
        axes[x].set_ylim(0, 35); axes[x].set_xlim(0, df_all['ntop_cur'].max())          
        axes[x].set_xticks(index+bar_width)
        axes[x].set_xticklabels( df_fin['ntop_cur'].tolist() )
        if x == 0:
            axes[x].set_xlabel('ntop_cur'); axes[x].set_ylabel('area (cm2)')
        axes[x].text(0.4, 33, 'HS : '+HS, bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10})
        fig.suptitle('Surface obs/sim pour '+name, fontsize=10)
   
    #axes[x].legend( (rects1[0], rects2[0], rects3[0]), ('Obs', 'Sim Slv', 'Sim Slvgreen') )
    axes[x].legend((rects1[0], rects2[0]), ('Obs', 'Sim') )
    
    
    return df_all


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
    gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = 2144, image_height = 1424, getImages=True, replicate=rep)

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
    sim_green.plot('HS','TCtot',style='-ob', label='TC total sim '+name)
    sim_green.plot('HS','TCgreen',style='-og', label = 'TC green sim '+name)
    sim_green.plot('HS','TCsen',style='-oy', label = 'TC senescent sim '+name)
    
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
    
    # PARTIE OBS ------------------------------------------
    if name is 'Mercia' or name is 'Rht3': 
        data_file = 'METEO_stationINRA_20102011.csv'
        # correspondance entre date et TT
        met = Boigneville_2010_2011()
        seq = pandas.date_range(start="2010-10-15", end="2011-07-18", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2010-10-15", end="2011-07-18")
        bid = bid[seq]
    if name is 'Tremie12': 
        data_file = 'METEO_stationINRA_20112012.csv'
        # correspondance entre date et TT
        met = Boigneville_2011_2012()
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18")
        bid = bid[seq]
    if name is 'Tremie13': 
        data_file = 'METEO_stationINRA_20122013.csv'
        # correspondance entre date et TT
        met = Boigneville_2012_2013()
        seq = pandas.date_range(start="2012-10-29", end="2013-08-01", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2012-10-29", end="2013-08-01")
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
    
    dd = range(500,2500,100) #dd = range(400,2600,300)
    
    # PARTIE DONNEES ----------------------------------------------------------------------------
    #obs tous les points
    dfa = dfa.merge(bid); dfa['HS'] = conv(dfa.TT); dfa = dfa.sort(['HS'])
    if name is 'Mercia':
        for col in [1,2,3,4]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+g', label = '_nolegend_') #markersize=4
        tab0['HS'] = conv(tab0.TT); tab0.plot('HS','%',color='g', label = 'Moyenne capteurs PAR sol')
        list_level = [[0,'g']]
    elif name is 'Rht3':
        for col in [5,6,7,8]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+g', label = '_nolegend_') #markersize=4
        tab20['HS'] = conv(tab20.TT); tab20.plot('HS','%',color='g', label = 'Moyenne capteurs PAR sol')
        list_level = [[0,'g']]
    elif name is 'Tremie12' or name is 'Tremie13':
        for col in [1,2,3,4]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+g', label = '_nolegend_') #markersize=4
        for col in [5,6,7,8]:
            colf = '%SE'+str(col)
            dfa.plot('HS', colf, style='+c', label = '_nolegend_') #markersize=4
        tab0['HS'] = conv(tab0.TT); tab0.plot('HS','%',color='g', label = 'Moyenne capteurs PAR sol')
        tab20['HS'] = conv(tab20.TT); tab20.plot('HS','%',color='c', label = 'Moyenne capteurs PAR 20cm')
        list_level = [[0,'g'],[20,'c']]
    
    # BOITE DE PETRI ----------------------------------------------------------------------------
    if name is 'Tremie12' or name is 'Tremie13':    
        if name is 'Tremie12':
            data_file_T1 = 'T1_20112012.csv'; data_file_T2 = 'T2_20112012.csv'
            date_T1 = 1277.8; date_T2 = 1589.4
        if name is 'Tremie13':
            data_file_T1 = 'T1_20122013.csv'; data_file_T2 = 'T2_20122013.csv'
            date_T1 = 1041.8; date_T2 = 1309.4
        header_row=['PETRI','Volume','Niveau','Bloc','ABSORBANCE','DILUTION','concentration(mg/l)','ConcentrationArrondie(mg/l)','quantiteRetenue(mg)','quantite(g/ha)','rapportPoucentage(sol/emis)']
        df1 = pandas.read_csv(data_file_T1, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=0, decimal=',')
        df2 = pandas.read_csv(data_file_T2, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=0, decimal=',')
        # bon format pourcentage
        df1['rapportPoucentage(sol/emis)'] = df1['rapportPoucentage(sol/emis)'] / 100.
        df2['rapportPoucentage(sol/emis)'] = df2['rapportPoucentage(sol/emis)'] / 100.
        # ajout TT et conversion en HS
        df1['TT'] = date_T1; df2['TT'] = date_T2
        df1['HS'] = conv(df1.TT); df2['HS'] = conv(df2.TT)
        # on filtre sur niveau = sol
        grouped1 = df1.groupby('Niveau'); grouped2 = df2.groupby('Niveau')
        petri1 = grouped1.get_group('sol'); petri2 = grouped2.get_group('sol')
        # calcul mean et IC
        lst1 = petri1['rapportPoucentage(sol/emis)']; lst2 = petri2['rapportPoucentage(sol/emis)']
        mean1 = mean(lst1); mean2 = mean(lst2)
        IC1 = conf_int(lst1, perc_conf=95); IC2 =conf_int(lst2, perc_conf=95)
        plt.errorbar(df1['HS'][1], mean1, yerr=IC1, fmt='or', label = 'Petri T1, mean+IC, '+name)
        plt.errorbar(df2['HS'][1], mean2, yerr=IC2, fmt='^r', label = 'Petri T2, mean+IC, '+name)
    else :
        print '--Pas de donnees de boites de Petri pour '+name+'--'
        
    # PARTIE SIM --------------------------------------------------------------------------------
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, **kwds)
    sim = [adel.setup_canopy(age) for age in dd]

    for level,c in list_level : 
        light_sim = 'light_sim_'+str(level); res_sim = 'sim_'+str(level)
        light_sim = [draft_light(g, adel, domain, z_level=level) for g in sim]
        res_sim = pandas.DataFrame({'TT':dd, 'light':light_sim})
        res_sim['HS'] = conv(res_sim.TT)
        #plot
        res_sim.plot('HS', 'light', style = '-'+c, linewidth=2, label = 'Light sim level = '+str(level))
 
    plt.ylim(ymax=1.); plt.xlabel("HS"); plt.ylabel("%")
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})
