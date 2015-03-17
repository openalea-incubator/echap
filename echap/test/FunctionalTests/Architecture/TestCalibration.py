# -*- coding: latin1 -*- 
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
import alinea.echap.interception_data as idata

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
     
    dd = range(500,2500,100)
    
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res

def compare_LAI(name='Mercia', n_sim=1, n_plt=30, aborting_tiller_reduction=1, seed=1, density=1, **kwds): 

    if name is 'Mercia':
        color='r'; T1=9.74; T2=12.8
    elif name is 'Rht3':
        color='g'; T1=9.15; T2=12.48
    elif name is 'Tremie12':
        color='b'; T1=10.98; T2=12.63
    else:
        color='m'; T1=8.7; T2=11.04
        
    if density==0.1:
        picto = '-+'
    elif density==1:
        picto = '-o'
    elif density==5:
        picto = '-:'
    elif density==8:
        picto = '-^'
    elif density==10:
        picto='-<'
       
    conv = HSconv[name]
    
    x=1; sim_all = pandas.DataFrame()
    while x <= n_sim:
        adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n_plt, aborting_tiller_reduction = aborting_tiller_reduction, seed=seed, stand_density_factor = {name:density}, **kwds)   
        df_sim = simLAI(adel, domain_area, convUnit, nplants)
        df_sim['numero_sim']=str(x)
        sim_all = sim_all.append(df_sim)
        x += 1
      
    sim_std = sim_all.groupby(['ThermalTime']).std()
    sim_std = sim_std.reset_index()
    sim = sim_all.groupby(['ThermalTime']).mean()
    sim = sim.reset_index()
    sim['std_LAI_vert'] = sim_std['LAI_vert']
      
    sim['HS'] = conv(sim.ThermalTime)
    
    if n_sim>1:
        plt.errorbar(sim['HS'], sim['LAI_vert'], yerr=sim['std_LAI_vert'], fmt='--o'+color, label = '5x30plts + IC '+name)
    else :
        sim.plot('HS','LAI_vert',style='-'+color, label='LAI vert simule '+name+', density = '+str(density))
    
    #sim.plot('HS','LAI_tot',color=color, label='LAI tot simule '+name+', density = '+str(density))
    
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
        #print 'HS = ', obs_HS
        IC.append(conf_int(obs_HS['PAI_vert_photo'], perc_conf=95))
    obs_new['IC'] = IC
        #plot mean with IC
    #plt.errorbar(obs_new['HS'], obs_new['PAI_vert_photo'], yerr=obs_new['IC'], fmt='--o'+color, label = 'PAI vert photo mean+IC '+name)
    
    #traitement T1 et T2
    plt.axvline(x=T1, color=color, alpha=0.4)
    plt.axvline(x=T2, color=color, alpha=0.4)
    
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
    
    plt.xlabel("haun stage")
    plt.ylabel("LAI")
    plt.grid(False) # enlever la grille
    #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
    return sim
    
#-----------------------------------------------------------------------------
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
    
    h.plot('HS', 'height', style='-'+color, label = name)
    plt.xlabel("haun stage"); plt.ylabel("hauteur")
    plt.grid(False) # enlever la grille
    #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
    return h
    
#-----------------------------------------------------------------------------
# Image plante
def silhouette(name='Mercia', n=1, aborting_tiller_reduction=1, seed=1, **kwds):   
    # visualiser une seule plante aux differents HS correspondants à ceux de la fonction plot_HS dans dye_interception.py et enregistrer image
    # !!! pas de positionnement de la camera : il faut regler la vue sur plantGL et relancer la fonction pour avoir de belles images et pouvoir notamment monter le film de la croissance de la plante
    if name is 'Mercia':
        HS = [9.34,9.54,9.74,9.94,10.14,12.4,12.6,12.8,13,13.2,13.3,13.8,14.3,14.8,15.3]
    elif name is 'Rht3':
        HS = [8.75,8.95,9.15,9.35,9.55,12.08,12.28,12.48,12.68,12.88,12.98,13.48,13.98,14.48,14.98]
    elif name is 'Tremie12':
        HS = [10.58,10.78,10.98,11.18,11.38,12.23,12.43,12.63,12.83,13.03,13.13,13.63,14.13,14.63,15.13]
    else :
        HS = [8.3,8.5,8.7,8.9,9.1,10.64,10.84,11.04,11.24,11.44,11.54,12.04,12.54,13.04,13.54]
    
    # conversion HS en TT
    conv = HSconv[name]
    TT = conv.TT(HS)

    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, seed=seed, **kwds)
    
    for age in TT :
        g = adel.setup_canopy(age)
        HS = conv(age)
        gms = adel.get_axis(g, 'plant1')
        scene = adel.plot(gms)
        fname = 'plot_'+name+'_'+str(HS)+'_seed'+str(seed)+'.png'
        pgl.Viewer.display(scene)
        pgl.Viewer.frameGL.saveImage(fname)
        
def silhouette_dyn(name='Mercia', n=30, HS=7, aborting_tiller_reduction=1, seed=1, density=1, **kwds):   
    '''visualiser les 30 plantes, verifier le changement de density dans la reconstruction et enregistrer image
    '''
    #HS = numpy.arange(7,8,1)
    
    # conversion HS en TT
    conv = HSconv[name]
    TT = conv.TT(HS)

    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, seed=seed, stand_density_factor = {name:density}, **kwds)
    
    #for age in TT :
    g = adel.setup_canopy(TT)
    #HS = conv(TT)
    scene = plot3d(g)
    fname = 'plot_'+name+'_'+str(HS)+'_seed'+str(seed)+'.png'
    pgl.Viewer.display(scene)
    pgl.Viewer.frameGL.saveImage(fname)
        
    return adel, g
        
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

def plot_scan_obs_dimfactor(name='Tremie12', n=30): # Pour Tremie12 et Tremie13
    #fichier obs
    df_obs_all = pandas.DataFrame()
    for name in ['Tremie12','Tremie13']:
        df_obs = archidb.treatment_scan(name)
        df_obs['var'] = name
        df_obs_all = df_obs_all.append(df_obs)
    
    #simulation
    '''if name is 'Tremie12':
        HS = [7.55, 10.15, 10.98, 12.63]
    else:
        HS = [8.36, 8.7, 9.7, 11.04]'''
    #conversion HS en TT
    '''conv = HSconv[name]
    dd = conv.TT(HS)  '''

    df_dim_all = pandas.DataFrame()
    # dim factor
    d = rec.leaf_fits(); surf_dim = pandas.DataFrame()
    #if name is 'Tremie12':
    for name in ['Tremie12','Tremie13']:   
        dim_fits = archidb.dimension_fits()
        if name=='Tremie12':
            df_dim = dim_fits[name]
        else:
            df_dim = dim_fits[name]
        data = d[name]
        fact_haut = data.form_factor()[2]
        fact_bas = data.form_factor()[1]
        if name=='Tremie12':
            nff_lst=[12,13]
        else:
            nff_lst=[11,12]
        for nff in nff_lst:
            df_dim[nff] = df_dim[nff][df_dim[nff]['index_phytomer']>=8]
            df_dim[nff]['var'] = name; df_dim[nff]['nff'] = nff
            df_dim[nff]['surface_bas'] = df_dim[nff]['L_blade']*df_dim[nff]['W_blade']*fact_bas
            df_dim[nff]['surface_haut'] = df_dim[nff]['L_blade']*df_dim[nff]['W_blade']*fact_haut
            df_dim_all = df_dim_all.append(df_dim[nff])
            
    bar_width = 0.4; opacity = 0.5
    fig, axes = plt.subplots(nrows=1, ncols=5)
    val = [[0,'Tremie12',None],[1,'Tremie12',12],[2,'Tremie12',13],[3,'Tremie13',None],[4,'Tremie13',12]]
    for x, var, nff in val :
        if x==1 or x==2 or x==4 :
            df_dim_val = df_dim_all[df_dim_all['var']==var]
            df_dim_val = df_dim_val[df_dim_val['nff']==nff]
            df_dim_val = df_dim_val.sort(['index_phytomer'])
            n_groups = len(df_dim_val)
            index = numpy.arange(n_groups)
            if nff==13:
                rects1 = axes[x].bar(index, df_dim_val['surface_bas'], bar_width, alpha=opacity, color=['c','y','m','g','r','b'])
                rects2 = axes[x].bar(index + bar_width, df_dim_val['surface_haut'], bar_width, color=['c','y','m','g','r','b'])
            else:
                rects1 = axes[x].bar(index, df_dim_val['surface_bas'], bar_width, alpha=opacity, color=['c','y','m','g','r'])
                rects2 = axes[x].bar(index + bar_width, df_dim_val['surface_haut'], bar_width, color=['c','y','m','g','r'])
            # Mise en forme
            axes[x].set_ylim(0, 40)
            axes[x].set_xlim(0, len(df_dim_val))                 
            axes[x].set_xticks(index+bar_width)
            df_dim_val['label'] = df_dim_val['index_phytomer'].astype(str)
            axes[x].set_xticklabels( df_dim_val['label'].tolist(), rotation=90, fontsize='small' )
            if x == 1 or x==2 or x==3 :
                axes[x].set_xlabel('index_phytomer')
            if x == 1:
                axes[x].text(0.4, 38, 'TREMIE 12 - NFF = '+str(nff), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
            elif x == 2:
                axes[x].text(0.4, 38, 'TREMIE 12 - NFF = '+str(nff), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
            elif x == 4:
                axes[x].text(0.4, 38, 'TREMIE 13 - NFF = '+str(nff), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
        if x==0 or x==3:
            df_obs_val = df_obs_all[df_obs_all['var']==var]
            df_obs_val = df_obs_val.sort(['ntop_cur', 'HS'], ascending=False)
            n_groups = len(df_obs_val)
            index = numpy.arange(n_groups)
            if x == 0 :
                rects1 = axes[x].bar(index, df_obs_val['Area A_bl'], bar_width, alpha=opacity, color=['c','y','y','m','m','m','m','g','g','g','g','r','r','r'])
            if x == 3 :
                rects1 = axes[x].bar(index, df_obs_val['Area A_bl'], bar_width, alpha=opacity, color=['c','c','y','y','m','m','g','g','r','r'])
            # Mise en forme
            axes[x].set_ylim(0, 40)
            axes[x].set_xlim(0, len(df_obs_val))                 
            axes[x].set_xticks(index+bar_width)
            df_obs_val['label'] = df_obs_val['HS'].astype(str) + ' / '+ df_obs_val['ntop_cur'].astype(str)
            axes[x].set_xticklabels( df_obs_val['label'].tolist(), rotation=90, fontsize='small' )
            axes[x].set_xlabel('HS / ntop_cur')
            if x == 0:
                axes[x].text(0.4, 38, 'TREMIE 12 - SCAN', bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
            if x == 3:
                axes[x].text(0.4, 38, 'TREMIE 13 - SCAN', bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            
#------------------------------------------------------------------------------------- 
# Taux de couverture

def draft_TC(g, adel, domain, zenith, rep, scale = 1):
    from alinea.adel.postprocessing import ground_cover
    
    #modelisation afin de voir si erreur
    #scene=plot3d(g)
    #pgl.Viewer.display(scene)
    
    echap_top_camera =  {'type':'perspective', 'distance':200., 'fov':50., 'azimuth':0, 'zenith':zenith}
    #high resolution
    #gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = 4288, image_height = 2848, getImages=True, replicate=rep)
    gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = int(2144 * scale), image_height = int(1424*scale), getImages=True, replicate=rep)

    return gc
    
def comp_TC(name='Mercia', n=30, zenith=0, dd = range(400,2600,100), scale = 1, aborting_tiller_reduction=1, seed=1, density=1, **kwds): #zenith = 0 or 57
    conv = HSconv[name]

    if zenith==0:
        zen='0'; rep=1
    else:
        zen='57'; rep=2
        
    if name is 'Mercia':
        color='r'; T1=9.74; T2=12.8
    elif name is 'Rht3':
        color='g'; T1=9.15; T2=12.48
    elif name is 'Tremie12':
        color='b'; T1=10.98; T2=12.63
    else:
        color = 'm'; T1=8.7; T2=11.04
        
    if density==0.5:
        picto = '+'
    elif density==0.625:
        picto = 'o'
    elif density==0.75:
        picto = ':'
    elif density==0.875:
        picto = '^'
    elif density==1:
        picto = '-'

    
    
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, seed=seed, stand_density_factor = {name:density}, **kwds)
    sim = (adel.setup_canopy(age) for age in dd)
    TC_sim = [draft_TC(g, adel, domain, zenith, rep, scale=scale) for g in sim]

    n=0; tc=[]; tc_sen=[]
    while n<len(TC_sim):
        tc.append(TC_sim[n]['green'])
        tc_sen.append(TC_sim[n]['senescent'])
        n = n + 1

    sim_green = pandas.DataFrame({'tt':dd, 'TCgreen':tc, 'TCsen':tc_sen})
    sim_green['HS'] = conv(sim_green.tt)

    # plot sim TC_tot, TC_green et senescence sur le meme graph
    #sim_green['TCtot'] = sim_green['TCgreen'] + sim_green['TCsen']
    #sim_green.plot('HS','TCtot',style='--'+color, label='TC total sim '+name+', density = '+str(density))
    sim_green.plot('HS','TCgreen',style=picto+color, label = 'TC green sim '+name+', density = '+str(density))
    #sim_green.plot('HS','TCsen',style='-y', label = 'TC senescent sim '+name)
    
    #traitement T1 et T2 + ajout trait horizontal y=0.95
    plt.axvline(x=T1, color=color, alpha=0.4)
    plt.axvline(x=T2, color=color, alpha=0.4)
    #if zenith==57:
    #    plt.axhline(y=0.95, color='k')
    
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
    plt.errorbar(obs_new['HS'], obs_new['TC'], yerr=obs_new['IC'], fmt='--o'+color, label = 'TC mean+IC '+name)
    
    plt.xlabel("haun stage")
    plt.ylabel("TC")
    plt.grid(False) # enlever la grille
    #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
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
def graph_meteo(name='Mercia', level=5, n=30, aborting_tiller_reduction=1, seed=1, **kwds):
    # level prend en arg 0 5 20 25 cm
    from alinea.astk.TimeControl import thermal_time # Attention la fonction marche seulement avec freq='H' !!! (pas en jour)
    conv = HSconv[name]
    
    # PARTIE OBS ------------------------------------------
    if name is 'Mercia' or name is 'Rht3': 
        # correspondance entre date et TT
        met = Boigneville_2010_2011()
        seq = pandas.date_range(start="2010-10-15", end="2011-07-18", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2010-10-15", end="2011-07-18")
        bid = bid[seq]
    if name is 'Tremie12': 
        # correspondance entre date et TT
        met = Boigneville_2011_2012()
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18")
        bid = bid[seq]
    if name is 'Tremie13': 
        # correspondance entre date et TT
        met = Boigneville_2012_2013()
        seq = pandas.date_range(start="2012-10-29", end="2013-08-01", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2012-10-29", end="2013-08-01")
        bid = bid[seq]
        
    dfa, tab0, tab20 = mat_ray_obs(name)
        
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
    
    # PARTIE DONNEES ----------------------------------
    #obs tous les points
    dfa = dfa.merge(bid); dfa['HS'] = conv(dfa.TT); dfa = dfa.sort(['HS'])
    if name is 'Mercia':
        for col in [1,2,3,4]:
            colf = '%SE'+str(col)
            #dfa.plot('HS', colf, style='+g', label = '_nolegend_') #markersize=4
        tab0['HS'] = conv(tab0.TT)
        tab0['1-%'] = 1 - tab0['%']
        if level==0 or level==5:
            tab0.plot('HS','1-%',color='r', label = 'Moyenne capteurs PAR sol')
        list_level = [[level,'r']]
    elif name is 'Rht3':
        for col in [5,6,7,8]:
            colf = '%SE'+str(col)
            #dfa.plot('HS', colf, style='+g', label = '_nolegend_') #markersize=4
        tab20['HS'] = conv(tab20.TT)
        tab20['1-%'] = 1 - tab20['%']
        if level==0 or level==5:
            tab20.plot('HS','1-%',color='g', label = 'Moyenne capteurs PAR sol')
        list_level = [[level,'g']]
    elif name is 'Tremie12' or name is 'Tremie13':
        if name=='Tremie12':
            color='b'
        else:
            color='m'
        for col in [1,2,3,4]:
            colf = '%SE'+str(col)
            #dfa.plot('HS', colf, style='+g', label = '_nolegend_') #markersize=4
        for col in [5,6,7,8]:
            colf = '%SE'+str(col)
            #dfa.plot('HS', colf, style='+c', label = '_nolegend_') #markersize=4
        tab0['HS'] = conv(tab0.TT)
        tab0['1-%'] = 1 - tab0['%']
        if level==0 or level==5:
            tab0.plot('HS','1-%',color=color, label = 'Moyenne capteurs PAR sol')
        tab20['HS'] = conv(tab20.TT)
        tab20['1-%'] = 1 - tab20['%']
        if level==20 or level==25:
            tab20.plot('HS','1-%',color=color, label = 'Moyenne capteurs PAR 20cm')
        list_level = [[level,color]]
        
    if name is 'Mercia':
        color='r'; T1=9.74; T2=12.8
    elif name is 'Rht3':
        color='g'; T1=9.15; T2=12.48
    elif name is 'Tremie12':
        color='b'; T1=10.98; T2=12.63
    else:
        color = 'm'; T1=8.7; T2=11.04
        
    #traitement T1 et T2
    plt.axvline(x=T1, color=color, alpha=0.4)
    plt.axvline(x=T2, color=color, alpha=0.4)
    
    # BOITE DE PETRI ----------------------------------------------------------------------------
    if name is 'Tremie12' or name is 'Tremie13':    
        petri_T1, petri_T2 = idata.Petri_data(name)
        # calcul mean et IC
        lst1 = petri_T1['rapportPoucentage(sol/emis)']
        lst2 = petri_T2['rapportPoucentage(sol/emis)']
        mean1 = mean(lst1); mean2 = mean(lst2)
        IC1 = conf_int(lst1, perc_conf=95)
        IC2 =conf_int(lst2, perc_conf=95)
        #plt.errorbar(petri_T1['HS'][1], mean1, yerr=IC1, fmt='or', label = 'Petri T1, mean+IC, '+name)
        #plt.errorbar(petri_T2['HS'][1], mean2, yerr=IC2, fmt='^r', label = 'Petri T2, mean+IC, '+name)
    else :
        print '--Pas de donnees de boites de Petri pour '+name+'--'
        
    # PARTIE SIM --------------------------------------------------------------------------------
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, aborting_tiller_reduction = aborting_tiller_reduction, seed=seed, **kwds)
    sim = (adel.setup_canopy(age) for age in dd)

    for level,c in list_level : 
        light_sim = 'light_sim_'+str(level); res_sim = 'sim_'+str(level)
        light_sim = [draft_light(g, adel, domain, z_level=level) for g in sim]
        res_sim = pandas.DataFrame({'TT':dd, 'light':light_sim})
        res_sim['HS'] = conv(res_sim.TT)
        res_sim['1-light'] = 1 - res_sim['light']
        #plot
        res_sim.plot('HS', '1-light', style = '-'+c, linewidth=2, label = '1 - light sim level = '+str(level))
 
    plt.xlim(xmin=0); plt.ylim(ymin=0)
    plt.ylim(ymax=1.); plt.xlabel("haun stage"); plt.ylabel("1-%")
    plt.grid(False) # enlever la grille
    #plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})

    
def plot_sup(graph1=False, graph2=False, graph3=False, graph4=False, graph5=True, stade='inf'): #stade = inf ou sup au max de LAI pour les graphs 1 et 2 seulement
    ''' 
    graph 1 = TCvert obs / LAIvert obs + TC0 vert sim/ LAI vert sim
    graph 2 = TCtot sim / LAItot sim
    graph 3 = TCvert obs / 1_rayonnement data 
    graph 4 = LAIvert obs / 1_rayonnement data
    graph 5 = TC 57 vert / LAI vert
    '''
    var_lst=['Mercia', 'Rht3', 'Tremie12', 'Tremie13']

    if graph1==True:
        for var in var_lst:
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            
            #LAI sim
            sim_lai = 'files/sim_lai_'+var+'.csv'
            sim_lai = pandas.read_csv(sim_lai, decimal='.', sep=',')
            #TC sim
            sim_tc = 'files/sim_tc_'+var+'.csv'
            sim_tc = pandas.read_csv(sim_tc, decimal='.', sep=',')
            #merge LAI sim + TC sim
            sim = sim_tc.merge(sim_lai)
            #decoupage selon max value du LAI
            max_lai = sim['LAI_vert'].max()
            max_tt = int(sim.TT[sim['LAI_vert']==max_lai])
           
            #LAI vert data
            obs_lai = archidb.PAI_photo_data()[var]    
            obs_lai = obs_lai.rename(columns={'TT_date': 'TT'}) #rename column for merge
            obs_new_lai = obs_lai.groupby('TT').mean() #mean by TT
            obs_new_lai = obs_new_lai.reset_index()
            #TC vert data
            obs_tc = archidb.TC_data()[var+'_0']
            obs_new_tc = obs_tc.groupby('TT').mean() #mean by TT
            obs_new_tc = obs_new_tc.reset_index()
            # merge LAI vert obs + TC vert obs
            obs = obs_new_tc.merge(obs_new_lai)
            #LAI vert biomasse
            if var=='Tremie12':
                obs_bio = archidb.LAI_biomasse_data()[var]
                obs_bio_lai = obs_bio.rename(columns={'TT_date': 'TT'})
                obs_bio_lai = obs_bio_lai.groupby('TT').mean() #mean by TT
                obs_bio_lai = obs_bio_lai.reset_index()
                obs = obs.merge(obs_bio_lai)
            if var=='Tremie13':
                obs_bio = archidb.LAI_biomasse_data()[var]
                obs_bio_lai = obs_bio.rename(columns={'TT_date': 'TT'})
                obs_bio_lai = obs_bio_lai.groupby('TT').mean() #mean by TT
                obs_bio_lai = obs_bio_lai.reset_index()
                add = float((994 * obs_bio_lai.LAI_vert_biomasse[obs_bio_lai['TT']==941]) / 941)
                obs_bio_lai.ix[0,0]=994 # grosse bidouille
                obs_bio_lai.ix[0,1]=add
                obs = obs.merge(obs_bio_lai, how='outer')
                        
            #plot <= max_lai
            if stade=='inf':
                sim_min = sim[sim['TT']<=max_tt]
                sim_min.plot('LAI_vert', 'TCgreen', style='-^'+color, label='Sim '+var, alpha=0.4)
                obs_min = obs[obs['TT']<=max_tt]
                obs_min.plot('PAI_vert_photo', 'TC', style='--o'+color, label='Obs '+var) 
                if var=='Tremie12':
                    obs_min.plot('LAI_vert_biomasse', 'TC', style='v'+color, label='Biomasse '+var) 
                elif var=='Tremie13':
                    obs_min.plot('LAI_vert_biomasse', 'TC', style='vy', label='Biomasse '+var)
            #plot >= max_lai
            else:
                sim_max = sim[sim['TT']>=max_tt]
                sim_max.plot('LAI_vert', 'TCgreen', style='-^'+color, label='Sim '+var, alpha=0.4)
                obs_max = obs[obs['TT']>=max_tt]
                obs_max.plot('PAI_vert_photo', 'TC', style='--o'+color, label='Obs '+var)
                if var=='Tremie12':
                    obs_max.plot('LAI_vert_biomasse', 'TC', style='v'+color, label='Biomasse '+var) 
            
            #mise en forme
            #plt.ylim(ymin=0); plt.xlim(xmin=0); plt.ylim(ymax=1); plt.xlim(xmax=5)
            plt.xlabel("LAI vert"); plt.ylabel("TC0 vert")
            plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
           
    if graph2==True:
        for var in var_lst:
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            
            #LAI sim            
            adel, domain, domain_area, convUnit, nplants = get_reconstruction(name=var, nplants = 30, aborting_tiller_reduction = 1, seed=1)   
            sim_lai = simLAI(adel, domain_area, convUnit, nplants)
            sim_lai = sim_lai.rename(columns={'ThermalTime' : 'TT'})            
            #TC sim, ajout colonne TC total    
            dd = range(400,2600,100)    
            adel, domain, domain_area, convUnit, nplants = get_reconstruction(name=var, nplants = 30, aborting_tiller_reduction = 1, seed=1)
            sim = (adel.setup_canopy(age) for age in dd)
            TC_sim = [draft_TC(g, adel, domain, zenith=0, rep=1) for g in sim]
            n=0; tc=[]; tc_sen=[]
            while n<len(TC_sim):
                tc.append(TC_sim[n]['green'])
                tc_sen.append(TC_sim[n]['senescent'])
                n = n + 1
            sim_tc = pandas.DataFrame({'TT':dd, 'TCgreen':tc, 'TCsen':tc_sen})
            sim_tc['TCtot'] = sim_tc['TCgreen'] + sim_tc['TCsen']
            
            #merge LAI sim + TC sim
            sim = sim_tc.merge(sim_lai)
            #decoupage selon max value du LAI
            max_lai = sim['LAI_tot'].max()
            max_tt = int(sim.TT[sim['LAI_tot']==max_lai])
            
            #plot <= max_lai
            if stade=='inf':
                sim_min = sim[sim['TT']<=max_tt]
                sim_min.plot('LAI_tot', 'TCtot', style='-^'+color, label='Sim '+var, alpha=0.4) 
            #plot >= max_lai
            else:
                sim_max = sim[sim['TT']>=max_tt]
                sim_max.plot('LAI_tot', 'TCtot', style='-^'+color, label='Sim '+var, alpha=0.4)
            
            #mise en forme
            #plt.ylim(ymin=0); plt.xlim(xmin=0); plt.ylim(ymax=1); plt.xlim(xmax=5)
            plt.xlabel("LAI total"); plt.ylabel("TC0 total")
            plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
            
    if graph3==True:
        for var in var_lst:
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
                
            #observation
            #TC vert obs
            obs_tc = archidb.TC_data()[var+'_0']
            obs_new = obs_tc.groupby('TT').mean() #mean
            obs_new = obs_new.reset_index()            
            #mat rayonnement data
            tab0 = 'files/ray_obs_sol_'+var+'.csv'
            tab0 = pandas.read_csv(tab0, decimal='.', sep=',')
            
            obs_new['%'] = 0.0; obs_new['% exact'] = 0.0; obs_new['1-% exact'] = 0.0
            n=0
            while n < len(obs_new.index):
                pt = pandas.DataFrame()
                pt = tab0[(tab0['TT'] >= obs_new['TT'][n])]
                pt = pt.reset_index()               
                if pt['TT'][0] == obs_new['TT'][n] :
                    obs_new['% exact'][n] = pt['%'][0]
                    obs_new['1-% exact'] = 1 - obs_new['% exact']
                    obs_new['%'][n] = None; obs_new['1-%'][n] = None
                    n = n + 1
                else :
                    obs_new['%'][n] = (  pt['%'][0] * obs_new['TT'][n] ) / pt['TT'][0]
                    obs_new['1-%'] = 1 - obs_new['%']
                    obs_new['% exact'][n] = None; obs_new['1-% exact'][n] = None
                    n = n + 1   
               
            obs_new.plot('TC', '1-% exact', style='^'+color, label='pt exact '+var, alpha=1)
            obs_new.plot('TC', '1-%', style='o'+color, label='pt extrapole '+var, alpha=1)
            obs_new['1-% all']=pandas.concat([obs_new['1-%'].dropna(), obs_new['1-% exact'].dropna()]).reindex_like(obs_new)
            obs_new.plot('TC', '1-% all', style='--'+color, label='_nolegend_', alpha=0.5)
            
            #simulation
            #TC
            sim_tc = 'files/sim_tc_'+var+'.csv'
            sim_tc = pandas.read_csv(sim_tc, decimal='.', sep=',')
            #meteo
            sim_ray = 'files/ray_sim_sol_'+var+'.csv'
            sim_ray = pandas.read_csv(sim_ray, decimal='.', sep=',')
            #merge
            sim = sim_tc.merge(sim_ray)
            #plot
            sim.plot('TCgreen', '1-light', style='-x'+color, label='simulation '+var, alpha=0.5)
                       
        #mise en forme
        plt.ylim(ymin=0.7); plt.xlim(xmin=0); plt.ylim(ymax=1); plt.xlim(xmax=1)
        plt.xlabel("TC vert"); plt.ylabel("1 - rayonnement obs")
        plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})

    if graph4==True:
        for var in var_lst:
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'

            #observation
            #LAI vert obs
            obs_lai = archidb.PAI_photo_data()[var]    
            obs_lai = obs_lai.rename(columns={'TT_date': 'TT'}) #rename column
            obs_new = obs_lai.groupby('TT').mean() #mean by TT
            obs_new = obs_new.reset_index()
            #mat rayonnement data
            tab0 = 'files/ray_obs_sol_'+var+'.csv'
            tab0 = pandas.read_csv(tab0, decimal='.', sep=',')
            
            obs_new['%'] = 0.0; obs_new['% exact'] = 0.0; obs_new['1-% exact'] = 0.0
            n=0
            while n < len(obs_new.index):
                pt = pandas.DataFrame()
                pt = tab0[(tab0['TT'] >= obs_new['TT'][n])]
                pt = pt.reset_index()               
                if pt['TT'][0] == obs_new['TT'][n] :
                    obs_new['% exact'][n] = pt['%'][0]
                    obs_new['1-% exact'] = 1 - obs_new['% exact']
                    obs_new['%'][n] = None; obs_new['1-%'][n] = None
                    n = n + 1
                else :
                    obs_new['%'][n] = (  pt['%'][0] * obs_new['TT'][n] ) / pt['TT'][0]
                    obs_new['1-%'] = 1 - obs_new['%']
                    obs_new['% exact'][n] = None; obs_new['1-% exact'][n] = None
                    n = n + 1   
                  
            obs_new.plot('PAI_vert_photo', '1-% exact', style='^'+color, label='pt exact '+var, alpha=1)
            obs_new.plot('PAI_vert_photo', '1-%', style='o'+color, label='pt extrapole '+var, alpha=1)
                  
            obs_new['1-% all']=pandas.concat([obs_new['1-%'].dropna(), obs_new['1-% exact'].dropna()]).reindex_like(obs_new)
            obs_new.plot('PAI_vert_photo', '1-% all', style='--'+color, label='_nolegend_', alpha=0.5)
            
            #simulation
            #LAI
            sim_lai = 'files/sim_lai_'+var+'.csv'
            sim_lai = pandas.read_csv(sim_lai, decimal='.', sep=',')
            #meteo
            sim_ray = 'files/ray_sim_sol_'+var+'.csv'
            sim_ray = pandas.read_csv(sim_ray, decimal='.', sep=',')
            #merge
            sim = sim_lai.merge(sim_ray)
            #plot
            sim.plot('LAI_vert', '1-light', style='-x'+color, label='simulation '+var, alpha=0.5)

        #mise en forme
        plt.ylim(ymin=0); plt.xlim(xmin=0); plt.ylim(ymax=1); plt.xlim(xmax=10)
        plt.xlabel("LAI vert"); plt.ylabel("1 - rayonnement obs au sol")
        plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
        
    if graph5==True:
        for var in var_lst:
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            
            #LAI sim           
            adel, domain, domain_area, convUnit, nplants = get_reconstruction(name=var, nplants = 30, aborting_tiller_reduction = 1, seed=1)   
            sim_lai = simLAI(adel, domain_area, convUnit, nplants)
            sim_lai = sim_lai.rename(columns={'ThermalTime' : 'TT'})
            
            #TC sim
            #sim_tc = 'files/sim_tc_'+var+'.csv'
            #sim_tc = pandas.read_csv(sim_tc, decimal='.', sep=',')            
            dd = range(400,2600,100)    
            adel, domain, domain_area, convUnit, nplants = get_reconstruction(name=var, nplants = 30, aborting_tiller_reduction = 1, seed=1)
            sim = (adel.setup_canopy(age) for age in dd)
            TC_sim = [draft_TC(g, adel, domain, zenith=57, rep=2) for g in sim]
            n=0; tc=[]; tc_sen=[]
            while n<len(TC_sim):
                tc.append(TC_sim[n]['green'])
                tc_sen.append(TC_sim[n]['senescent'])
                n = n + 1
            sim_tc = pandas.DataFrame({'TT':dd, 'TCgreen':tc, 'TCsen':tc_sen})
            #sim_tc['TCtot'] = sim_tc['TCgreen'] + sim_tc['TCsen']
            
            #merge LAI sim + TC sim
            sim = sim_tc.merge(sim_lai)
            #decoupage selon max value du LAI
            max_lai = sim['LAI_vert'].max()
            max_tt = int(sim.TT[sim['LAI_vert']==max_lai])
            
            #plot <= max_lai
            if stade=='inf':
                sim_min = sim[sim['TT']<=max_tt]
                sim_min.plot('LAI_vert', 'TCgreen', style='-^'+color, label='Sim '+var, alpha=0.4) 
            #plot >= max_lai
            else:
                sim_max = sim[sim['TT']>=max_tt]
                sim_max.plot('LAI_vert', 'TCgreen', style='-^'+color, label='Sim '+var, alpha=0.4)
            
            #mise en forme
            #plt.ylim(ymin=0); plt.xlim(xmin=0); plt.ylim(ymax=1); plt.xlim(xmax=5)
            plt.xlabel("LAI vert"); plt.ylabel("TC 57 vert")
            plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
        
        return sim
