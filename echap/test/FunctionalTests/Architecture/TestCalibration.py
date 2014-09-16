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
    pgen, adel, domain, domain_area, convUnit, nplants = fun(**args)
    return pgen, adel, domain, domain_area, convUnit, nplants

def get_pgen(name='Mercia', original = False, dTT_stop = 0):
    fun = reconst_db[name]
    pgen, _, _, _, _, _ = fun(nplants=1,nsect=1,as_pgen=original, dTT_stop=dTT_stop)
    return pgen

# Creation des 2 fichiers pris en entree par curve()
def prep_curve(name='Tremie20122013'):
    
    if name is 'Tremie20112012':
        data_file_xydb = shared_data(alinea.echap, 'prepa_xydb_GrignonTremie2011.csv')
        data_file_srdb = shared_data(alinea.echap, 'prepa_srdb_GrignonTremie2011.csv')
        #out_xy_csv = "C://Users//Administrateur//openaleapkg//echap//src//alinea//echap//share//data" + "/" + 'xydb_GrignonTremie2011.csv'
        #out_sr_csv = "C://Users//Administrateur//openaleapkg//echap//src//alinea//echap//share//data" + "/" + 'srdb_GrignonTremie2011.csv'
        header_row_xydb = ['Image','ID_Plant','ID_Axis','Organ','ID_Metamer','HS','XY','Pt1','Pt2','Pt3','Pt4','Pt5','Pt6','Pt7','Pt8','Pt9','Pt10','Pt11','Pt12','Pt13','Pt14','Pt15']									
        dfxy = pandas.read_csv(data_file_xydb, names=header_row_xydb, sep=';', index_col=0, skiprows=1, decimal=',')
        header_row_srdb = ['operateur','MethodeAnalyse','Date d\'analyse','fichier','Var','N_prelev','date prelev','rep','N_plante','id_Axe','id_Feuille','nmax','HS','Llimbe','Lgaine','Lentrenoeud','%vert','Llimbe2','Wlimbe','A_bl','A_bl_green','statut','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00']									
        dfsr = pandas.read_csv(data_file_srdb, names=header_row_srdb, sep=';', index_col=0, skiprows=1, decimal=',')
    if name is 'Tremie20122013':
        data_file_xydb = shared_data(alinea.echap, 'prepa_xydb_GrignonTremie2012.csv')
        data_file_srdb = shared_data(alinea.echap, 'prepa_srdb_GrignonTremie2011.csv')
        #out_xy_csv = "C://Users//Administrateur//openaleapkg//echap//src//alinea//echap//share//data" + "/" + 'xydb_GrignonTremie2012.csv'
        #out_sr_csv = "C://Users//Administrateur//openaleapkg//echap//src//alinea//echap//share//data" + "/" + 'srdb_GrignonTremie2012.csv'
        header_row_xydb = ['Image','ID_Plant','Stem','Organ_type','Organ','HS','XY','Pt1','Pt2','Pt3','Pt4','Pt5','Pt6','Pt7','Pt8','Pt9','Pt10','Pt11','Pt12','Pt13','Pt14','Pt15','Pt16','Pt17','Pt18']									
        dfxy = pandas.read_csv(data_file_xydb, names=header_row_xydb, sep=';', index_col=0, skiprows=1, decimal=',')
        header_row_srdb = ['operateur','MethodeAnalyse','Date d\'analyse','fichier','Var','N_prelev','date prelev','rep','N_plante','id_Axe','id_Feuille','nmax','HS','Llimbe','Lgaine','Lentrenoeud','%vert','Llimbe2','Wlimbe','A_bl','A_bl_green','statut','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00']									
        dfsr = pandas.read_csv(data_file_srdb, names=header_row_srdb, sep=';', index_col=0, skiprows=1, decimal=',')
    
    # preparation fichier pr methode R Christian
    dfxy = dfxy.reset_index()
    #traitement des x et y
    #x=[]; y=[]; fin=[]; n = 0; row_tot = dfxy['XY'].count(); chif = range(7,21,1)
    x=[]; y=[]; fin=[]; n = 0; row_tot = dfxy['XY'].count(); chif = range(7,24,1)
    
    c = csv.writer(open("test.csv", "wb"))
    c.writerow(["ID_plant","Stem","Organ_type","Organ","HS","stat_curve","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19","x20","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"])
    
    for row in dfxy.loc[n]:
        while n < row_tot :
            if dfxy['XY'][n]==0 :
                # initialisation liste
                x=[]; y=[]; fin=[]
                # ajout valeurs premieres colonnes
                plant = (dfxy.loc[n,1]); fin.append(plant)
                axis = (dfxy.loc[n,2]); fin.append(axis)
                organ = (dfxy.loc[n,3]); fin.append(organ)
                metamer = (dfxy.loc[n,4]); fin.append(metamer)
                HS = (dfxy.loc[n,5]); fin.append(HS)
                stat = 3; fin.append(stat)
                # remplissage de x et y
                for e in chif:
                    x.append(dfxy.loc[n,e])
                    y.append(dfxy.loc[n+1,e])
                    # enlever valeur nan
                    xe = [value for value in x if not math.isnan(value)]
                    ye = [value for value in y if not math.isnan(value)]
                # obtenir 20 valeurs de x en allant de xmin a xmax
                xmin = min(xe); xmax = max(xe)
                pdt = (xmax-xmin)/19
                if pdt == 0:
                    # !!! CAS PARTICULIER, il y a tjrs au minimun 2 valeurs mais il se peut qu'elles soient identiques > pdt=0 donc pb dans le arange !!!
                    newx = [xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin,xmin]
                else:
                    newx = numpy.arange(xmin, xmax+pdt, pdt)
                # interpolation de y
                interpy = numpy.interp(newx, xe, ye)
                # concatenation x et y pour ecriture dans csv
                ran = range(0,20,1)
                for r in ran:
                    fin.append(newx[r])
                for r in ran:
                    fin.append(interpy[r])
                # ecriture dans csv de sortie
                c.writerow(fin)
                # passage a la ligne suivante
                n = n + 1
            else:
                n = n + 1
                
    '''
    # TRAITEMENT X ET Y
    # obtenir un fichier de la forme :
    # plant,rank,ranktop,HS,x,y
    # 1,7,13,7.6,1.32,5.82

    dfxy = dfxy.reset_index()
    dfxy = dfxy[dfxy['ID_Metamer']>0]
    dfxy['rank'] = dfxy['ID_Metamer']
    #dfxy['ranktop']
    dfxy = dfxy.reset_index()
    dfxy['ranktop']=13
    s_list=['ID_Plant','rank','ranktop','HS','XY','Pt1','Pt2','Pt3','Pt4','Pt5','Pt6','Pt7','Pt8','Pt9','Pt10','Pt11','Pt12','Pt13','Pt14','Pt15']
    dfxy = dfxy[s_list]
    
    plant=[]; rank=[]; HS=[]; x=[]; y=[]; ranktop=[]; n = 0; e = 5; row_tot = dfxy['rank'].count(); col_tot = 20
    
    for row in dfxy.loc[n]:
        while n < row_tot :
            if dfxy['XY'][n]==0 :
                if pandas.isnull(dfxy.loc[n,e]) == False :
                    # recup valeur de la colonne rank et HS
                    plant.append(dfxy.loc[n,0])
                    rank.append(dfxy.loc[n,1])
                    ranktop.append(dfxy.loc[n,2])
                    HS.append(dfxy.loc[n,3])
                    # remplissage de x et y
                    x.append(dfxy.loc[n,e])
                    y.append(dfxy.loc[n+1,e])
                    e = e + 1
                else:
                    n = n + 2
                    e = 5
            if e==col_tot :
                e = 5
                n = n + 2
    plant = map(int, plant); rank = map(int, rank); ranktop = map(int, ranktop)
    DataXY = zip(plant,rank,ranktop,HS,x,y)
    dfxy_fin = pandas.DataFrame(data = DataXY, columns=['plant', 'rank', 'ranktop', 'HS', 'x', 'y'])
    
    # TRAITEMENT S ET R 
    # on ne garde que si id_Axe='MB' et statut=1
    sortdfsr = dfsr[dfsr['id_Axe']=='MB']
    sortdfsr = sortdfsr[sortdfsr['statut']==1]
    # tableau simplifie
    s_list=['id_Feuille','N_plante','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00']
    dfsr = sortdfsr[s_list]
    # enlever les lignes avec des nan
    dfsr = dfsr[dfsr['N_plante'].notnull()]
    dfsr = dfsr[dfsr['id_Feuille'].notnull()]
    dfsr['N_plante'] = dfsr['N_plante'].astype(int); dfsr['id_Feuille'] = dfsr['id_Feuille'].astype(int)
    dfsr['ranktop'] = dfsr['N_plante'] - dfsr['id_Feuille'] + 1
    #dfsr = dfsr.groupby('ranktop')
    # cle = rankclass
    dfsr['rankclass'] = 1
    dfsr['rankclass'][dfsr['ranktop'] <= 4] = 2
    dfsr = dfsr.groupby('rankclass').mean()
    dfsr = dfsr.reset_index()
    s_list=['rankclass','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00']
    dfsr = dfsr[s_list]

    rankclass=[]; s=[]; r=[]; n = 0; e = 1; row_tot = dfsr['rankclass'].count(); col_tot = 22
    for row in dfsr.loc[n]:
        while n < row_tot :
            while e < col_tot :
                rankclass.append(dfsr.loc[n,0])
                s.append(s_list[e])
                r.append(dfsr.loc[n,e])
                e = e + 1
            if e==col_tot:
                e = 1
            n = n + 1
    rankclass = map(int, rankclass)
    DataSR = zip(rankclass,s,r)
    dfsr_fin = pandas.DataFrame(data = DataSR, columns=['rankclass', 's', 'r'])
    
    # export xy et sr en format csv
    #dfxy_fin.to_csv(out_xy_csv, index=False)
    #dfsr_fin.to_csv(out_sr_csv, index=False)'''

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
    
def test_axis_dynamics(name='Mercia'):
# todo : add comparison to primary emission

    pgen = get_pgen(name, original=True)
    newpgen = get_pgen(name, original=False)
    
    if name is 'Mercia':
        obs = archidb.Plot_data_Mercia_Rht3_2010_2011()[name]
    elif name is 'Rht3':
        obs = archidb.Plot_data_Mercia_Rht3_2010_2011()[name]
    elif name is 'Tremie':
        obs = archidb.Plot_data_Tremie_2011_2012()[name]
    elif name is 'Tremie2':
        obs = archidb.Plot_data_Tremie_2012_2013()[name]
    
    def _getfit(pgen):
        primary_proba = pgen['decide_child_axis_probabilities']
        plant_density = pgen['plants_density']
        ears_density = pgen['ears_density']
        ears_per_plant = float(ears_density) / plant_density
        v=list(pgen['MS_leaves_number_probabilities'].values())
        k=list(pgen['MS_leaves_number_probabilities'].keys())
        nff = int(k[v.index(max(v))]) # key of maximal value
        m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
        fit = m.axis_dynamics(plant_density = plant_density)
        return fit
    
    #fit by Mariem
    fit = _getfit(pgen)
    #newfit
    new_fit = _getfit(newpgen)

    fit.plot('HS', ['total','primary','others'], style=['--r','--g','--b'])
    new_fit.plot('HS', ['total','primary','others'],style=['-r','-g','-b'])
    plt.plot([1,13],[obs['plant_density_at_emergence'],obs['ear_density_at_harvest']] ,'ob') 
    
    return fit, new_fit

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

def compare_LAI(name='Mercia_maq1', name_obs='Mercia', dTT_stop=0, original=False, n=30, **kwds): #name='Mercia_comp'

    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original, **kwds)
    sim = simLAI(adel, domain_area, convUnit, nplants)
    sim['HS'] = (sim.ThermalTime - pgen['dynT_user'].TT_col_0[0]) * pgen['dynT_user'].a_cohort[0]
    sim.plot('HS','LAI_vert',color='r')
    
    #donnees obs
    obs = archidb.PAI_data()[name_obs]
    d = {'TT':pandas.Series(obs['TT_date']), 'PAI_vert':pandas.Series(obs['PAI_vert_moy_photo'])}
    df = pandas.DataFrame(d)
    df['HS'] = (df.TT - pgen['dynT_user'].TT_col_0[0]) * pgen['dynT_user'].a_cohort[0]
    df.plot('HS','PAI_vert',style='or'); df.plot('HS','PAI_vert',style='--r')
    
    #Graph LAI_vert et LAI_tot en fonction des TT pour Corinne
    '''sim.plot('ThermalTime','LAI_vert', color='g')
    sim.plot('ThermalTime','LAI_tot', color='r')
    plt.title('Mercia 2010/2011')
    plt.legend(("LAI_vert", "LAI_tot"), 'best')
    plt.show()'''
    
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
    gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = 4288, image_height = 2848, getImages=True, replicate=rep)
    
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

