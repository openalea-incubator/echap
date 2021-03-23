import os
import pandas 
from new_annual_loop_pest import pesticide_loop, repartition_at_application,repartition_at_applicationArch,get_reconstruction
import pesticide_protocol as proto
import numpy as np
import pylab as p
import matplotlib.pyplot as plt
from alinea.echap.weather_data import *
import time 
from alinea.echap.recorder import *
import alinea.echap.architectural_reconstructions as reconstructions

from itertools import cycle, islice
from pylab import *

plt.ion()


# from openalea.deploy.shared_data import shared_data
# from alinea.astk.Weather import Weather
# from alinea.astk.TimeControl import *
# from alinea.caribu.caribu_star import rain_and_light_star
# from alinea.alep.protocol import infect
# from alinea.echap.interception_leaf import pesticide_applications
# from alinea.echap.microclimate_leaf import microclimate_leaf
# from alinea.echap.recorder import LeafElementRecorder
# from alinea.echap.interfaces import record as do_record #to avoid confusion with numpy record

# from macros_annual_loop import setup_canopy, update_lesions, update_pesticides, dispersion, contamination, pesticide_intercept
# import alinea.septo3d

# from alinea.echap.tests_nodes import plot_pesticide

#g, recorder = pesticide_loop(meteo_file='meteo00-01.txt', applications=proto.Mercia2011_treat1, start="2011-04-11", periods=24, freq='H', TB=0, delayH=1, delayDD=15)
# recorder.save_records('test.csv')   
# recorder.save_plot(what='area',t = 'dd',  by='ntop', axe = 'MS', plant='all', fig_name='test.png')   


# test
"""
g, recorder = pesticide_loop(meteo_file='meteo00-01.txt', start="2001-04-25", periods=50, freq='H', TB=0, delayH=1, delayDD=15, applications=proto.Test2001_treat1)
recorder.save_records('testMAI2014.csv')   
recorder.save_plot(what='surfacic_doses_Epoxiconazole', t = 'iter',  by='ntop', axe = 'MS', plant='all', fig_name='test2.png')   
recorder.plt.show()
"""

# Appel a la fonction repartition_at_application
# df = repartition_at_application(appdate = '2011-04-19', dose = 0.5, age = 1166)
#print 'df.columns after = ', df.columns

# -------------------------------------------------
# VERSION SIMPLE (prend un HS en entree)
HSconv = reconstructions.HS_converter

def g_constr(name, sim, nplants):
    # Appel a la fonction get_reconstruction
    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name=name, nplants=nplants)
    
    if sim == 'Mercia_T1':
        HS = 9.74 #HS mesure au 18/04 au lieu du 19/04
    elif sim == 'Mercia_T2':
        HS = 12.8 #HS convert car pas de HS mesure +/- 15j
    elif sim == 'Rht3_T1':
        HS = 9.15 #HS mesure au 18/04 au lieu du 19/04
    elif sim == 'Rht3_T2':
        HS = 12.48 #HS convert car pas de HS mesure +/- 15j
    elif sim == 'Tremie12_T1':
        HS = 10.98 #HS mesure
    elif sim == 'Tremie12_T2':
        HS = 12.63 #HS mesure
    elif sim == 'Tremie13_T1':
        # 22/04
        # HS = 8.36 #HS convert (HS mesure 8.33 au 19/04 au lieu du 22/04)
        # 25/04
        HS = 8.7
    elif sim == 'Tremie13_int':
        HS = 9.7 #HS mesure au 02/05 au lieu du 03/05 (date intermediaire pour surfaces Steph)
    elif sim == 'Tremie13_T2':
        HS = 11.04 #HS convert
    else :
        print('Warning = PB sim - pas de HS')
    
    #conversion du HS en TT necessaire a setup_canopy
    conv = HSconv[name]
    age = conv.TT(HS)
    g = adel.setup_canopy(age)
    
    if sim == 'Mercia_T1' or sim == 'Rht3_T1':
        df = repartition_at_applicationArch(appdate = '2011-04-19', dose = 1, g=g)
    elif sim == 'Mercia_T2' or sim == 'Rht3_T2':
        df = repartition_at_applicationArch(appdate = '2011-05-11', dose = 1, g=g)
    elif sim == 'Tremie12_T1':
        df = repartition_at_applicationArch(appdate = '2012-04-11', dose = 1, g=g)
    elif sim == 'Tremie12_T2':
        df = repartition_at_applicationArch(appdate = '2012-05-09', dose = 1, g=g)
    elif sim == 'Tremie13_T1':
        df = repartition_at_applicationArch(appdate = '2013-04-25', dose = 1, g=g)
    elif sim == 'Tremie13_int':
        df = repartition_at_applicationArch(appdate = '2013-05-03', dose = 1, g=g)
    elif sim == 'Tremie13_T2':
        df = repartition_at_applicationArch(appdate = '2013-05-17', dose = 1, g=g)
    else :
        print('Warning = PB sim - pas de return df')
        
    return df
    
#script calcul intervalle de confiance
def mean(lst):
    return sum(lst) / float(len(lst))

def variance(lst):
    """
    #Uses standard variance formula (sum of each (data point - mean) squared)
    #all divided by number of data points
    """
    mu = mean(lst)
    return 1.0/(len(lst)-1) * sum([(i-mu)**2 for i in lst])

def conf_int(lst, perc_conf=95):
    """
    #Confidence interval - given a list of values compute the square root of
    #the variance of the list (v) divided by the number of entries (n)
    #multiplied by a constant factor of (c). This means that I can
    #be confident of a result +/- this amount from the mean.
    #The constant factor can be looked up from a table, for 95pcent confidence
    #on a reasonable size sample (>=500) 1.96 is used.
    """
    from scipy.stats import t
    
    n, v = len(lst), variance(lst)
    c = t.interval(perc_conf * 1.0 / 100, n-1)[1]
    
    return math.sqrt(v/n) * c
    
    
def save_csv(name, data, data2, dataSFsum, dataCF, grCF, grCFm, grCFmean):
        data.to_csv('dataALL.csv')
        data2.to_csv('data_ntop_cur4.csv')
        dataSFsum.to_csv('dataSFsum.csv')
        dataCF.to_csv('dataCF.csv')
        grCF.to_csv('grCF.csv')
        grCFm.to_csv('grCFm.csv')
        grCFmean.to_csv('grCFmean'+name+'.csv')
        if os.path.isfile('grCFmean_global.csv'):
            grCFmean_global = pandas.read_csv('grCFmean_global.csv', index_col=False)
            grCFmean_global = pandas.concat([grCFmean_global, grCFmean], ignore_index=True)
        else:
            grCFmean_global = grCFmean
        grCFmean_global.to_csv('grCFmean_global.csv', index=False)

    
def treatment(name='Tremie12', sim='T1', nplants=30, save=True): 

    simul = name+'_'+sim
    df = g_constr(name, simul, nplants)
    gr=df.groupby(['plant', 'date', 'axe'], group_keys=False)
    
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub        
    data = gr.apply(_fun)
    
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1
    data = data.convert_objects() # conversion des objets en types numpy
    data2 = data[(data.ntop_cur <= 4)&(data.axe=='MS')]

    #Calcul surface foliaire : est egale à la sommme de toutes les surfaces par metamer
    dataSF=data2.groupby(['plant', 'date', 'metamer','ntop_cur'],as_index=False)
    dataSFsum=dataSF['area'].agg({'area_sum':np.sum,'area_std':np.std})

    #calcul pour sortir donnees interception sur leaf 1 et metamer A CONFIRMER
    dataCF = data2.groupby(['plant', 'ntop_cur', 'date','metamer'], as_index=False)
    grCF = dataCF.sum()
    grCFm = grCF.groupby(['ntop_cur','date'],as_index=False)
    
    grCFmeanSDE = grCFm['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std':np.std})

    grCFmeanSF = grCFm['area'].agg({'area_mean' : np.mean, 'area_std':np.std})
    grCFmean = pandas.merge(grCFmeanSDE,grCFmeanSF)
    grCFmean['SDESF']= grCFmean['surfacic_doses_Epoxiconazole_mean']*grCFmean['area_mean']

    grCFmean['name'] = name

    data3 = data[data.ntop_cur <= 4]
    gr=data3.groupby(['plant', 'ntop_cur', 'date', 'axe'], as_index=False)
    dmp=gr.sum() # dmp contient, pour chaque quadruplet ('plant', 'ntop_cur', 'date', 'axe'), 
                 # la somme des area, green_area, id, length, metamer, ntop, penetrated_doses_Epoxiconazole, 
                 # senesced_area, surfacic_doses_Epoxiconazole, n_max

    gr = dmp.groupby(['ntop_cur','plant','date'], as_index=False)
    # calcul, pour chaque triplet ('ntop_cur','plant','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest = gr['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})

    gr2 = dmp.groupby(['ntop_cur','date'],as_index=False)
    # calcul, pour chaque doublon ('ntop_cur','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest4 = gr2['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})

    gr3=dmp.groupby(['date'],as_index=False)
    ## calcul, pour chaque date et ntop_cur<4 de la surface totale verte et senecente et de l'interception totale
    grtest3 = gr3.sum()
          
    if save is True:
        save_csv(name, data, data2, dataSFsum, dataCF, grCF, grCFm, grCFmean)
    else:
        print('Pas de sauvegarde de fichiers sous forme csv')
  
    return grtest4, grtest3

#plot des donnees obs
def plot(name='Mercia'):

    if name is 'Mercia' or 'Rht3' :
#        data_file = 'ObsData2011.csv'
         data_file = 'ObsdataFinal.csv'
    df_obs = pandas.read_csv(data_file, sep=';')
    #df_obs = df_obs[df_obs['feuille']<5]
    df_obs['var_stade'] = df_obs['name'] + "_" + df_obs['date']

    cols = df_obs.var_stade.value_counts().shape[0]
    fig, axes = plt.subplots(1, cols, figsize=(8, 8))
    
    # avoir les graphes par ordre alphabetique
    val = df_obs.var_stade.value_counts().index.values
    val.sort()

    for x, var_stade in enumerate(val):
        dat = df_obs[(df_obs['var_stade'] == var_stade)]
        data = dat.groupby(['feuille']).moy_intercepte.sum()
        print (data)
        print(type(data.index))
        left = [k[0] for k in enumerate(data)]
        right = [k[1] for k in enumerate(data)]

        #gestion des 4 couleurs differentes des barplot
        my_colors = list(islice(cycle(['#8181F7', '#FA5858', '#BEF781', '#F3F781','#765081']), None, 5))
        SD = dat.groupby(['feuille']).ecart_type_intercepte.sum()
        axes[x].bar(left,right,label="%s" % (var_stade), width=0.5, color=my_colors, alpha=0.4, yerr=SD)
        axes[x].set_xticks(left, minor=False)
        # modification de la legende axe des abcisses
        data = data.reset_index()
        #axes[x].set_xticklabels(data['dose'].values)
        labels = [item.get_text() for item in axes[x].get_xticklabels()]
#        labels[1] = '|   40'; labels[5] = '|   80'; labels[9] = '|   150'; labels[13] = '|   200'
        axes[x].set_xticklabels(labels)
        #axe des ordonnees
        axes[x].set_ylim(0, 35)
        axes[x].grid(True)
        #changement couleur autour graph
        fig.patch.set_facecolor('#FCF8F8')
        #changement couleur fond graph
        axes[x].patch.set_facecolor('#D1D1D1')
        #titre
        fig.suptitle('Dose par Feuille par Stade par Variete', fontsize=15)
        axes[x].text(0, 23, var_stade, bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10})

    axes[x].text(6, 20, '1', style='italic', bbox={'facecolor':'#8181F7', 'alpha':0.6, 'pad':10})
    axes[x].text(6, 18.5, '2', style='italic', bbox={'facecolor':'#FA5858', 'alpha':0.6, 'pad':10})
    axes[x].text(6, 17, '3', style='italic', bbox={'facecolor':'#BEF781', 'alpha':0.6, 'pad':10})
    axes[x].text(6, 15.5, '4', style='italic', bbox={'facecolor':'#F3F781', 'alpha':0.6, 'pad':10})
    axes[x].text(6, 14, '5', style='italic', bbox={'facecolor':'#765081', 'alpha':0.6, 'pad':10})
    plt.show()

'''
# -----------------------------------------------------------------------------------------------------------------------------------
# IDEM mais gestion d'un range de date  
def g_constr_age(name='Mercia',nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, dose=0.5, age=1166):

    # Conversion age (TT) en date ('aaaa-mm-jj')
    from alinea.astk.TimeControl import thermal_time
    if name is 'Mercia' or 'Rht3':
        data_file = 'METEO_stationINRA_20102011.csv'
        # correspondance entre date et TT
        met = Boigneville_2010_2011()
        seq = pandas.date_range(start="2010-10-15", end="2011-06-20", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2010-10-15", end="2011-06-20")
        bid = bid[seq]
        #renommer colonnes
        bid = bid.reset_index()
        bid.columns = ['datetime','TT']
        bid['datetime'] = bid['datetime'].map(lambda x: x.strftime('%Y-%m-%d'))
        bid['TT'] = bid['TT'].astype(int)
        # trouver la correspondance en notre age TT et la date dont a besoin repartition_at_applicationArch()
        da=0; appdate=0
        while appdate==0:
            if da<len(bid):
                if bid['TT'][da]==age:
                    appdate=bid['datetime'][da]
                else:
                    da = da+1
            if da==len(bid):
                da=0; age=age+1
    # Appel a la fonction get_reconstruction
    pgen,adel,domain, domain_area, convUnit, nplants = get_reconstruction(name=name, nplants=nplants, nsect=nsect, seed=seed, sample=sample, as_pgen=as_pgen, dTT_stop=dTT_stop)
    # Appel a la fonction repartition_at_application lié à l'architecture
    g = adel.setup_canopy(age)
    df = repartition_at_applicationArch(appdate = appdate, dose = 40, g = g)
    return df

def treatment_age(name='Mercia',nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, dose=0.5): 

    dd = range(1000,1400,300)
    #df_sim = [g_constr_age(name, nplants, nsect, seed, sample, as_pgen, dTT_stop, dose, age) for age in dd]
    
    # lancer la reconstruction pour chaque age du range et concatener les dataframe de sortie
    dl = len(dd); dlc = 1
    df = g_constr_age(name, nplants, nsect, seed, sample, as_pgen, dTT_stop, dose, dd[0])
    while dlc < dl :
        df_other = g_constr_age(name, nplants, nsect, seed, sample, as_pgen, dTT_stop, dose, dd[dlc])
        dlc = dlc+1
        df = df.append(df_other, ignore_index=True)

    gr = df.groupby(['plant', 'date', 'axe'], group_keys=False)
    
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
    
    data = gr.apply(_fun)    
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1    
    data = data.convert_objects() # conversion des objets en types numpy
    
    #test groupby by multiple apres avoir fait calcul sur data
    """data = data[data.ntop_cur <= 4]
    f = {'surfacic_doses_Epoxiconazole':['sum','mean']}
    dftest = fdf.groupby(['plant', 'ntop_cur', 'date', 'axe']).agg(f)
    """

    data = data[data.ntop_cur <= 4]
    gr=data.groupby(['plant', 'ntop_cur', 'date', 'axe'], as_index=False)
    dmp=gr.sum() # dmp contient, pour chaque quadruplet ('plant', 'ntop_cur', 'date', 'axe'), 
                 # la somme des area, green_area, id, length, metamer, ntop, penetrated_doses_Epoxiconazole, 
                 # senesced_area, surfacic_doses_Epoxiconazole, n_max
     
    gr = dmp.groupby(['ntop_cur','plant','date'],as_index=False)
    # calcul, pour chaque triplet ('ntop_cur','plant','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest = gr['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})

    gr2 = dmp.groupby(['ntop_cur','date'],as_index=False)
    # calcul, pour chaque doublon ('ntop_cur','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest4 = gr2['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})
    
    grtest4.to_csv('grtest_treatment_age')

    return grtest4   											 
'''