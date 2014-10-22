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
    adel,domain, domain_area, convUnit, nplants = get_reconstruction(name=name, nplants=nplants)
    
    if sim is 'Mercia_T1':
        HS = 9.74 #HS mesure au 18/04 au lieu du 19/04
    elif sim is 'Mercia_T2':
        HS = 12.8 #HS convert car pas de HS mesure +/- 15j
    elif sim is 'Rht3_T1':
        HS = 9.15 #HS mesure au 18/04 au lieu du 19/04
    elif sim is 'Rht3_T2':
        HS = 12.48 #HS convert car pas de HS mesure +/- 15j
    elif sim is 'Tremie12_T1':
        HS = 10.98 #HS mesure
    elif sim is 'Tremie12_T2':
        HS = 12.63 #HS mesure
    elif sim is 'Tremie13_T1':
        HS = 8.33 #HS mesure au 19/04 au lieu du 22/04 (8.36 avec HS convert)
    elif sim is 'Tremie13_T2':
        HS = 11.04 #HS convert
    else :
        print 'Warning = PB name or sim'
    
    #conversion du HS en TT necessaire a setup_canopy
    conv = HSconv[name]
    age = conv.TT(HS)
    g = adel.setup_canopy(age)
    
    if sim is 'Mercia_T1' or name is 'Rht3_T1':
        df = repartition_at_applicationArch(appdate = '2011-04-19', dose = 1, g=g)
    elif sim is 'Mercia_T2' or name is 'Rht3_T2':
        df = repartition_at_applicationArch(appdate = '2011-05-11', dose = 1, g=g)
    elif sim is 'Tremie12_T1':
        df = repartition_at_applicationArch(appdate = '2012-04-11', dose = 1, g=g)
    elif sim is 'Tremie12_T2':
        df = repartition_at_applicationArch(appdate = '2012-05-09', dose = 1, g=g)
    elif sim is 'Tremie13_T1':
        df = repartition_at_applicationArch(appdate = '2013-04-25', dose = 1, g=g)
    elif sim is 'Tremie13_T2':
        df = repartition_at_applicationArch(appdate = '2013-05-17', dose = 10000, g=g)
    else :
        print 'Warning = PB name or sim'
    return df
        
    # A supprimer si tt fonctionne au dessus
    # stade 2N et DFE pour Mercia et Rht3
    #df = repartition_at_applicationArch(appdate = '2011-04-19', dose = 1, g=g)
    #df = repartition_at_applicationArch(appdate = '2011-05-11', dose = 1, g=g)
    # stade 2N et DFE pour Tremie12
    #df = repartition_at_applicationArch(appdate = '2012-04-11', dose = 1, g=g)
    #df = repartition_at_applicationArch(appdate = '2012-05-09', dose = 1, g=g)
    # stade 2N et DFE pour Tremie13
    #df = repartition_at_applicationArch(appdate = '2013-04-25', dose = 1, g=g)
    #df = repartition_at_applicationArch(appdate = '2013-05-17', dose = 10000, g=g)

    
def treatment(name='Tremie12', sim='T1', nplants=30): 

    simul = name+'_'+sim
    df = g_constr(name = name, sim = simul, nplants = nplants)
    gr=df.groupby(['plant', 'date', 'axe'], group_keys=False)
    
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
        
    data = gr.apply(_fun)
    
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1

    data = data.convert_objects() # conversion des objets en types numpy

    data.to_csv('dataALL.csv')
    data = data[(data.ntop_cur <= 4)&(data.axe=='MS')]
    data.to_csv('data_ntop_cur4.csv')
    
    #Calcul surface foliaire : est egale à la sommme de toutes les surfaces par metamer
    dataSF=data.groupby(['plant', 'date', 'metamer','ntop_cur'],as_index=False)
#    dataSF['Foliar_Surface']=dataSF['green_area']+dataSF['senesced_area']
    dataSFsum=dataSF['area'].agg({'area_sum':np.sum,'area_std':np.std})

#    dataSFsum = dataSF['area'].agg({'area_sum' : np.sum})
    dataSFsum.to_csv('dataSFsum.csv')
#    if data['axe'] is 'MS' : 

#    dataCB=data.groupby(['date'],as_index=False)
    # calcul, pour chaque date de la surface totale verte et senecente et de l'interception totale
#    grCB= dataCB.sum()
#    grCB.to_csv('grCB.csv')
    #grCBm=grCB.groupby(['ntop_cur', 'date', 'axe'],as_index=False)
    #grCFmean = grCBm['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})
    #grCBmean.to_csv('grCBmean.csv')
    

    #calcul pour sortir donnees interception sur leaf 1 et metamer A CONFIRMER
    dataCF=data.groupby(['plant', 'ntop_cur', 'date','metamer'], as_index=False)
    dataCF.to_csv('dataCF.csv')
#    grCFmeanSF = dataCF['Foliar_Surface'].agg({'Foliar_Surface_sum' : np.sum, 'Foliar_Surface_std':np.max})
    grCF= dataCF.sum()
    grCF.to_csv('grCF.csv')
    grCFm=grCF.groupby(['ntop_cur','date'],as_index=False)
    grCFm.to_csv('grCFm.csv')

# essai ajout intervalle de confiance (, 'surfacic_doses_IC':np.(mean+((1.96*np.std)** n))
    grCFmeanSDE = grCFm['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std':np.std})
#    grCFarea= grCFm['area']

    grCFmeanSF = grCFm['area'].agg({'area_mean' : np.mean, 'area_std':np.std})
#    grCFmean = pandas.concat([grCFmeanSDE, grCFmeanSF])
    grCFmean = pandas.merge(grCFmeanSDE,grCFmeanSF)
#    grCFmean = pandas.merge(grCFmeanSDE,grCFarea)
#    grCFmean['SDESF']= grCFmean['surfacic_doses_Epoxiconazole_mean']*dataSFsum['area_sum']
#    grCFmean['SDESF']= grCFmean['surfacic_doses_Epoxiconazole_mean']*grCFm['area']
    grCFmean['SDESF']= grCFmean['surfacic_doses_Epoxiconazole_mean']*grCFmean['area_mean']

#    grCFmean = grCFm.agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std,'Foliar_Surface_mean' : np.mean,'Foliar_Surface_std':np.std})
    grCFmean['name'] = name

#old
    grCFmean.to_csv('grCFmean'+name+'.csv')
    
    global_grCFmean_filename = 'grCFmean_global.csv'
    
    if os.path.isfile(global_grCFmean_filename):
        grCFmean_global = pandas.read_csv(global_grCFmean_filename, index_col=False)
        grCFmean_global = pandas.concat([grCFmean_global, grCFmean], ignore_index=True)
    else:
        grCFmean_global = grCFmean
    
    grCFmean_global = grCFmean_global.to_csv(global_grCFmean_filename, index=False)

    #ficname = list(grCFmean.to_csv('grCFmean'+name+'.csv'))
    #concatenation des fichiers plutot que de les lire/ecrire (pb car beoins de liste or ici traitement de dataframe
    #import shutil
    #fichier_final = open("ficname[0]",'w')
    #for i in ficname.count:
    #      shutil.copyfileobj(open(i, 'r'), fichier_final)
    #fichier_final.close()


    #new ci-dessous a peaufiner
    #
    #import csv
    #c = csv.writer(open('grCFmean.csv',"w"))
    #list= list(
    #n=0
    #while grCFmean.loc[n+1]<grCFmean.len:
    #    c.writerow(grCFmean.loc[n])
    
    #grCFmean.plot('ntop_cur','surfacic_doses_Epoxiconazole_mean' )
    

    #test groupby by multiple apres avoir fait calcul sur data
    """data = data[data.ntop_cur <= 4]
    f = {'surfacic_doses_Epoxiconazole':['sum','mean']}
    dftest = fdf.groupby(['plant', 'ntop_cur', 'date', 'axe']).agg(f)
    """
    """script calcul intervalle de confiance
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

    gr3=dmp.groupby(['date'],as_index=False)
    ## calcul, pour chaque date et ntop_cur<4 de la surface totale verte et senecente et de l'interception totale
    grtest3 = gr3.sum()
       
    return grtest4,grtest3

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
        print type(data.index)
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
	
#grtest = dmp.groupby('ntop_cur','date').apply(lambda x: 
#             Series({'r': (x.y + x.z).sum() / x.z.sum(), 
#                     's': (x.y + x.z ** 2).sum() / x.z.sum()}))

#gr=dmp.groupby(['ntop_cur','date'])
#dms=gr.agg(numpy.mean)


#NG ajout ecart type
#dmserr=gr.agg(numpy.std)
#dms=dms.reset_index()
#dmserr=dmserr.reset_index()

#plt.figure()
#plt = grtest.plot('ntop_cur','surfacic_doses_Epoxiconazole_mean',kind="bar")
#plt.show()

'''mise en forme sous format barplot
fig = p.figure()
ax = fig.add_subplot(1,1,1)
#N = len(ntop_cur)

#err = [1.2, 1.5, 2.5, 1.2]
#ax.bar('ntop_cur', 'n_max', facecolor='#777777',
#       align='center', yerr=err, ecolor='black')
ax.bar(dms.ntop_cur, dms.n_max, facecolor='#777777',
       align='center', ecolor='black')
#Create a y label
ax.set_ylabel('Intercepted dose (?/?)')
 
# Create a title, in italics
ax.set_title('Intercepted dose, by ntop_cur',fontstyle='italic')
 
# This sets the ticks on the x axis to be exactly where we put
# the center of the bars.
ax.set_xticks(dms.ntop_cur)
 
# Labels for the ticks on the x axis.  It needs to be the same length
# as y (one label for each bar)
group_labels = ['1', '2',
                 '3', '4']
 
# Set the x tick labels to the group_labels defined above.
ax.set_xticklabels(group_labels)
 
# Extremely nice function to auto-rotate the x axis labels.
# It was made for dates (hence the name) but it works
# for any long x tick labels
fig.autofmt_xdate()
 
p.show()

#
# OK autre solution pour les barplots
#

#definition du nombre de groupe : correspond ici à la dose de pesticide epandu
n_groups = 3

fig, ax = plt.subplots()

#index = np.arange(n_groups)
bar_width = 1

opacity = 0.4
error_config = {'ecolor': '0.3'}

#for i in range (0,n_groups):
rects1 = plt.bar(grtest2.ntop_cur[0], grtest2.surfacic_doses_Epoxiconazole_mean[0], bar_width,
                 alpha=opacity,
                 color='b',
                 yerr=grtest2.surfacic_doses_Epoxiconazole_std[0],
                 error_kw=error_config,
                 label='leaf 1')




                 
plt.xlabel('Dose 0.5')
plt.ylabel('Quantite retenue')
plt.title('Quantite de pesticide par feuille')
plt.xticks(grtest2.ntop_cur, (' '))
plt.legend()

plt.tight_layout()
plt.show()

#solution 3 : probleme avec non reconnaissance DataFrame

#df2 = DataFrame(dms.n_max[1], columns=['1', '2', '3', '4'])
plt.figure();
#grplot = DataFrame(grtest2.surfacic_doses_Epoxiconazole_mean,columns=['1','2','3','4'])
grplot = grtest2.surfacic_doses_Epoxiconazole_mean
grplot.plot(kind='bar',yerr=grtest2.surfacic_doses_Epoxiconazole_std,grid=False)
plt.show()

"""autre solution"""
fig = plt.figure()
ax1 = plt.subplot(111,ylabel='A')
#ax2 = gcf().add_axes(ax1.get_position(), sharex=ax1, frameon=False, ylabel='axes2')
#ax2 =ax1.twinx()
#ax2.set_ylabel('B')
ax1.bar(grtest2.ntop_cur.values,grtest2.surfacic_doses_Epoxiconazole_mean.values, width =0.4, color ='g', align = 'center')
#ax2.bar(drtest2.index,df.B.values, width = 0.4, color='r', align = 'edge')
ax1.legend(['A'], loc = 'upper left')
#ax2.legend(['B'], loc = 'upper right')
fig.show()

N = 1
ind = np.arange(N)  # the x locations for the groups dans echap = nombre de doses = N
width = 0.2    # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)

yvals = [grtest2.surfacic_doses_Epoxiconazole_mean[0]]
rects1 = ax.bar(ind, yvals, width, color='r')
zvals = [grtest2.surfacic_doses_Epoxiconazole_mean[1]]
rects2 = ax.bar(ind+width, zvals, width, color='g')
kvals = [grtest2.surfacic_doses_Epoxiconazole_mean[2]]
rects3 = ax.bar(ind+width*2, kvals, width, color='b')
mvals =  [grtest2.surfacic_doses_Epoxiconazole_mean[3]]
rects4 = ax.bar(ind+width*3, mvals, width, color='y',yerr=grtest2.surfacic_doses_Epoxiconazole_std[3])

ax.set_ylabel('Depot pesticides par etage foliaire')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('etage 1', 'etage 2', 'etage 3','etage 4') )
#ax.legend( (rects1[0], rects2[0], rects3[0],rects4[0]), ('1', '2', '3','4') )
ax.legend( (yvals, zvals, kvals,mvals), ('1', '2', '3','4') )

def autolabel(rects):
    for rect in rects:
        h = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)

plt.show()

'''
