import os
import pandas
import numpy
from itertools import cycle, islice

import alinea.echap.architectural_reconstructions as reconstructions
from alinea.echap.interception_leaf import InterceptModel, pesticide_applications
from alinea.echap.interfaces import pesticide_interception
import alinea.echap.interception_data as idata

from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record #to avoid confusion with numpy record

import matplotlib.pyplot as plt
plt.ion()

# populate with reconstructions and data
HSconv = reconstructions.HS_converter
Reconstructions = reconstructions.EchapReconstructions()
obs = idata.dye_interception()
HS_applications = idata.dye_applications()
interceptor = InterceptModel({'Tartrazine':{'Tartrazine': 1}}) # consider a 1g.l-1 tartrazine solution

def repartition_at_application(name, appdate='T1', dose=1e4, nplants=30):
    """ 10000 l ha-1 is for 1 l/m2
    """
    
    adel = Reconstructions.get_reconstruction(name=name, nplants=nplants)
    date, hs = HS_applications[name][appdate]
    conv = HSconv[name] 
    age = conv.TT(hs)
    g = adel.setup_canopy(age)
    
    recorder = LeafElementRecorder()
    applications= 'date,dose, product_name\n%s 10:00:00, %f, Tartrazine'%(date, dose)
    application_data = pesticide_applications(applications)
    g,_ = pesticide_interception(g, interceptor, application_data)
    do_record(g, application_data, recorder)
    df =  recorder.get_records()
    print 'repartition_at_application df.columns before ', df.columns
    return df


    
def treatment(name='Tremie12', sim='T1', nplants=30): 

    df = repartition_at_application(name,sim,nplants=nplants)
    gr=df.groupby(['plant', 'date', 'axe'], group_keys=False)
    
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
        
    data = gr.apply(_fun)
    
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1

    data = data.convert_objects() # conversion des objets en types numpy

    data.to_csv('dataALL.csv')
    data = data[(data.ntop_cur <= 5)&(data.axe=='MS')]
    data.to_csv('data_ntop_cur4.csv')
    
    #Calcul surface foliaire : est egale à la sommme de toutes les surfaces par metamer
    dataSF=data.groupby(['plant', 'date', 'metamer','ntop_cur'],as_index=False)
#    dataSF['Foliar_Surface']=dataSF['green_area']+dataSF['senesced_area']
    dataSFsum=dataSF['area'].agg({'area_sum':numpy.sum,'area_std':numpy.std})

#    dataSFsum = dataSF['area'].agg({'area_sum' : numpy.sum})
    dataSFsum.to_csv('dataSFsum.csv')
#    if data['axe'] is 'MS' : 

#    dataCB=data.groupby(['date'],as_index=False)
    # calcul, pour chaque date de la surface totale verte et senecente et de l'interception totale
#    grCB= dataCB.sum()
#    grCB.to_csv('grCB.csv')
    #grCBm=grCB.groupby(['ntop_cur', 'date', 'axe'],as_index=False)
    #grCFmean = grCBm['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : numpy.mean, 'surfacic_doses_Epoxiconazole_std' : numpy.std})
    #grCBmean.to_csv('grCBmean.csv')
    

    #calcul pour sortir donnees interception sur leaf 1 et metamer A CONFIRMER
    dataCF=data.groupby(['plant', 'ntop_cur', 'date','metamer'], as_index=False)
    dataCF.to_csv('dataCF.csv')
#    grCFmeanSF = dataCF['Foliar_Surface'].agg({'Foliar_Surface_sum' : numpy.sum, 'Foliar_Surface_std':numpy.max})
    grCF= dataCF.sum()
    grCF.to_csv('grCF.csv')
    grCFm=grCF.groupby(['ntop_cur','date'],as_index=False)
    grCFm.to_csv('grCFm.csv')

# essai ajout intervalle de confiance (, 'surfacic_doses_IC':numpy.(mean+((1.96*numpy.std)** n))
    grCFmeanSDE = grCFm['surfacic_doses_Tartrazine'].agg({'surfacic_doses_Tartrazine_mean' : numpy.mean, 'surfacic_doses_Tartrazine_std':numpy.std})
#    grCFarea= grCFm['area']

    grCFmeanSF = grCFm['area'].agg({'area_mean' : numpy.mean, 'area_std':numpy.std})
#    grCFmean = pandas.concat([grCFmeanSDE, grCFmeanSF])
    grCFmean = pandas.merge(grCFmeanSDE,grCFmeanSF)
#    grCFmean = pandas.merge(grCFmeanSDE,grCFarea)
#    grCFmean['SDESF']= grCFmean['surfacic_doses_Epoxiconazole_mean']*dataSFsum['area_sum']
#    grCFmean['SDESF']= grCFmean['surfacic_doses_Epoxiconazole_mean']*grCFm['area']
    grCFmean['SDESF']= grCFmean['surfacic_doses_Tartrazine_mean']*grCFmean['area_mean']
    print grCFmean['SDESF']
#    grCFmean = grCFm.agg({'surfacic_doses_Epoxiconazole_mean' : numpy.mean, 'surfacic_doses_Epoxiconazole_std' : numpy.std,'Foliar_Surface_mean' : numpy.mean,'Foliar_Surface_std':numpy.std})
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
    """
    data = data[data.ntop_cur <= 5]
    gr=data.groupby(['plant', 'ntop_cur', 'date', 'axe'], as_index=False)
    dmp=gr.sum() # dmp contient, pour chaque quadruplet ('plant', 'ntop_cur', 'date', 'axe'), 
                 # la somme des area, green_area, id, length, metamer, ntop, penetrated_doses_Epoxiconazole, 
                 # senesced_area, surfacic_doses_Epoxiconazole, n_max

    gr = dmp.groupby(['ntop_cur','plant','date'],as_index=False)
    # calcul, pour chaque triplet ('ntop_cur','plant','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest = gr['surfacic_doses_Tartrazine'].agg({'surfacic_doses_Tartrazine_mean' : numpy.mean, 'surfacic_doses_Tartrazine_std' : numpy.std})

    gr2 = dmp.groupby(['ntop_cur','date'],as_index=False)
    # calcul, pour chaque doublon ('ntop_cur','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest4 = gr2['surfacic_doses_Tartrazine'].agg({'surfacic_doses_Tartrazine_mean' : numpy.mean, 'surfacic_doses_Tartrazine_std' : numpy.std})

    gr3=dmp.groupby(['date'],as_index=False)
    ## calcul, pour chaque date et ntop_cur<4 de la surface totale verte et senecente et de l'interception totale
    grtest3 = gr3.sum()
    """   
    return grCFmean_global

#plot des donnees obs
def plot(name='Mercia'):

    if name is 'Mercia' or 'Rht3' :
#        data_file = 'ObsData2011.csv'
         data_file = 'grCFmean_global.csv'
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
    #    data = dat.groupby(['feuille']).moy_intercepte.sum()
        data = dat.groupby(['ntop_cur']).SDESF()
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

