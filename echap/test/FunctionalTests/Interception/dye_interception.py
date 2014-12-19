import os
import pandas
import numpy
import math
from itertools import cycle, islice

import alinea.echap.architectural_reconstructions as reconstructions
import alinea.echap.architectural_data as archidb
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

def repartition_at_application(name, appdate='T1', dose=1e4, nplants=30, density=1, dimension=1):
    """ 10000 l ha-1 is for 1 l/m2 """
    
    adel = Reconstructions.get_reconstruction(name=name, nplants=nplants, stand_density_factor = {name:density} , dimension=dimension)
    
    appdate_ref = ['T1-0.4', 'T1-0.2', 'T1', 'T1+0.2', 'T1+0.4',
            'T2-2.5', 'T2-2', 'T2-1.5', 'T2-1', 'T2-0.5', 'T2-0.4', 'T2-0.2', 'T2', 'T2+0.2', 'T2+0.4', 'T2+0.5', 'T2+1', 'T2+1.5', 'T2+2', 'T2+2.5']
    if appdate in appdate_ref:
        date, hs = HS_applications[name][appdate]
    else :
        hs = float(appdate); date = ''
    conv = HSconv[name] 
    age = conv.TT(hs)
    g = adel.setup_canopy(age)
    
    recorder = LeafElementRecorder()
    applications= 'date,dose, product_name\n%s 10:00:00, %f, Tartrazine'%(date, dose)
    application_data = pesticide_applications(applications)
    g,_ = pesticide_interception(g, interceptor, application_data)
    do_record(g, application_data, recorder, header={'TT':age, 'HS':hs})
    df =  recorder.get_records()
    print 'repartition_at_application df.columns before ', df.columns
    return adel, g, df

def aggregate_by_leaf(df):
    """
    Aggregate interceptioin data by leaf and add colmun 'deposits' (= area * surfacic doses)
    """
    def first_val(x):
        return x[x.first_valid_index()]
        
    df['deposit_Tartrazine'] = df['surfacic_doses_Tartrazine'] * df['area']
    df = df.convert_objects()
    gr = df.groupby(['HS', 'plant','axe','metamer'], as_index=False)
    return gr.agg({'TT': first_val, 
                   'area': numpy.sum,
                   'date': first_val,
                   'green_area': numpy.sum,
                   'senesced_area': numpy.sum,
                   'id':lambda(x): '_'.join(map(str,x)),
                   'length':numpy.sum,
                   'ntop':first_val,
                   'organ':first_val,
                   'penetrated_doses_Tartrazine': numpy.mean,
                   'surfacic_doses_Tartrazine':numpy.mean,
                   'deposit_Tartrazine': numpy.sum,
                   'lifetime': first_val,
                   'exposition': first_val,
                   })


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
    
   
def treatment(name='Tremie13', sim='T1', nplants=200, axis='MS', to_csv=False, density=1, dimension=1): 

    adel, g, df = repartition_at_application(name, sim, nplants=nplants, density=density, dimension=dimension)
    
    #Calcul age depuis emergence de la feuille    
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
    df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
    df_phenT['ntop_cur'] = df_phenT['nff'] - df_phenT['n']
    dates = [numpy.mean(df_phenT[df_phenT['ntop_cur']==lf]['tip']) for lf in range(1, int(max(df_phenT['n'])+2))]
    df_dates = pandas.DataFrame(dates, index = range(1, len(dates)+1), columns = ['leaf_emergence'])
    df_dates.index.name = 'ntop_cur'
    df_dates = df_dates.reset_index()
    
    if axis == 'MS':
        df = df[(df['axe'] == 'MS') & (df['hasEar'] == 1)] # do not consider aborting axes
    else:
        df = df[(df['hasEar'] == 1)]
    
    data = aggregate_by_leaf(df)
        
    # compute ntop_cur
    gr = data.groupby(['HS','plant','axe'], group_keys=False)
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
    data = gr.apply(_fun)
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1

    # aggregation by ntop cur
    sub = data.ix[:,('HS','TT','axe','metamer','ntop','ntop_cur','length','area','green_area','deposit_Tartrazine','exposition','lifetime')]
    dfmoy = sub.groupby(['HS','axe','ntop_cur']).mean().reset_index()
    dfsd = sub.groupby(['HS','axe','ntop_cur']).std().reset_index()
    
    #ajout colonne leaf_emergence dans dfmoy
    dfmoy = dfmoy.merge(df_dates, how='inner')

    if to_csv:
        data.to_csv(''.join(('data_',name,'_',sim,'_',axis,'.csv')))
        dfmoy.to_csv(''.join(('moy_',name,'_',sim,'_',axis,'.csv')))
        dfsd.to_csv(''.join(('std_',name,'_',sim,'_',axis,'.csv')))
             
    return adel.nplants, dfmoy, dfsd
    
#sensibilite aux nombres de plantes et de simulation
#ici : choix de faire 5 simulations a 30, 60, 100 et 200 plantes
def sensibilite_nplants(var_lst = ['Mercia','Rht3','Tremie12','Tremie13'], axis='MS', csv=True): #axis = MS or all 
    df_all = pandas.DataFrame()
    for var in var_lst: 
        df_sim = pandas.DataFrame()
        for stade in ['T1','T2']:
            if var=='Mercia':
                nplants_lst = [43,94,149,279] # soit n=31, 61, 93, 205
            elif var=='Rht3':
                nplants_lst = [49,100,150,285] # soit n=31, 58, 107, 201
            elif var=='Tremie12':
                nplants_lst = [35] # seulement n=25 pour Tremie12 car bug 
            elif var=='Tremie13':
                nplants_lst = [42,82,130,235] # soit n=28, 60, 100, 200  
            for n in nplants_lst : 
                x=1
                while x<=5: 
                    npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis=axis, to_csv=False)
                    print 'var = '+var+' stade = '+stade+' - nplants = '+str(npl)+' - simulation num = '+str(x)+'/5 ok' #verif etat avancement
                    dfmoy['var'] = var
                    dfmoy['dim'] = stade
                    dfmoy['nb_plantes_sim'] = npl
                    dfmoy['numero_sim'] = str(x)
                    df_sim = df_sim.append(dfmoy)
                    x+=1
        df_sim_gr = df_sim.groupby(['var', 'HS', 'dim', 'nb_plantes_sim', 'numero_sim', 'ntop_cur']).mean()
        df_sim_gr = df_sim_gr.reset_index()
        if csv == True:
            df_sim.to_csv('sensibilite_'+axis+'_all_'+var+'.csv')
            df_sim_gr.to_csv('sensibilite_'+axis+'_mean_'+var+'.csv')
        df_all = df_all.append(df_sim_gr)
        
    df_all = df_all.sort(['var', 'HS', 'ntop_cur'])
    df_all.to_csv('analyse_'+axis+'.csv', index=False)
    return df_all
 
#plot publication data+boite de petri 
def simulation_diff(var_lst = ['Mercia','Rht3','Tremie12','Tremie13'], axis='MS', csv=True): #axis = MS or all 
    df_sim = pandas.DataFrame()
    for var in var_lst:        
        for stade in ['T1','T2']:
            if var=='Tremie12':
                nplants_lst = [30] # !!!bug nplants
            else :
                nplants_lst = [200]
            for n in nplants_lst : 
                x=1
                while x<=5: 
                    npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis=axis, to_csv=False)
                    print 'var = '+var+' stade = '+stade+' - nplants = '+str(npl)+' - simulation num = '+str(x)+'/5 ok' #verif etat avancement
                    dfmoy['var'] = var
                    dfmoy['treatment'] = stade
                    dfmoy['nb_plantes_sim'] = npl
                    dfmoy['numero_sim'] = str(x)
                    df_sim = df_sim.append(dfmoy)
                    x+=1
    return df_sim
        
def plot_data(plot1=True, plot2=True):
    '''Plot1 = obs
    Plot2 = obs/area'''
    #obs + preparation de obs
    df_obs = idata.dye_interception() 
    df_obs.rename(columns={'name':'var'}, inplace=True) 
    df_obs.rename(columns={'N feuille':'ntop_cur'}, inplace=True)
    #sim
    df_sim = simulation_diff()
    df_sim = df_sim.groupby(['var', 'treatment', 'ntop_cur']).mean()
    df_sim = df_sim.reset_index()
    #merge
    df_all = df_obs.merge(df_sim)
    #column tartrazine/area
    df_all['mean/area'] = df_all['mean']/df_all['area']
    
    if plot1==True:
        bar_width = 0.2; opacity = 0.4
        val=[[0, 'T1'], [1, 'T2']]
        fig, axes = plt.subplots(nrows=1, ncols=2) 
    
        for x, date in val:        
            df_fin = df_all[df_all['treatment']==date] 
            
            df_obs1= df_all[df_all['feuille']=='F1'] 
            df_fin1 = df_obs1[df_obs1['treatment']==date] 
            
            df_obs2 = df_all[df_all['feuille']=='F2'] 
            df_fin2 = df_obs2[df_obs2['treatment']==date] 
            
            df_obs3 = df_all[df_all['feuille']=='F3'] 
            df_fin3 = df_obs3[df_obs3['treatment']==date] 
            
            df_obs4 = df_all[df_all['feuille']=='F4'] 
            df_fin4 = df_obs4[df_obs4['treatment']==date] 
            
            n_groups = len(df_fin1)
            index = numpy.arange(n_groups) 

            rects1 = axes[x].bar(index, df_fin1['mean'], bar_width, alpha=opacity, color='c', yerr=df_fin1['IC'], error_kw=dict(ecolor='c'))
            rects2 = axes[x].bar(index + bar_width, df_fin2['mean'], bar_width, alpha=opacity, color='m', yerr=df_fin2['IC'], error_kw=dict(ecolor='m'))
            rects3 = axes[x].bar(index + 2*bar_width, df_fin3['mean'], bar_width, alpha=opacity, color='g', yerr=df_fin3['IC'], error_kw=dict(ecolor='g'))
            rects4 = axes[x].bar(index + 3*bar_width, df_fin4['mean'], bar_width, alpha=opacity, color='r', yerr=df_fin4['IC'], error_kw=dict(ecolor='r'))   
            
            # Mise en forme
            axes[x].set_ylim(0, 6.5)
            axes[x].set_xlim(0, len(df_fin1))                 
            axes[x].set_xticks(index+bar_width)
            label = df_fin['var'].unique()
            axes[x].set_xticklabels( label.tolist(), rotation=90, fontsize='small' )
            if x == 0:
                axes[x].set_xlabel('varieties')
                axes[x].set_ylabel('g per g.cm-2')
            axes[x].text(0.2, 6.1, ''+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
            
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('F1', 'F2', 'F3', 'F4'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    if plot2==True:
        bar_width = 0.2; opacity = 0.4
        val=[[0, 'T1'], [1, 'T2']]
        fig, axes = plt.subplots(nrows=1, ncols=2) 
    
        for x, date in val:

            df_fin = df_all[df_all['treatment']==date] 
            
            df_obs1= df_all[df_all['feuille']=='F1'] 
            df_fin1 = df_obs1[df_obs1['treatment']==date] 
            
            df_obs2 = df_all[df_all['feuille']=='F2'] 
            df_fin2 = df_obs2[df_obs2['treatment']==date] 
            
            df_obs3 = df_all[df_all['feuille']=='F3'] 
            df_fin3 = df_obs3[df_obs3['treatment']==date] 
            
            df_obs4 = df_all[df_all['feuille']=='F4'] 
            df_fin4 = df_obs4[df_obs4['treatment']==date] 
            
            n_groups = len(df_fin1)
            index = numpy.arange(n_groups) 

            rects1 = axes[x].bar(index, df_fin1['mean/area'], bar_width, alpha=opacity, color='c')
            rects2 = axes[x].bar(index + bar_width, df_fin2['mean/area'], bar_width, alpha=opacity, color='m')
            rects3 = axes[x].bar(index + 2*bar_width, df_fin3['mean/area'], bar_width, alpha=opacity, color='g')
            rects4 = axes[x].bar(index + 3*bar_width, df_fin4['mean/area'], bar_width, alpha=opacity, color='r')   
            
            # Mise en forme
            axes[x].set_ylim(0, 1.2)
            axes[x].set_xlim(0, len(df_fin1))                 
            axes[x].set_xticks(index+bar_width)
            label = df_fin['var'].unique()
            axes[x].set_xticklabels( label.tolist(), rotation=90, fontsize='small' )
            if x == 0:
                axes[x].set_xlabel('varieties')
                axes[x].set_ylabel('g per g.cm-2 / area')
            axes[x].text(0.2, 1.1, ''+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
            
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('F1', 'F2', 'F3', 'F4'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
  
def plot_petri():
    #boite de petri, seulement Tremie12 et Tremie13
    petriT1_Tremie12, petriT2_Tremie12 = idata.Petri_data('Tremie12') 
    petriT1_Tremie13, petriT2_Tremie13 = idata.Petri_data('Tremie13')
    #traitement boite de Petri
    pdata_Tremie12 = archidb.Plot_data_Tremie_2011_2012()
    sowing_density_Tremie12 = pdata_Tremie12['sowing_density']
    lst1_Tremie12 = petriT1_Tremie12['rapportPoucentage(sol/emis)']*(10**4/sowing_density_Tremie12)
    lst2_Tremie12 = petriT2_Tremie12['rapportPoucentage(sol/emis)']*(10**4/sowing_density_Tremie12)
    mean1_Tremie12 = mean(lst1_Tremie12); mean2_Tremie12 = mean(lst2_Tremie12)
    IC1_Tremie12 = conf_int(lst1_Tremie12, perc_conf=95)
    IC2_Tremie12 =conf_int(lst2_Tremie12, perc_conf=95)
    
    pdata_Tremie13 = archidb.Plot_data_Tremie_2012_2013()
    sowing_density_Tremie13 = pdata_Tremie13['sowing_density']
    lst1_Tremie13 = petriT1_Tremie13['rapportPoucentage(sol/emis)']*(10**4/sowing_density_Tremie13)
    lst2_Tremie13 = petriT2_Tremie13['rapportPoucentage(sol/emis)']*(10**4/sowing_density_Tremie13)
    mean1_Tremie13 = mean(lst1_Tremie13); mean2_Tremie13 = mean(lst2_Tremie13)
    IC1_Tremie13 = conf_int(lst1_Tremie13, perc_conf=95)
    IC2_Tremie13 =conf_int(lst2_Tremie13, perc_conf=95)
    
    bar_width = 0.45; opacity = 0.4
    val=[[0, 'T1'], [1, 'T2']]

    fig, axes = plt.subplots(nrows=1, ncols=2) 
    for x, date in val:

        if x==0:
            df_petri = pandas.DataFrame({'var':['Mercia','Rht3','Tremie12','Tremie13'], 'mean':[None, None, mean1_Tremie12, mean1_Tremie13], 'IC':[None, None, IC1_Tremie12, IC1_Tremie13]})
        else:
            df_petri = pandas.DataFrame({'var':['Mercia','Rht3','Tremie12','Tremie13'], 'mean':[None, None, mean2_Tremie12, mean2_Tremie13], 'IC':[None, None, IC2_Tremie12, IC2_Tremie13]})
        
        n_groups = len(df_petri)
        index = numpy.arange(n_groups) 

        rects1 = axes[x].bar(index + bar_width, df_petri['mean'], bar_width, alpha=opacity, color='b', yerr=df_petri['IC'], error_kw=dict(ecolor='b'))    
        
        # Mise en forme
        axes[x].set_ylim(0, 15)
        axes[x].set_xlim(0, len(df_petri))                 
        axes[x].set_xticks(index+bar_width)
        label = df_petri['var'].unique()
        axes[x].set_xticklabels( label.tolist(), rotation=90, fontsize='small' )
        if x == 0:
            axes[x].set_xlabel('varieties')
            axes[x].set_ylabel('g per g.cm-2')
        axes[x].text(0.2, 14, ''+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
        
        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    axes[x].legend((rects1[0], ''), ('Petri',''), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
 
# simulation interception des talles et des MB (4 valeurs des simulations cote a cote pour les maquettes) 
def plot_diff(varieties=['Mercia']):
    for var in varieties: 
        file = simulation_diff(var_lst=[var], axis='all')
        file.to_csv('plot_diff_'+var+'.csv', index=False)
        file = file[file['area']>0]
        file.to_csv('plot_diff_withoutArea0_'+var+'.csv', index=False)
        df_MS = file[file['axe']=='MS']
        df_MS = df_MS.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        #df_MS_count = df_MS.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        #df_MS = df_MS.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_MS = df_MS.reset_index(); df_MS = df_MS[df_MS['ntop_cur']<=5]
        df_all = file
        df_all = df_all.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        #df_all_count = df_all.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        #df_all = df_all.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_all = df_all.reset_index(); df_all = df_all[df_all['ntop_cur']<=5]
        df_talles = file[file['axe']!='MS']
        df_talles = df_talles.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        #df_talles_count = df_talles.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        #df_talles = df_talles.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_talles = df_talles.reset_index(); df_talles = df_talles[df_talles['ntop_cur']<=5]
        df_T1 = file[file['axe']=='T1']
        df_T1 = df_T1.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        #df_T1_count = df_T1.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        #df_T1 = df_T1.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_T1 = df_T1.reset_index(); df_T1 = df_T1[df_T1['ntop_cur']<=5]
        
        bar_width = 0.2; opacity = 0.4
        if var=='Mercia':
            val = [[0,9.74],[1,12.8]]
        elif var == 'Rht3' :
            val = [[0,9.15],[1,12.48]]
        elif var == 'Tremie12' :
            val = [[0,10.98],[1,12.63]]
        else:
            val = [[0,8.7],[1,11.04]]
        fig, axes = plt.subplots(nrows=1, ncols=2) 
        for x, date in val:
            df_MS_date = df_MS[df_MS['HS']==date] 
            df_all_date = df_all[df_all['HS']==date] 
            df_talles_date = df_talles[df_talles['HS']==date] 
            df_T1_date = df_T1[df_T1['HS']==date] 
            
            n_groups = len(df_MS_date)
            index = numpy.arange(n_groups) 

            rects1 = axes[x].bar(index, df_MS_date['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c'])
            rects2 = axes[x].bar(index + bar_width, df_all_date['deposit_Tartrazine'], bar_width, alpha=opacity, color=['m'])
            rects3 = axes[x].bar(index + 2*bar_width, df_talles_date['deposit_Tartrazine'], bar_width, alpha=opacity, color=['g'])
            rects4 = axes[x].bar(index + 3*bar_width, df_T1_date['deposit_Tartrazine'], bar_width, alpha=opacity, color=['y'])
            '''rects1 = axes[x].bar(index, df_MS_date['area'], bar_width, alpha=opacity, color=['c'])
            rects2 = axes[x].bar(index + bar_width, df_all_date['area'], bar_width, alpha=opacity, color=['m'])
            rects3 = axes[x].bar(index + 2*bar_width, df_talles_date['area'], bar_width, alpha=opacity, color=['g'])
            rects4 = axes[x].bar(index + 3*bar_width, df_T1_date['area'], bar_width, alpha=opacity, color=['y'])'''
            
            #Mise en forme
            axes[x].set_ylim(0, 6)
            axes[x].set_xlim(0, len(df_MS_date))                 
            axes[x].set_xticks(index+bar_width)
            df_all_date['label'] = df_all_date['ntop_cur'].astype(str)
            axes[x].set_xticklabels( df_all_date['label'].tolist(), rotation=90, fontsize='small' )
            if x == 0:
                axes[x].set_xlabel('ntop')
            axes[x].text(0.4, 5.5, var+' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

            fig.suptitle('Simulation somme surfaces des talles et des MB', fontsize=10)
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Interception MB', 'Interception MB + talles', 'Interception talles', 'Interception T1'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})

#effet dimension
def simulation_dimension(name='Tremie12', n_sim=5, n_plt=200, csv=True): 
    if name=='Tremie12':
        n = 30
    else :
        n = n_plt
        
    lst = [[name,'T1'], [name,'T2']]
        
    df_sim = pandas.DataFrame()
    x=0
    while x<=n_sim:
        for var, stade in lst :      
            dim_lst = [[0.5,'/2'], [0.67,'/1.5'], [1,'1'], [1.5,'x1.5'], [2,'x2']]
            for dim, dim_label in dim_lst:
                npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis='MS', dimension=dim, to_csv=False)
                print 'var = '+var+' stade = '+dim_label+' - nbre de plantes sim = '+str(npl)
                    
                dfmoy['var'] = var; dfmoy['dim'] = dim_label
                df_sim = df_sim.append(dfmoy)
        x+=1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'dim', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    if csv == True:
        df_sim.to_csv('dim_all_'+name+'.csv')
        df_sim_gr.to_csv('dim_synth_'+name+'.csv')
    return df_sim_gr
    
def plot_dimension(plot1=True, plot2=False, plot3=True, varieties=['Tremie13'], n_sim=1, n_plt=30):
    '''
    Plot 1 =  Barplot observation / simulation avec ntop_cur en abscisse et (dfmoy['deposit_Tartrazin']  versus obs 'moyenne' par genotype et par date interception
    Plot 2 = Barplot de controle 'area' versus scan data (avec la colone area de dfmoy)
    Plot 3 = Barplot observation/simulation de deposit/area et obs/area 
    
    !!! On plotte seulement les 3 feuilles du haut !!!
    '''
    
    #base
    df_scan_tremie12 = archidb.treatment_scan('Tremie12'); df_scan_tremie12['var']='Tremie12'
    df_scan_tremie13 = archidb.treatment_scan('Tremie13'); df_scan_tremie13['var']='Tremie13'
    df_scan = df_scan_tremie12.append(df_scan_tremie13)
    df_obs = idata.dye_interception()
    df_sim = pandas.DataFrame()
    for n in varieties:
        df_sim_var = simulation_dimension(name=n, n_sim=n_sim, n_plt=n_plt)
        df_sim = df_sim.append(df_sim_var)
    #df_sim = 'dim_tremie13.csv'
    #df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')
    
    # PLOT 1
    if plot1 is True :
        for var in varieties:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=3]
            #plot
            bar_width = 0.4; opacity = 0.4
            if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]
            elif var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date] 
                #index_dim pour classer les dim dans l ordre souhaite
                df_fin['index_dim'] = df_fin['dim']
                df_fin.ix[df_fin.dim.isin(['/2']), 'index_dim'] = 1
                df_fin.ix[df_fin.dim.isin(['/1.5']), 'index_dim'] = 2
                df_fin.ix[df_fin.dim.isin(['1']), 'index_dim'] = 3
                df_fin.ix[df_fin.dim.isin(['x1.5']), 'index_dim'] = 4
                df_fin.ix[df_fin.dim.isin(['x2']), 'index_dim'] = 5
                #rects2 pour ne tracer les obs que quand dim = normale
                df_fin['mean_dim'] = df_fin['mean']
                df_fin.ix[~df_fin.dim.isin(['1']), 'mean_dim'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur','HS','index_dim'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups) 

                rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g'])
                rects2 = axes[x].bar(index + bar_width, df_fin['mean_dim'], bar_width, alpha=opacity, color='y')
                # Mise en forme
                axes[x].set_ylim(0, 17)
                axes[x].set_xlim(0, len(df_fin))                 
                axes[x].set_xticks(index+bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['label'] = df_fin['dim'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                if x == 0:
                    axes[x].set_xlabel('dim/ntop')
                axes[x].text(0.4, 16.2, var+' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

                fig.suptitle('Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur', fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    # PLOT 2
    if plot2 is True :
        for var in varieties:
            if var=='Tremie12' or var=='Tremie13':
                df_scan_var = df_scan[df_scan['var']==var]
                df_sim_var = df_sim[df_sim['var']==var]
                df_scan_var.HS=map(str,df_scan_var.HS)
                df_sim_var.HS=map(str,df_sim_var.HS)
                df_all = df_sim_var.merge(df_scan_var.ix[:,['HS','ntop_cur','Area A_bl']], how='outer')
                df_all = df_all[df_all['ntop_cur']<=3]
                #plot
                bar_width = 0.4; opacity = 0.4
                fig, axes = plt.subplots(nrows=1, ncols=2) 
                if var=='Mercia':
                    val = [[0,9.74],[1,12.8]]
                elif var == 'Rht3' :
                    val = [[0,9.15],[1,12.48]]
                elif var == 'Tremie12' :
                    val = [[0,10.98],[1,12.63]]
                else:
                    val = [[0,8.7],[1,11.04]]
                for x, date in val:
                    df_all.HS=map(float,df_all.HS)
                    df_fin = df_all[df_all['HS']==date]
                    #index_dim
                    df_fin['index_dim'] = df_fin['dim']
                    df_fin.ix[df_fin.dim.isin(['/2']), 'index_dim'] = 1
                    df_fin.ix[df_fin.dim.isin(['/1.5']), 'index_dim'] = 2
                    df_fin.ix[df_fin.dim.isin(['1']), 'index_dim'] = 3
                    df_fin.ix[df_fin.dim.isin(['x1.5']), 'index_dim'] = 4
                    df_fin.ix[df_fin.dim.isin(['x2']), 'index_dim'] = 5
                    #rects2 
                    df_fin['Area A_bl_dim'] = df_fin['Area A_bl']
                    df_fin.ix[~df_fin.dim.isin(['1']), 'Area A_bl_dim'] = 0
                    #---
                    df_fin = df_fin.sort(['ntop_cur','HS','index_dim'])
                    n_groups = len(df_fin)
                    index = numpy.arange(n_groups) 
                    rects1 = axes[x].bar(index, df_fin['Area A_bl_dim'], bar_width, alpha=opacity, color='y')
                    rects2 = axes[x].bar(index + bar_width, df_fin['area'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g'])
                    
                    # Mise en forme
                    axes[x].set_ylim(0, 50)
                    axes[x].set_xlim(0, len(df_fin))          
                    axes[x].set_xticks(index+bar_width)
                    df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                    df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                    df_fin['label'] = df_fin['dim'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                    axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                    if x==0:
                        axes[x].set_xlabel('dim/ntop'); axes[x].set_ylabel('area (cm2)')
                    axes[x].text(0.4, 48, var+' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                    
                    fig.suptitle('Surface scan/sim par ntop_cur', fontsize=10)
                    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                axes[x].legend((rects1[0], rects2[0], rects2[5], rects2[10]), ('Scan', 'Sim F1', 'Sim F2', 'Sim F3'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )
            
    # PLOT 3
    if plot3 is True :
        for var in varieties:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=3]
            df_all['tartrazine/area'] = df_all['deposit_Tartrazine']/df_all['area']
            df_all['mean/area'] = df_all['mean']/df_all['area']
            #plot
            bar_width = 0.4; opacity = 0.4
            if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]
            elif var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date]
                #index_dim
                df_fin['index_dim'] = df_fin['dim']
                df_fin.ix[df_fin.dim.isin(['/2']), 'index_dim'] = 1
                df_fin.ix[df_fin.dim.isin(['/1.5']), 'index_dim'] = 2
                df_fin.ix[df_fin.dim.isin(['1']), 'index_dim'] = 3
                df_fin.ix[df_fin.dim.isin(['x1.5']), 'index_dim'] = 4
                df_fin.ix[df_fin.dim.isin(['x2']), 'index_dim'] = 5
                #rects2
                df_fin['mean_dim'] = df_fin['mean/area']
                df_fin.ix[~df_fin.dim.isin(['1']), 'mean_dim'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur','HS','index_dim'])      
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)   
                
                rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g'])
                rects2 = axes[x].bar(index + bar_width, df_fin['mean_dim'], bar_width, alpha=opacity, color='y')   
                # Mise en forme
                axes[x].set_ylim(0, 0.35)
                axes[x].set_xlim(0, len(df_fin)) 
                axes[x].set_xticks(index+bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['dim'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                if x == 0:
                    axes[x].set_xlabel('dim/ntop')
                axes[x].text(0.4, 0.33, var+' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                    
                fig.suptitle('Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur', fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )

    return df_scan, df_obs, df_sim
 
#effet densite 
def simulation_density(name='Tremie12', n_sim=5, n_plt=200, csv=True):
    if name=='Tremie12':
        n = 30
    else :
        n = n_plt
        
    lst = [[name,'T1'], [name,'T2']]
        
    df_sim = pandas.DataFrame()
    x=1
    while x<=n_sim:
        for var, stade in lst :  
            dens_lst = [0.1,1,2.5,5,8,9,10] #tjrs mettre 7 valeurs sinon decalage dans les couleurs du plot, penser a modifier argument color dans rects de plot_density
            for dens in dens_lst:         
                npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis='MS', to_csv=False, density=dens)
                print 'var = '+var+' stade = '+stade+' density = '+str(dens)+' - nbre de plantes sim = '+str(npl)
                dfmoy['var'] = var
                dfmoy['density'] = dens
                df_sim = df_sim.append(dfmoy)
        x+=1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'density', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    if csv==True:
        df_sim.to_csv('density_all_'+name+'.csv')
        df_sim_gr.to_csv('density_synth_'+name+'.csv')
    return df_sim_gr
  
def plot_density(plot1=True, plot2=False, plot3=True, plot4=False, varieties=['Tremie12','Tremie13'], n_sim=1, n_plt=30):
    '''
    Plot 1 =  Barplot observation / simulation avec ntop_cur en abscisse et (dfmoy['deposit_Tartrazin']  versus obs 'moyenne' par genotype et par date interception
    Plot 2 = Barplot de controle 'area' versus scan data (avec la colone area de dfmoy)
    !!! plot 2 = PAS DE SCAN POUR MERCIA ET RHT3 donc pas de graph genere (test si plot2=True, ne posera pas de pb)
    Plot 3 = Barplot observation/simulation de deposit/area et obs/area
    Plot 4 = 2 graphes supplémentaires récapitulatifs avec en x = âge des feuilles (en dd depuis la ligulation, quelque soit le numéro des feuilles) ET Ntop (2 graphes donc) et en y = interception feuille par cm2 
    '''
    
    #base
    df_scan_tremie12 = archidb.treatment_scan('Tremie12'); df_scan_tremie12['var']='Tremie12'
    df_scan_tremie13 = archidb.treatment_scan('Tremie13'); df_scan_tremie13['var']='Tremie13'
    df_scan = df_scan_tremie12.append(df_scan_tremie13)
    df_obs = idata.dye_interception()
    df_sim = pandas.DataFrame()
    for n in varieties:
        df_sim_var = simulation_density(name=n, n_sim=n_sim, n_plt=n_plt)
        df_sim = df_sim.append(df_sim_var)
    # test
    #df_sim = 'density_synth_tremie12.csv'
    #df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')
    
    # PLOT 1
    if plot1 is True :
        for var in varieties:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=5]
            #plot
            bar_width = 0.4; opacity = 0.4
            if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]
            elif var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date] 
                #mean_dens pour dessiner les obs seulement qd densite=1
                df_fin['mean_dens'] = df_fin['mean']
                df_fin.ix[~df_fin.density.isin([1]), 'mean_dens'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur','HS','density'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups) 

                rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','m','m','m','m','m','m','m','g','g','g','g','g','g','g','r','r','r','r','r','r','r','b','b','b','b','b','b','b'])
                rects2 = axes[x].bar(index + bar_width, df_fin['mean_dens'], bar_width, alpha=opacity, color='y')
                # Mise en forme
                axes[x].set_ylim(0, 8)
                axes[x].set_xlim(0, len(df_fin))                 
                axes[x].set_xticks(index+bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['density'].astype(str) + ' - ' + df_fin['TT'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                if x == 0:
                    axes[x].set_xlabel('densite/TT/ntop')
                axes[x].text(0.4, 7.4, var + ' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

                fig.suptitle('Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[7], rects1[14], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    # PLOT 2 - !!! PAS DE SCAN POUR MERCIA ET RHT3 donc pas de graph genere
    if plot2 is True :
        for var in varieties:
            if var=='Tremie12' or var=='Tremie13':
                df_scan_var = df_scan[df_scan['var']==var]
                df_sim_var = df_sim[df_sim['var']==var]
                df_scan_var.HS=map(str,df_scan_var.HS)
                df_sim_var.HS=map(str,df_sim_var.HS)
                df_all = df_sim_var.merge(df_scan_var.ix[:,['HS','ntop_cur','Area A_bl']], how='outer')
                df_all = df_all[df_all['ntop_cur']<=5]
                #plot
                bar_width = 0.4; opacity = 0.4
                fig, axes = plt.subplots(nrows=1, ncols=2) 
                if var == 'Tremie12' :
                    val = [[0,10.98],[1,12.63]]
                else:
                    val = [[0,8.7],[1,11.04]]
                for x, date in val:
                    df_all.HS=map(float,df_all.HS)
                    df_fin = df_all[df_all['HS']==date]
                    #Area A_bl dens pour dessiner les obs seulement qd densite=1
                    df_fin['Area A_bl dens'] = df_fin['Area A_bl']
                    df_fin.ix[~df_fin.density.isin([1]), 'Area A_bl dens'] = 0
                    df_fin = df_fin.sort(['ntop_cur', 'HS','density'])
                    n_groups = len(df_fin)
                    index = numpy.arange(n_groups) 

                    rects1 = axes[x].bar(index, df_fin['Area A_bl dens'], bar_width, alpha=opacity, color='y')
                    rects2 = axes[x].bar(index + bar_width, df_fin['area'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','m','m','m','m','m','m','m','g','g','g','g','g','g','g','r','r','r','r','r','r','r','b','b','b','b','b','b','b'])
                    
                    # Mise en forme
                    axes[x].set_ylim(0, 30)
                    axes[x].set_xlim(0, len(df_fin))          
                    axes[x].set_xticks(index+bar_width)
                    df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                    df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                    df_fin['label'] = df_fin['density'].astype(str) + ' - ' + df_fin['TT'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                    axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='x-small' )
                    if x==0:
                        axes[x].set_xlabel('densite/TT/ntop'); axes[x].set_ylabel('area (cm2)')
                    axes[x].text(0.4, 28.5, var + ' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                    fig.suptitle('Surface scan/sim par ntop_cur pour '+var, fontsize=10)
                    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                axes[x].legend((rects1[0], rects2[0], rects2[7], rects2[14]), ('Scan', 'Sim F1', 'Sim F2', 'Sim F3'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )
            else :
                print 'Pas de scan pour variete : '+var
            
    # PLOT 3
    if plot3 is True :
        for var in varieties:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=5]
            df_all['tartrazine/area'] = df_all['deposit_Tartrazine']/df_all['area']
            df_all['mean/area'] = df_all['mean']/df_all['area']
            #plot
            bar_width = 0.4; opacity = 0.4
            if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]
            elif var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date]
                #mean_dens pour dessiner les obs seulement qd densite=1
                df_fin['mean/area_dens'] = df_fin['mean/area']
                df_fin.ix[~df_fin.density.isin([1]), 'mean/area_dens'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur','HS','density'])      
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)   
                
                rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','m','m','m','m','m','m','m','g','g','g','g','g','g','g','r','r','r','r','r','r','r','b','b','b','b','b','b','b'])
                rects2 = axes[x].bar(index + bar_width, df_fin['mean/area_dens'], bar_width, alpha=opacity, color='y')   
                # Mise en forme
                axes[x].set_ylim(0, 0.5)
                axes[x].set_xlim(0, len(df_fin)) 
                axes[x].set_xticks(index+bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['density'].astype(str) + ' - ' + df_fin['TT'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='x-small' )
                if x == 0:
                    axes[x].set_xlabel('densite/TT/ntop')
                axes[x].text(0.4, 0.45, var + ' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                fig.suptitle('Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[7], rects1[14], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )

    # PLOT 4
    if plot4 is True :
        for var in varieties:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=5]
            #plot
            bar_width = 0.4; opacity = 0.4
            if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]
            elif var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date] 
                #mean_dens pour dessiner les obs seulement qd densite=1
                df_fin['mean_dens'] = df_fin['mean']
                df_fin.ix[~df_fin.density.isin([1]), 'mean_dens'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur', 'leaf_emergence'], ascending=[1,0])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
                rects2 = axes[x].bar(index + bar_width, df_fin['mean_dens'], bar_width, alpha=opacity, color='y')
                # Mise en forme
                axes[x].set_ylim(0, 8)
                axes[x].set_xlim(0, len(df_fin))                 
                axes[x].set_xticks(index+bar_width)
                df_fin['leaf_emergence'] = numpy.round(df_fin['leaf_emergence'], decimals=1)
                df_fin['label'] = df_fin['leaf_emergence'].astype(str) + ' - ' + df_fin['density'].astype(str) 
                axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                if x == 0:
                    axes[x].set_xlabel('leaf_since_emergence/densite')
                axes[x].text(0.4, 7.4, var + ' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

                fig.suptitle('Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    return df_scan, df_obs, df_sim  
 
#effet HS 
def simulation_HS(name='Tremie12', n_sim=5, n_plt=200, csv=True):
    if name=='Tremie12':
        n = 30
    else :
        n = n_plt
        
    lst = [[name,'T1-0.4'], [name,'T1-0.2'], [name,'T1'], [name,'T1+0.2'], [name,'T1+0.4'],
            [name,'T2-2.5'], [name,'T2-2'], [name,'T2-1.5'], [name,'T2-1'], [name,'T2-0.5'], [name,'T2-0.4'], [name,'T2-0.2'], [name,'T2'], [name,'T2+0.2'], [name,'T2+0.4'], [name,'T2+0.5'], [name,'T2+1'], [name,'T2+1.5'], [name,'T2+2'], [name,'T2+2.5']]


    df_sim = pandas.DataFrame()
    x=1
    while x<=n_sim:
        for var, stade in lst :          
            npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis='MS', to_csv=False)
            print 'var = '+var+' stade = '+stade+' - nbre de plantes sim = '+str(npl)
            dfmoy['var'] = var
            df_sim = df_sim.append(dfmoy)
        x+=1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    if csv==True:
        df_sim.to_csv('HS_all_'+name+'.csv')
        df_sim_gr.to_csv('HS_synth_'+name+'.csv')
    return df_sim_gr
    
def plot_HS(plot1=True, plot2=False, plot3=True, varieties=['Mercia', 'Rht3', 'Tremie12', 'Tremie13'], n_sim=1, n_plt=30):
    '''
    Plot 1 =  Barplot observation / simulation avec ntop_cur en abscisse et (dfmoy['deposit_Tartrazin']  versus obs 'moyenne' par genotype et par date interception
    Plot 2 = Barplot de controle 'area' versus scan data (avec la colone area de dfmoy)
    Plot 3 = Barplot observation/simulation de deposit/area et obs/area 
    '''
    
    #base
    df_scan_tremie12 = archidb.treatment_scan('Tremie12'); df_scan_tremie12['var']='Tremie12'
    df_scan_tremie13 = archidb.treatment_scan('Tremie13'); df_scan_tremie13['var']='Tremie13'
    df_scan = df_scan_tremie12.append(df_scan_tremie13)
    df_obs = idata.dye_interception()
    df_sim = pandas.DataFrame()
    for n in varieties:
        df_sim_var = simulation_HS(name=n, n_sim=n_sim, n_plt=n_plt)
        df_sim = df_sim.append(df_sim_var)
    # test 
    #df_sim = 'HS_synth_tremie13.csv'
    #df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')
    
    # PLOT 1
    if plot1 is True :
        for var in varieties:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=5]
            #plot
            bar_width = 0.4; opacity = 0.4
            if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]
            elif var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date] 
                if x == 0:
                    for f in [-0.4,-0.2,0.2,0.4]:
                        df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                else :
                    for f in [-2.5,-2,-1.5,-1,-0.5,-0.4,-0.2,0.2,0.4,0.5,1,1.5,2,2.5]:
                        df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                df_fin = df_fin.sort(['ntop_cur','HS'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups) 

                if date == 8.7 or date == 10.98 or date==9.15 or date==9.74 :
                    rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
                else :
                    rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','g','g','g','g','g','g','g','g','g','g','g','g','g','g','g','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b'])
   
                rects2 = axes[x].bar(index + bar_width, df_fin['mean'], bar_width, alpha=opacity, color='y')
                # Mise en forme
                axes[x].set_ylim(0, 7)
                axes[x].set_xlim(0, len(df_fin))                 
                axes[x].set_xticks(index+bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['HS'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                if x == 0:
                    axes[x].set_xlabel('HS/ntop')
                axes[x].text(0.4, 6.5, var + ' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                fig.suptitle('Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[15], rects1[30], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    # PLOT 2
    if plot2 is True :
        for var in varieties:
            if var=='Tremie12' or var=='Tremie13':
                df_scan_var = df_scan[df_scan['var']==var]
                df_sim_var = df_sim[df_sim['var']==var]
                df_scan_var.HS=map(str,df_scan_var.HS)
                df_sim_var.HS=map(str,df_sim_var.HS)
                df_all = df_sim_var.merge(df_scan_var.ix[:,['HS','ntop_cur','Area A_bl']], how='outer')
                df_all = df_all[df_all['ntop_cur']<=5]
                #plot
                bar_width = 0.4; opacity = 0.4
                fig, axes = plt.subplots(nrows=1, ncols=2) #fig, axes = plt.subplots(nrows=1, ncols=len(df_all['HS'].unique()))
                if var == 'Tremie12':
                    val = [[0,10.98],[1,12.63]]
                elif var == 'Tremie13' :
                    val = [[0,8.7],[1,11.04]]
                for x, date in val:
                    df_all.HS=map(float,df_all.HS)
                    df_fin = df_all[df_all['HS']==date] 
                    if x == 0:
                        for f in [-0.4,-0.2,0.2,0.4]:
                            df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                    else :
                        for f in [-0.4,-0.2,0.2,0.4,0.5,1,1.5,2,2.5]:
                            df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                    df_fin = df_fin.sort(['ntop_cur', 'HS'])
                    n_groups = len(df_fin)
                    index = numpy.arange(n_groups)   
                    rects1 = axes[x].bar(index, df_fin['Area A_bl'], bar_width, alpha=opacity, color='y')
                    
                    if date == 10.98 or date == 8.7  :
                        rects2 = axes[x].bar(index + bar_width, df_fin['area'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
                    else :
                        rects2 = axes[x].bar(index + bar_width, df_fin['area'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','g','g','g','g','g','g','g','g','g','g','g','g','g','g','g','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b'])
                    
                    # Mise en forme
                    axes[x].set_ylim(0, 30)
                    axes[x].set_xlim(0, len(df_fin))          
                    axes[x].set_xticks(index+bar_width)
                    df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                    df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                    df_fin['label'] = df_fin['HS'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                    axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                    if x==0:
                        axes[x].set_xlabel('HS/ntop'); axes[x].set_ylabel('area (cm2)')
                    axes[x].text(0.4, 28.5, var + ' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                    fig.suptitle('Surface scan/sim par ntop_cur pour '+var, fontsize=10)
                    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                axes[x].legend((rects2[0], rects2[15], rects2[30], rects1[0]), ( 'Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )
            else :
                print 'Pas de scan pour variete : '+var
            
    # PLOT 3
    if plot3 is True :
        for var in varieties:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean','sd']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=5]
            df_all['tartrazine/area'] = df_all['deposit_Tartrazine']/df_all['area']
            df_all['mean/area'] = df_all['mean']/df_all['area']
            #plot
            bar_width = 0.4; opacity = 0.4
            if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]
            elif var == 'Tremie12':
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date] 
                if x == 0:
                    for f in [-0.4,-0.2,0.2,0.4]:
                        df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                else :
                    for f in [-2.5,-2,-1.5,-1,-0.5,-0.4,-0.2,0.2,0.4,0.5,1,1.5,2,2.5]:
                        df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                df_fin = df_fin.sort(['ntop_cur','HS'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)   
                
                #rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','m','g','r','b'])
                if date == 10.98 or date == 8.7 or date==9.15 or date==9.74:
                    rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
                else :
                    rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','g','g','g','g','g','g','g','g','g','g','g','g','g','g','g','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b'])
                    
                rects2 = axes[x].bar(index + bar_width, df_fin['mean/area'], bar_width, alpha=opacity, color='y')   
                # Mise en forme
                axes[x].set_ylim(0, 0.35)
                axes[x].set_xlim(0, len(df_fin)) 
                axes[x].set_xticks(index+bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['HS'].astype(str) + ' - ' + df_fin['metamer'].astype(str)
                axes[x].set_xticklabels( df_fin['label'].tolist(), rotation=90, fontsize='small' )
                if x == 0:
                    axes[x].set_xlabel('HS/ntop')
                axes[x].text(0.4, 0.32, var + ' - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                fig.suptitle('Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[15], rects1[30], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )
            
    return df_scan, df_obs, df_sim
    
#comparer ntop_cur=1 pour des ages de feuilles identiques
def sim_HScom(var='Mercia', stade='T1', axis='MS'):
    df_sim = pandas.DataFrame()
    x=1
    while x<=1: 
        if var=='Tremie12':
            n=30
        else:
            n=30
        npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis=axis, to_csv=False)
        print 'var = '+var+' stade = '+stade+' - nbre de plantes sim = '+str(npl)+ ' - nbr simul = '+str(x)+'/5'
        dfmoy['var'] = var
        if stade=='T2':
            dfmoy['stade'] = 'T2_'+var
        else :
            dfmoy['stade'] = stade
        df_sim = df_sim.append(dfmoy)
        x+=1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'stade', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr
    
def plot_HScom(varieties=['Mercia','Rht3','Tremie12','Tremie13'], appdate_lst = ['T1','T2'], delta_opt=True, axis='MS', nplants=30, plot1=True, plot2=True): #delta = FALSE (6) ou TRUE
    '''
    Comparer les F1 T1 et F1 T2 pour les 4 varietes au HS de Mercia, au HS de Rht3, au HS de Tremie12 et au HS de Tremie13
    Methode : (ex. pour Mercia, T1)
        delta_HS_Mercia = HS T1 Mercia - nff moyen de Mercia
    Pour les autres varietes, faire :
        HS = nff moyen + delta_HS_Mercia
        
    Plot 1 = Tartrazine
    Plot 2 = Tartrazine / area
    '''
    df_all = pandas.DataFrame()
    #Ref Mercia
    for name in varieties :
        for appdate in appdate_lst:
            # HS
            date, hs = HS_applications[name][appdate]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=name, nplants=nplants)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            # calcul du delta
            #delta = hs - mean['nff']
            if delta_opt==True:
                delta = hs - mean['nff']
            else:
                delta = 6
            #verif
            print 'nff moyen ref '+name+' = '+str(mean['nff'])
            print 'delta ref '+name+' = '+str(delta)
            #simulation
            #hs_ref = hs+delta
            df_ref = sim_HScom(var=name, stade=str(hs), axis=axis)
            df_ref['dim'] = appdate+'_'+name
            df_ref['delta'] = delta
            df_all = df_all.append(df_ref)
            
            if name == 'Mercia':
                name_other = ['Rht3', 'Tremie12', 'Tremie13']
            elif name == 'Rht3':
                name_other = ['Mercia', 'Tremie12', 'Tremie13']
            elif name == 'Tremie12':
                name_other = ['Mercia', 'Rht3', 'Tremie13']
            elif name == 'Tremie13':
                name_other = ['Mercia', 'Rht3', 'Tremie12']
            for nameo in name_other:
                # nff moyen
                adel = Reconstructions.get_reconstruction(name=nameo, nplants=nplants)
                df_phenT2 = adel.phenT()
                nffs = df_phenT2.set_index('plant').groupby(level=0).count()['n']
                df_phenT2['nff'] = [nffs[v] for v in df_phenT2['plant']]
                mean_other = df_phenT2.mean()
                #HS
                HS_other = mean_other['nff'] + delta
                #verif
                print 'nff moyen '+nameo+' = '+str(mean_other['nff'])
                print 'HS other '+nameo+' = '+str(HS_other)
                #simulation
                df_other = sim_HScom(var=nameo, stade=str(HS_other), axis=axis)
                df_other['dim'] = appdate+'_'+name
                df_all = df_all.append(df_other)
    
    if plot1==True:
        bar_width = 0.2; opacity = 0.4
        fig, axes = plt.subplots(nrows=1, ncols=2) 

        for x, date in [[0,'T1'], [1,'T2']]:
        
            df_sim = df_all[df_all['ntop_cur']==1]
            
            df_sim_HSMercia = df_sim[df_sim['dim'] == date+'_Mercia']; df_sim_HSMercia = df_sim_HSMercia.reset_index()  

            df_sim_HSRht3 = df_sim[df_sim['dim'] == date+'_Rht3']; df_sim_HSRht3 = df_sim_HSRht3.sort(['var']); df_sim_HSRht3 = df_sim_HSRht3.reset_index()
            
            df_sim_HSTremie12 = df_sim[df_sim['dim'] == date+'_Tremie12']; df_sim_HSTremie12 = df_sim_HSTremie12.sort(['var']); df_sim_HSTremie12 = df_sim_HSTremie12.reset_index()
            
            df_sim_HSTremie13 = df_sim[df_sim['dim'] == date+'_Tremie13']; df_sim_HSTremie13 = df_sim_HSTremie13.sort(['var']); df_sim_HSTremie13 = df_sim_HSTremie13.reset_index()
                        
            n_groups = len(df_sim_HSMercia)
            index = numpy.arange(n_groups) 

            rects1 = axes[x].bar(index, df_sim_HSMercia['deposit_Tartrazine'], bar_width, alpha=opacity, color=['r'])
            rects2 = axes[x].bar(index + bar_width, df_sim_HSRht3['deposit_Tartrazine'], bar_width, alpha=opacity, color=['g'])
            rects3 = axes[x].bar(index + 2*bar_width, df_sim_HSTremie12['deposit_Tartrazine'], bar_width, alpha=opacity, color=['b'])
            rects4 = axes[x].bar(index + 3*bar_width, df_sim_HSTremie13['deposit_Tartrazine'], bar_width, alpha=opacity, color=['m'])
            
            #Mise en forme
            axes[x].set_ylim(0, 6.5)
            axes[x].set_xlim(0, len(df_sim_HSMercia)) 
            axes[x].set_xticks(index+bar_width)
            #Label au dessus des colonnes
            df_sim_HSMercia['leaf_emergence'] = numpy.round(df_sim_HSMercia['leaf_emergence'], decimals=1)
            df_sim_HSRht3['leaf_emergence'] = numpy.round(df_sim_HSRht3['leaf_emergence'], decimals=1)
            df_sim_HSTremie12['leaf_emergence'] = numpy.round(df_sim_HSTremie12['leaf_emergence'], decimals=1)
            df_sim_HSTremie13['leaf_emergence'] = numpy.round(df_sim_HSTremie13['leaf_emergence'], decimals=1)
            df_sim_HSMercia['HS'] = numpy.round(df_sim_HSMercia['HS'], decimals=1)
            df_sim_HSRht3['HS'] = numpy.round(df_sim_HSRht3['HS'], decimals=1)
            df_sim_HSTremie12['HS'] = numpy.round(df_sim_HSTremie12['HS'], decimals=1)
            df_sim_HSTremie13['HS'] = numpy.round(df_sim_HSTremie13['HS'], decimals=1)
            df_sim_HSMercia['label'] = df_sim_HSMercia['HS'].astype(str)+' / '+df_sim_HSMercia['leaf_emergence'].astype(str)           
            df_sim_HSRht3['label'] = df_sim_HSRht3['HS'].astype(str)+' / '+df_sim_HSRht3['leaf_emergence'].astype(str)          
            df_sim_HSTremie12['label'] = df_sim_HSTremie12['HS'].astype(str)+' / '+df_sim_HSTremie12['leaf_emergence'].astype(str)           
            df_sim_HSTremie13['label'] = df_sim_HSTremie13['HS'].astype(str)+' / '+df_sim_HSTremie13['leaf_emergence'].astype(str)
            
            list1 =  df_sim_HSMercia['label']
            list2 =  df_sim_HSRht3['label']
            list3 =  df_sim_HSTremie12['label']
            list4 =  df_sim_HSTremie13['label']
            
            def autolabel(rects):
                # attach some text labels
                if rects == rects1 :
                    list=list1
                elif rects == rects2 :
                    list=list2
                elif rects == rects3 :
                    list=list3
                elif rects == rects4 :
                    list=list4
                for rect in rects:
                    height = rect.get_height()
                    axes[x].text(rect.get_x()+rect.get_width()/2., 1.05*height, list[rect.get_x()], ha='center', va='bottom', fontsize=8, rotation=90)
            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)
            autolabel(rects4)
            #Label x
            df_label = df_sim_HSMercia.append(df_sim_HSRht3)
            df_label = df_label.append(df_sim_HSTremie12)
            df_label = df_label.append(df_sim_HSTremie13)
            delta = df_label['delta'].dropna()
            axes[x].set_xticklabels( delta.tolist(), rotation=90, fontsize='small')
            if x == 0:
                axes[x].set_xlabel('delta_HS')
            axes[x].text(0.4, 6.2, 'treatment = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

            fig.suptitle('', fontsize=10)
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Mercia', 'Rht3', 'Tremie12', 'Tremie13'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
        
    if plot2==True:
        bar_width = 0.2; opacity = 0.4
        fig, axes = plt.subplots(nrows=1, ncols=2) 

        for x, date in [[0,'T1'], [1,'T2']]:
        
            df_sim = df_all[df_all['ntop_cur']==1]
            df_sim['deposit_Tartrazine/area'] = df_sim['deposit_Tartrazine'] / df_sim['area']
            
            df_sim_HSMercia = df_sim[df_sim['dim'] == date+'_Mercia']; df_sim_HSMercia = df_sim_HSMercia.reset_index()            
            
            df_sim_HSRht3 = df_sim[df_sim['dim'] == date+'_Rht3']; df_sim_HSRht3 = df_sim_HSRht3.sort(['var']); df_sim_HSRht3 = df_sim_HSRht3.reset_index()
            
            df_sim_HSTremie12 = df_sim[df_sim['dim'] == date+'_Tremie12']; df_sim_HSTremie12 = df_sim_HSTremie12.sort(['var']); df_sim_HSTremie12 = df_sim_HSTremie12.reset_index()
            
            df_sim_HSTremie13 = df_sim[df_sim['dim'] == date+'_Tremie13']; df_sim_HSTremie13 = df_sim_HSTremie13.sort(['var']); df_sim_HSTremie13 = df_sim_HSTremie13.reset_index()
            
            n_groups = len(df_sim_HSMercia)
            index = numpy.arange(n_groups) 

            rects1 = axes[x].bar(index, df_sim_HSMercia['deposit_Tartrazine/area'], bar_width, alpha=opacity, color=['r'])
            rects2 = axes[x].bar(index + bar_width, df_sim_HSRht3['deposit_Tartrazine/area'], bar_width, alpha=opacity, color=['g'])
            rects3 = axes[x].bar(index + 2*bar_width, df_sim_HSTremie12['deposit_Tartrazine/area'], bar_width, alpha=opacity, color=['b'])
            rects4 = axes[x].bar(index + 3*bar_width, df_sim_HSTremie13['deposit_Tartrazine/area'], bar_width, alpha=opacity, color=['m'])
            
            #Mise en forme
            axes[x].set_ylim(0, 0.6)
            axes[x].set_xlim(0, len(df_sim_HSMercia))        
            axes[x].set_xticks(index+bar_width)
            #Label au dessus des colonnes
            df_sim_HSMercia['leaf_emergence'] = numpy.round(df_sim_HSMercia['leaf_emergence'], decimals=1)
            df_sim_HSRht3['leaf_emergence'] = numpy.round(df_sim_HSRht3['leaf_emergence'], decimals=1)
            df_sim_HSTremie12['leaf_emergence'] = numpy.round(df_sim_HSTremie12['leaf_emergence'], decimals=1)
            df_sim_HSTremie13['leaf_emergence'] = numpy.round(df_sim_HSTremie13['leaf_emergence'], decimals=1)
            df_sim_HSMercia['HS'] = numpy.round(df_sim_HSMercia['HS'], decimals=1)
            df_sim_HSRht3['HS'] = numpy.round(df_sim_HSRht3['HS'], decimals=1)
            df_sim_HSTremie12['HS'] = numpy.round(df_sim_HSTremie12['HS'], decimals=1)
            df_sim_HSTremie13['HS'] = numpy.round(df_sim_HSTremie13['HS'], decimals=1)
            
            df_sim_HSMercia['label'] = df_sim_HSMercia['HS'].astype(str)+' / '+df_sim_HSMercia['leaf_emergence'].astype(str) 
            
            df_sim_HSRht3['label'] = df_sim_HSRht3['HS'].astype(str)+' / '+df_sim_HSRht3['leaf_emergence'].astype(str)          
            df_sim_HSTremie12['label'] = df_sim_HSTremie12['HS'].astype(str)+' / '+df_sim_HSTremie12['leaf_emergence'].astype(str)      
            
            df_sim_HSTremie13['label'] = df_sim_HSTremie13['HS'].astype(str)+' / '+df_sim_HSTremie13['leaf_emergence'].astype(str)
            
            list1 =  df_sim_HSMercia['label']
            list2 =  df_sim_HSRht3['label']
            list3 =  df_sim_HSTremie12['label']
            list4 =  df_sim_HSTremie13['label']
            
            def autolabel(rects):
                # attach some text labels
                if rects == rects1 :
                    list=list1
                elif rects == rects2 :
                    list=list2
                elif rects == rects3 :
                    list=list3
                elif rects == rects4 :
                    list=list4
                for rect in rects:
                    height = rect.get_height()
                    axes[x].text(rect.get_x()+rect.get_width()/2., 1.05*height, list[rect.get_x()], ha='center', va='bottom', fontsize=8, rotation=90)
            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)
            autolabel(rects4)
            #Label x
            df_label = df_sim_HSMercia.append(df_sim_HSRht3)
            df_label = df_label.append(df_sim_HSTremie12)
            df_label = df_label.append(df_sim_HSTremie13)
            delta = df_label['delta'].dropna()
            axes[x].set_xticklabels( delta.tolist(), rotation=90, fontsize='small')
                
            if x == 0:
                axes[x].set_xlabel('delta_HS')
            axes[x].text(0.2, 0.56, 'treatment = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

            fig.suptitle('', fontsize=10)
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Mercia', 'Rht3', 'Tremie12', 'Tremie13'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
      
    df_sim.to_csv('HS_com_synth_all.csv')
    return df_sim

#plot tartrazine/area par HS (HSnff-6 à HSnff+6) pour ntop=1 a 3
def simulation_efficacy(name='Mercia', hs=12, n_sim=5, n_plt=200, axis='MS'):

    def frange(x, y, jump):
        while x < y:
            yield x
            x += jump
    
    hs_lst = frange(hs-6, hs+6, ((float((hs+6)-(hs-6)))/20))
    lst = []
    for hs_ in hs_lst :
        lst.append([name,str(hs_)])

    if name=='Tremie12': #cas particulier de tremie12 n=30 
        n = 30
    else:
        n = n_plt

    df_sim = pandas.DataFrame()
    x=1
    while x <= n_sim:
        for var, stade in lst :          
            npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis=axis, to_csv=False)
            print 'var = '+var+' stade = '+stade+' - nbre de plantes sim = '+str(npl)
            dfmoy['var'] = var
            df_sim = df_sim.append(dfmoy)
        x+=1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'ntop']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr
    
def plot_efficacy(varieties=['Mercia','Rht3','Tremie12','Tremie13'], n_sim=1, n_plt=30, axis='MS', plot_tartrazine=True, plot_intercept=True, plot_cov=True, plot_protect=True, plot_global=True, csv=True):
    '''
    !!! Tous les plots sont declines en 2 versions :
    HS et HS-nff moyen dans une fourchette de -6 a +6 par pas de temps calcule d'environ 0.5
    plot tartrazine : simulation 'deposit_Tratrazine' 
    plot intercept : simulation 'interception efficacy'
    RAPPEL : interception efficacy = deposit_Tartrazine/area
    plot cov : simulation 'coverage efficacy'
    coverage efficacy = interception efficacy * exposition
    plot protect : simulation 'protection efficacy'
    protection efficacy = (1-lifetime) * coverage efficacy
    plot global : simulation 'somme protection efficacy'
    '''
    df_sim = pandas.DataFrame()
    for var in varieties :
        # nff moyen
        adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
        df_phenT = adel.phenT()
        nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
        df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
        mean = df_phenT.mean()
        hs_moyen = mean['nff']
        #sim 
        df_sim_var = simulation_efficacy(name=var, hs=hs_moyen, n_sim=n_sim, n_plt=n_plt, axis=axis)
        df_sim_var['deposit_Tartrazine/area'] = df_sim_var['deposit_Tartrazine'] / df_sim_var['area']
        df_sim_var['coverage_efficacy'] = df_sim_var['deposit_Tartrazine/area'] * df_sim_var['exposition']
        df_sim_var['protection_efficacy'] = (1-df_sim_var['lifetime']) * df_sim_var['coverage_efficacy']
        df_sim_var['HS-nffmean'] = df_sim_var['HS'] - hs_moyen
        df_sim = df_sim.append(df_sim_var)
    if csv==True:
        df_sim.to_csv('efficacy_data.csv')
        
    if plot_tartrazine==True:
        #x = hs
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_tartra = df_sim[df_sim['var']==var]
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_tartra[df_sim_tartra['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['deposit_Tartrazine'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['deposit_Tartrazine'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['deposit_Tartrazine'], '-.'+color, label=label)
            #fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage'); plt.ylabel('deposit_tartrazine')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
            
        # x = hs - nff moyen   
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_tartra2 = df_sim[df_sim['var']==var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_tartra2[df_sim_tartra2['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['deposit_Tartrazine'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['deposit_Tartrazine'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['deposit_Tartrazine'], '-.'+color, label=label)
            #fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.7), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage - nff moyen'); plt.ylabel('deposit_tartrazine')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
        
    if plot_intercept==True: 
        #x = hs
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_intercept = df_sim[df_sim['var']==var]
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_intercept[df_sim_intercept['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['deposit_Tartrazine/area'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['deposit_Tartrazine/area'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['deposit_Tartrazine/area'], '-.'+color, label=label)
            #fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.03), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage'); plt.ylabel('interception efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
        
        #x = hs - nffmoyen
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_intercept2 = df_sim[df_sim['var']==var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_intercept2[df_sim_intercept2['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['deposit_Tartrazine/area'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['deposit_Tartrazine/area'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['deposit_Tartrazine/area'], '-.'+color, label=label)
            #fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage - nff moyen'); plt.ylabel('interception efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
    
    if plot_cov==True: 
        # x = hs
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_cov = df_sim[df_sim['var']==var]
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_cov[df_sim_cov['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['coverage_efficacy'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['coverage_efficacy'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['coverage_efficacy'], '-.'+color, label=label)
            #fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage'); plt.ylabel('coverage_efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
        
        # x = hs - nff moyen
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_cov2 = df_sim[df_sim['var']==var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_cov2[df_sim_cov2['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['coverage_efficacy'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['coverage_efficacy'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['coverage_efficacy'], '-.'+color, label=label)
            #fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage - nff moyen'); plt.ylabel('coverage_efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

    if plot_protect==True:
        # x = hs
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_protect = df_sim[df_sim['var']==var]
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_protect[df_sim_protect['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['protection_efficacy'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['protection_efficacy'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS'], df_sim_ntop['protection_efficacy'], '-.'+color, label=label)
            #fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage'); plt.ylabel('protection_efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
        
        # x = hs - nff moyen
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_protect2 = df_sim[df_sim['var']==var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            #plot
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_protect2[df_sim_protect2['ntop']==ntop]
                #gestion de la legende
                if var=='Mercia':
                    label = 'true F'+str(ntop)
                else:
                    label='_nolegend_'
                #plot par ntop
                if ntop==1:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['protection_efficacy'], '-'+color, label=label)
                elif ntop==2:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['protection_efficacy'], '--'+color, label=label)
                elif ntop==3:
                    plt.plot(df_sim_ntop['HS-nffmean'], df_sim_ntop['protection_efficacy'], '-.'+color, label=label)
            #fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage - nff moyen'); plt.ylabel('protection_efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
               
    if plot_global==True:
        # x = hs
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_protect = df_sim[df_sim['var']==var]
            #traitement dataframe
            df_sim_protect_all = pandas.DataFrame()
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_protect[df_sim_protect['ntop']==ntop]
                df_sim_protect_all = df_sim_protect_all.append(df_sim_ntop)
            df_sim_protect_all = df_sim_protect_all.groupby(['HS']).sum()  
            df_sim_protect_all = df_sim_protect_all.reset_index()#gestion de la legende
            if var=='Mercia':
                label = 'Somme des protections efficacy'
            else:
                label='_nolegend_'
            #plot
            plt.plot(df_sim_protect_all['HS'], df_sim_protect_all['protection_efficacy'], '-'+color, label=label)
            #fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage'); plt.ylabel('global_efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
        
        # x = hs - nff moyen
        plt.figure()
        for var in varieties:        
            if var=='Mercia':
                color='r'
            elif var=='Rht3':
                color='g'
            elif var=='Tremie12':
                color='b'
            elif var=='Tremie13':
                color='m'
            df_sim_protect = df_sim[df_sim['var']==var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            #traitement dataframe
            df_sim_protect_all = pandas.DataFrame()
            for ntop in [1,2,3]:
                df_sim_ntop = df_sim_protect[df_sim_protect['ntop']==ntop]
                df_sim_protect_all = df_sim_protect_all.append(df_sim_ntop)
            df_sim_protect_all = df_sim_protect_all.groupby(['HS-nffmean']).sum()  
            df_sim_protect_all = df_sim_protect_all.reset_index()#gestion de la legende
            if var=='Mercia':
                label = 'Somme des protections efficacy'
            else:
                label='_nolegend_'
            #plot
            plt.plot(df_sim_protect_all['HS-nffmean'], df_sim_protect_all['protection_efficacy'], '-'+color, label=label)
            #fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03), arrowprops=dict(color=color, arrowstyle="->", connectionstyle="arc3"))
        #mise en forme
        plt.xlabel('haun stage - nff moyen'); plt.ylabel('global_efficacy')
        plt.grid(False) # enlever la grille
        #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
                
    return df_sim
    
#plot data publi contre simulation au nff moyen    
def simulation_nffmoyen(var='Mercia', n_sim=5, n_plt=200, hs=12, axis='MS', csv=True): #axis = MS or all          
    if var=='Tremie12':
        nplants_lst = 30
    else :
        nplants_lst = n_plt
        
    x=1
    df_sim = pandas.DataFrame() 
    while x <= n_sim: 
        npl, dfmoy, dfsd = treatment(name=var, sim=str(hs), nplants=n_plt, axis=axis, to_csv=False)
        print 'var = '+var+' stade = '+str(hs)+' - nplants = '+str(npl)
        dfmoy['var'] = var
        df_sim = df_sim.append(dfmoy)
        x+=1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    if csv==True:
        df_sim.to_csv('sim_nffmoyen_all_'+var+'.csv')
        df_sim_gr.to_csv('sim_nffmoyen_synth_'+var+'.csv')
    return df_sim_gr  
    
def plot_data_simnffmoyen(varieties=['Mercia','Rht3','Tremie12','Tremie13'], n_sim=1, n_plt=30, axis='MS'):
    '''Plot1 = obs
    Plot2 = obs/area'''
    df_sim = pandas.DataFrame()
    #obs + preparation de obs
    df_obs = idata.dye_interception() 
    df_obs.rename(columns={'name':'var'}, inplace=True) 
    df_obs.rename(columns={'N feuille':'ntop_cur'}, inplace=True)
    df_obs.rename(columns={'HS':'HS_obs'}, inplace=True)
    #fichier de simulation
    for var in varieties :
        # nff moyen
        if var=='Mercia':
            delta = -0.66
        elif var=='Rht3':
            delta = -1.2
        elif var=='Tremie12':
            delta = -0.5
        elif var=='Tremie13':
            delta = -0.7
        adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
        df_phenT = adel.phenT()
        nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
        df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
        mean = df_phenT.mean()
        hs_moyen = mean['nff']
        hs_moyen = hs_moyen + delta
        #sim 
        df_sim_var = simulation_nffmoyen(var=var, hs=hs_moyen, n_sim=n_sim, n_plt=n_plt, axis=axis)
        df_sim_var['deposit_Tartrazine/area'] = df_sim_var['deposit_Tartrazine'] / df_sim_var['area']
        #df_sim_var['HS-nffmean'] = df_sim_var['HS'] - hs_moyen
        df_sim = df_sim.append(df_sim_var)
    #merge des 2 tableaux obs et sim
    df_obs = df_obs[df_obs['treatment']=='T2'] 
    df_all = df_sim.merge(df_obs)
    #column tartrazine/area
    df_all['meanObs/area'] = df_all['mean']/df_all['area']
    
    #plot
    bar_width = 0.3; opacity = 0.4
    val=[[0, 'mean', 'deposit_Tartrazine'], [1, 'meanObs/area', 'deposit_Tartrazine/area']]
    fig, axes = plt.subplots(nrows=1, ncols=2) 

    for x, col_data, col_sim in val:   
        #on ne veut que les 4 feuilles du haut
        df_all = df_all[df_all['ntop_cur']<=4] 
        
        n_groups = len(df_all)
        index = numpy.arange(n_groups) 

        #publi
        if x==0:
            rects1 = axes[x].bar(index, df_all[col_data], bar_width, alpha=opacity, color='y', yerr=df_all['IC'], error_kw=dict(ecolor='k'))
        else:
            rects1 = axes[x].bar(index, df_all[col_data], bar_width, alpha=opacity, color='y')
        #sim a hs nff moyen + delta
        rects2 = axes[x].bar(index + bar_width, df_all[col_sim], bar_width, alpha=opacity, color=['r','r','r','r','g','g','g','g','b','b','b','b','m','m','m','m'])
       
        # Mise en forme
        df_all['HS'] = numpy.round(df_all['HS'], decimals=1)
        label = ['','Mercia','','Rht3','','Tremie12','','Tremie13']
        axes[x].set_xticklabels(label, rotation=90, fontsize='small')
        if x == 0:
            axes[x].set_xlabel('var')
            axes[x].set_ylabel('deposit_tartrazine')
        else:
            axes[x].set_ylabel('deposit_tartrazine/area')
            
        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        
        def autolabel(rects):
            # attach some text labels
            if rects == rects2 :
                list=['F1','F2','F3','F4']
            else:
                print 'pb'
            for rect in rects:
                height = rect.get_height()
                if int(rect.get_x()) <4 :
                    if x==0:
                        print int(rect.get_x())
                        axes[x].text(rect.get_x()+rect.get_width()/2. -0.1, 1.05*height + 0.2, list[int(rect.get_x())], ha='center', va='bottom', fontsize=8)
        autolabel(rects2)
        
        
        '''axes[x].text(0.2, 6.1, ''+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
               
    axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('F1', 'F2', 'F3', 'F4'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})'''
    
