import os
import pandas
import numpy
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

def repartition_at_application(name, appdate='T1', dose=1e4, nplants=30, density=1):
    """ 10000 l ha-1 is for 1 l/m2 """
    
    adel = Reconstructions.get_reconstruction(name=name, nplants=nplants, adjust_density = {name:None} , stand_density_factor = {name:density} )
    date, hs = HS_applications[name][appdate]
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
    return adel, df

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
                   'deposit_Tartrazine': numpy.sum
                   })
    
def treatment(name='Tremie12', sim='T1', nplants=200, axis='MS', to_csv=False, density=1): 

    adel, df = repartition_at_application(name,sim,nplants=nplants,density=density)
    
    if axis == 'MS':
        df = df[df['axe'] == 'MS']
    
    data = aggregate_by_leaf(df)
        
    # compute ntop_cur
    gr = data.groupby(['HS','plant','axe'], group_keys=False)
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
    data = gr.apply(_fun)
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1

    # aggregation by ntop cur
    sub = data.ix[:,('HS','TT','axe','metamer','ntop','ntop_cur','length','area','green_area','deposit_Tartrazine')]
    dfmoy = sub.groupby(['HS','axe','ntop_cur']).mean().reset_index()
    dfsd = sub.groupby(['HS','axe','ntop_cur']).std().reset_index()

    if to_csv:
        data.to_csv(''.join(('data_',name,'_',sim,'_',axis,'.csv')))
        dfmoy.to_csv(''.join(('moy_',name,'_',sim,'_',axis,'.csv')))
        dfsd.to_csv(''.join(('std_',name,'_',sim,'_',axis,'.csv')))
             
    return adel.nplants, dfmoy, dfsd
   
#lancement simulation et concatenation des sorties pour obtenir un fichier de donnees simulees
def verif_sim(name='Mercia'):
    lst = [[name,'T1'], [name,'T2']]
    df_sim = pandas.DataFrame()
    for var, stade in lst :
        for nbr_plt in [30,60,100,200]:
            x = 1
            while x<=5 :
                print 'ETAPE = '+var+' '+stade+' '+str(nbr_plt)+' '+str(x)
                npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=nbr_plt, axis='MS', to_csv=False)
                dfmoy['var'] = var
                dfmoy['nbre plt simulees'] = npl
                df_sim = df_sim.append(dfmoy)
                x += 1
    df_sim.to_csv('verif_sim.csv')
    return df_sim

 
def simulation_dimension(name='Tremie12', stade_dim='div2', csv=True): 
    # pour faire graph T2 de Tremie12 contre T2 de Tremie13, moyenne de 5 simulations
    # pas de param dim => à changer directement dans architectural_data.py
    # stade_dim = div2 div1.5 normal mult1.5 mult2
    lst = [['Tremie13','T1'],['Tremie13','T2']]
    if name=='Tremie12':
        n = 30
    else:
        n = 200
    df_sim = pandas.DataFrame()
    x=0
    while x<5:
        for var, stade in lst :          
            npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis='MS', to_csv=False)
            #print 'var = '+var+' stade = '+stade+' - nbre de plantes sim = '+str(npl)
            dfmoy['var'] = var; dfmoy['dim'] = stade_dim
            df_sim = df_sim.append(dfmoy)
        x+=1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'dim', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    if csv == True:
        df_sim.to_csv('dim_all_'+stade_dim+'_'+name+'.csv')
        df_sim.to_csv('dim_synth_'+stade_dim+'_'+name+'.csv')
    return df_sim_gr
    
def plot_dimension(plot1=False, plot2=True, plot3=False, varieties=['Tremie13']):
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
    # ne marche pas pour le moment -> obligation de lancer simulation_dimension() autant de fois que de dimension differente pour chacune des varietes (a modifier dans architectural_data) puis concatener l'ensemble des sorties
    #for n in varieties:
    #    df_sim_var = simulation_dimension(name=n)
    #    df_sim = df_sim.append(df_sim_var)
    df_sim = 'dim_tremie13.csv'
    df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')
    
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
                df_fin.ix[df_fin.dim.isin(['norm']), 'index_dim'] = 3
                df_fin.ix[df_fin.dim.isin(['mult1.5']), 'index_dim'] = 4
                df_fin.ix[df_fin.dim.isin(['mult2']), 'index_dim'] = 5
                #rects2 pour ne tracer les obs que quand dim = normale
                df_fin['mean_dim'] = df_fin['mean']
                df_fin.ix[df_fin.dim.isin(['/2']), 'mean_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['/1.5']), 'mean_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['mult1.5']), 'mean_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['mult2']), 'mean_dim'] = 0
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
                if var=='Tremie12':
                    axes[x].text(0.4, 16.2, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                else:
                    axes[x].text(0.4, 16.2, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

                fig.suptitle('Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur', fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    # PLOT 2
    if plot2 is True :
        for var in varieties:
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
                df_fin.ix[df_fin.dim.isin(['norm']), 'index_dim'] = 3
                df_fin.ix[df_fin.dim.isin(['mult1.5']), 'index_dim'] = 4
                df_fin.ix[df_fin.dim.isin(['mult2']), 'index_dim'] = 5
                #rects2 
                df_fin['Area A_bl_dim'] = df_fin['Area A_bl']
                df_fin.ix[df_fin.dim.isin(['/2']), 'Area A_bl_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['/1.5']), 'Area A_bl_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['mult1.5']), 'Area A_bl_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['mult2']), 'Area A_bl_dim'] = 0
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
                if var=='Tremie12':
                    axes[x].text(0.4, 48, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                else :
                    axes[x].text(0.4, 48, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
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
                df_fin.ix[df_fin.dim.isin(['norm']), 'index_dim'] = 3
                df_fin.ix[df_fin.dim.isin(['mult1.5']), 'index_dim'] = 4
                df_fin.ix[df_fin.dim.isin(['mult2']), 'index_dim'] = 5
                #rects2
                df_fin['mean_dim'] = df_fin['mean/area']
                df_fin.ix[df_fin.dim.isin(['/2']), 'mean_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['/1.5']), 'mean_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['mult1.5']), 'mean_dim'] = 0
                df_fin.ix[df_fin.dim.isin(['mult2']), 'mean_dim'] = 0
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
                if var=='Tremie12':
                    axes[x].text(0.4, 0.33, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                else:
                    axes[x].text(0.4, 0.33, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                fig.suptitle('Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur', fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )

    return df_scan, df_obs, df_sim
    
def simulation_density(name='Tremie12', csv=True):
    # test date interception -0.4 -0.2 +0.2 +0.4 pour chaque date, +0.5 et +1 pour T2 seulement - !!! PB TREMIE12 SIM N=200 !!!
    if name=='Tremie12':
        lst = [['Tremie12','T1-0.4'], ['Tremie12','T1-0.2'], ['Tremie12','T1'], ['Tremie12','T1+0.2'], ['Tremie12','T1+0.4'],
            ['Tremie12','T2-0.4'], ['Tremie12','T2-0.2'], ['Tremie12','T2'], ['Tremie12','T2+0.2'], ['Tremie12','T2+0.4'], ['Tremie12','T2+0.5'], ['Tremie12','T2+1'], ['Tremie12','T2+1.5'],['Tremie12','T2+2'],['Tremie12','T2+2.5']]
        n = 30
    else:
        lst = [['Tremie13','T1-0.4'], ['Tremie13','T1-0.2'], ['Tremie13','T1'], ['Tremie13','T1+0.2'], ['Tremie13','T1+0.4'],
            ['Tremie13','T2-0.4'], ['Tremie13','T2-0.2'], ['Tremie13','T2'], ['Tremie13','T2+0.2'], ['Tremie13','T2+0.4'], ['Tremie13','T2+0.5'], ['Tremie13','T2+1'], ['Tremie13','T2+1.5'],['Tremie13','T2+2'],['Tremie13','T2+2.5']]
        n = 200

    df_sim = pandas.DataFrame()
    x=0
    while x<5:
        for var, stade in lst :  
            for dens in [0.5,0.625,0.75,0.875,1]:         
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
    
def plot_density(plot1=False, plot2=True, plot3=False, varieties=['Tremie12','Tremie13']):
    '''
    Plot 1 =  Barplot observation / simulation avec ntop_cur en abscisse et (dfmoy['deposit_Tartrazin']  versus obs 'moyenne' par genotype et par date interception
    Plot 2 = Barplot de controle 'area' versus scan data (avec la colone area de dfmoy)
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
        df_sim_var = simulation_density(name=n)
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
            '''if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]'''
            if var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date] 
                #mean_dens pour dessiner les obs seulement qd densite=1
                df_fin['mean_dens'] = df_fin['mean']
                df_fin.ix[df_fin.density.isin([0.5]), 'mean_dens'] = 0
                df_fin.ix[df_fin.density.isin([0.625]), 'mean_dens'] = 0
                df_fin.ix[df_fin.density.isin([0.74]), 'mean_dens'] = 0 # cas particulier de Tremie12 qui devrait disparaitre
                df_fin.ix[df_fin.density.isin([0.75]), 'mean_dens'] = 0
                df_fin.ix[df_fin.density.isin([0.875]), 'mean_dens'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur','HS','density'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups) 

                rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
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
                if var=='Tremie12':
                    axes[x].text(0.4, 7.4, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                else :
                    axes[x].text(0.4, 7.4, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

                fig.suptitle('Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    # PLOT 2
    if plot2 is True :
        for var in varieties:
            df_scan_var = df_scan[df_scan['var']==var]
            df_sim_var = df_sim[df_sim['var']==var]
            df_scan_var.HS=map(str,df_scan_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_scan_var.ix[:,['HS','ntop_cur','Area A_bl']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=5]
            #plot
            bar_width = 0.4; opacity = 0.4
            fig, axes = plt.subplots(nrows=1, ncols=2) 
            '''if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]'''
            if var == 'Tremie12' :
                val = [[0,10.98],[1,12.63]]
            else:
                val = [[0,8.7],[1,11.04]]
            for x, date in val:
                df_all.HS=map(float,df_all.HS)
                df_fin = df_all[df_all['HS']==date]
                #Area A_bl dens pour dessiner les obs seulement qd densite=1
                df_fin['Area A_bl dens'] = df_fin['Area A_bl']
                df_fin.ix[df_fin.density.isin([0.5]), 'Area A_bl dens'] = 0
                df_fin.ix[df_fin.density.isin([0.625]), 'Area A_bl dens'] = 0
                df_fin.ix[df_fin.density.isin([0.74]), 'Area A_bl dens'] = 0 # cas particulier de Tremie12 qui devrait disparaitre
                df_fin.ix[df_fin.density.isin([0.75]), 'Area A_bl dens'] = 0
                df_fin.ix[df_fin.density.isin([0.875]), 'Area A_bl dens'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur', 'HS','density'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups) 

                rects1 = axes[x].bar(index, df_fin['Area A_bl dens'], bar_width, alpha=opacity, color='y')
                rects2 = axes[x].bar(index + bar_width, df_fin['area'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
                
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
                if var=='Tremie12':
                    axes[x].text(0.4, 28.5, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                else :
                    axes[x].text(0.4, 28.5, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                fig.suptitle('Surface scan/sim par ntop_cur pour '+var, fontsize=10)
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
                df_fin.ix[df_fin.density.isin([0.5]), 'mean/area_dens'] = 0
                df_fin.ix[df_fin.density.isin([0.625]), 'mean/area_dens'] = 0
                df_fin.ix[df_fin.density.isin([0.74]), 'mean/area_dens'] = 0 # cas particulier de Tremie12 qui devrait disparaitre
                df_fin.ix[df_fin.density.isin([0.75]), 'mean/area_dens'] = 0
                df_fin.ix[df_fin.density.isin([0.875]), 'mean/area_dens'] = 0
                #---
                df_fin = df_fin.sort(['ntop_cur','HS','density'])      
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)   
                
                rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
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
                if var=='Tremie12':
                    axes[x].text(0.4, 0.45, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                else :
                    axes[x].text(0.4, 0.45, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)
                fig.suptitle('Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )
    '''           
    # PLOT 4
    if plot4 is True :
        for var in ['Mercia', 'Rht3', 'Tremie13']:
            df_obs_var = df_obs[df_obs['name']==var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var']==var]
            df_obs_var.HS=map(str,df_obs_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_obs_var.ix[:,['HS','ntop_cur','mean']], how='outer')
            df_all = df_all[df_all['ntop_cur']<=5]
            df_all['tartrazine/area'] = df_all['deposit_Tartrazine']/df_all['area']
            df_all['mean/area'] = df_all['mean']/df_all['area']
            if var=='Mercia' or var=='Rht3':
                TTlig = 1383.2
            elif var=='Tremie13':
                TTlig = 1230.21
            df_all['lig'] = df_all['TT'] - TTlig
            df_all['lig'] = numpy.round(df_all['lig'], decimals=2)
            df_all['ntop'] = numpy.round(df_all['ntop'], decimals=2)
            #plot
            n_groups = len(df_all['lig'])
            index = numpy.arange(n_groups)
            bar_width = 0.35; opacity = 0.4
            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2)
            # x = age feuille en dd depuis ligulation
            df_all = df_all.sort(['lig','ntop_cur'])
            rects1 = ax0.bar(index, df_all['tartrazine/area'], bar_width, alpha=opacity, color='b')
            rects2 = ax0.bar(index + bar_width, df_all['mean/area'], bar_width, alpha=opacity, color='r')
            ax0.set_xticks(index+bar_width)
            ax0.set_xticklabels( df_all['lig'].tolist(), rotation=90, fontsize='x-small' )
            ax0.set_xlabel('TTdepuisLig')
            # x = ntop (nb feuille depuis le bas)
            df_all = df_all.sort(['ntop'])
            rects1 = ax1.bar(index, df_all['tartrazine/area'], bar_width, alpha=opacity, color='b')
            rects2 = ax1.bar(index + bar_width, df_all['mean/area'], bar_width, alpha=opacity, color='r')
            ax1.set_xticks(index+bar_width)
            ax1.set_xticklabels( df_all['ntop'].tolist(), rotation=90, fontsize='x-small' )
            ax1.legend((rects1[0], rects2[0]), ('Sim', 'Obs'), bbox_to_anchor=[1.10, 1.12] )
            ax1.set_xlabel('ntop')
            fig.suptitle('Sim [deposit_Tartrazine/area] / Obs [mean/area] pour '+var, fontsize=10)
    '''
    return df_scan, df_obs, df_sim

    
    
    '''       
    lst_all = [['Mercia','T1-0.4'], ['Mercia','T1-0.2'], ['Mercia','T1'], ['Mercia','T1+0.2'], ['Mercia','T1+0.4'],
            ['Mercia','T2-0.4'], ['Mercia','T2-0.2'], ['Mercia','T2'], ['Mercia','T2+0.2'], ['Mercia','T2+0.4'], ['Mercia','T2+0.5'], ['Mercia','T2+1'],
            ['Rht3','T1-0.4'], ['Rht3','T1-0.2'], ['Rht3','T1'], ['Rht3','T1+0.2'], ['Rht3','T1+0.4'],
            ['Rht3','T2-0.4'], ['Rht3','T2-0.2'], ['Rht3','T2'], ['Rht3','T2+0.2'], ['Rht3','T2+0.4'], ['Rht3','T2+0.5'], ['Rht3','T2+1']]
            #['Tremie12','date1-0.4'], ['Tremie12','date1-0.2'], ['Tremie12','date1'], ['Tremie12','date1+0.2'], ['Tremie12','date1+0.4'],
            #['Tremie12','date2-0.4'], ['Tremie12','date2-0.2'], ['Tremie12','date2'], ['Tremie12','date2+0.2'], ['Tremie12','date2+0.4'],
            #['Tremie12','T1-0.4'], ['Tremie12','T1-0.2'], ['Tremie12','T1'], ['Tremie12','T1+0.2'], ['Tremie12','T1+0.4'],
            #['Tremie12','T2-0.4'], ['Tremie12','T2-0.2'],  ['Tremie12','T2'],  ['Tremie12','T2+0.2'], ['Tremie12','T2+0.4'], ['Tremie12','T2+0.5'], ['Tremie12','T2+1'],
            #['Tremie13','date1-0.4'], ['Tremie13','date1-0.2'], ['Tremie13','date1'], ['Tremie13','date1+0.2'], ['Tremie13','date1+0.4'],
            ['Tremie13','T1-0.4'], ['Tremie13','T1-0.2'], ['Tremie13','T1'], ['Tremie13','T1+0.2'], ['Tremie13','T1+0.4'],
            #['Tremie13','date2-0.4'], ['Tremie13','date2-0.2'], ['Tremie13','date2'], ['Tremie13','date2+0.2'], ['Tremie13','date2+0.4'],
            ['Tremie13','T2-0.4'], ['Tremie13','T2-0.2'], ['Tremie13','T2'], ['Tremie13','T2+0.2'], ['Tremie13','T2+0.4'], ['Tremie13','T2+0.5'], ['Tremie13','T2+1']]                   
    '''      

def simulation_HS(name='Tremie12', csv=True):
    # test date interception -0.4 -0.2 +0.2 +0.4 pour chaque date, +0.5 et +1 pour T2 seulement - !!! PB TREMIE12 SIM N=200 !!!
    lst = [['Tremie13','T1'],['Tremie13','T2']]
    if name=='Tremie12':
        lst = [['Tremie12','date1-0.4'], ['Tremie12','date1-0.2'], ['Tremie12','date1'], ['Tremie12','date1+0.2'], ['Tremie12','date1+0.4'],
            ['Tremie12','date2-0.4'], ['Tremie12','date2-0.2'], ['Tremie12','date2'], ['Tremie12','date2+0.2'], ['Tremie12','date2+0.4'],
            ['Tremie12','T1-0.4'], ['Tremie12','T1-0.2'], ['Tremie12','T1'], ['Tremie12','T1+0.2'], ['Tremie12','T1+0.4'],
            ['Tremie12','T2-0.4'], ['Tremie12','T2-0.2'],  ['Tremie12','T2'],  ['Tremie12','T2+0.2'], ['Tremie12','T2+0.4'], ['Tremie12','T2+0.5'], ['Tremie12','T2+1'],['Tremie12','T2+1.5'], ['Tremie12','T2+2'], ['Tremie12','T2+2.5']]
        n=30
    else:
        lst = [['Tremie13','date1-0.4'], ['Tremie13','date1-0.2'], ['Tremie13','date1'], ['Tremie13','date1+0.2'], ['Tremie13','date1+0.4'],
            ['Tremie13','date2-0.4'], ['Tremie13','date2-0.2'], ['Tremie13','date2'], ['Tremie13','date2+0.2'], ['Tremie13','date2+0.4'],
            ['Tremie13','T1-0.4'], ['Tremie13','T1-0.2'], ['Tremie13','T1'], ['Tremie13','T1+0.2'], ['Tremie13','T1+0.4'],
            ['Tremie13','T2-0.4'], ['Tremie13','T2-0.2'],  ['Tremie13','T2'],  ['Tremie13','T2+0.2'], ['Tremie13','T2+0.4'], ['Tremie13','T2+0.5'], ['Tremie13','T2+1'],['Tremie13','T2+1.5'], ['Tremie13','T2+2'], ['Tremie13','T2+2.5']]
        n=200
    df_sim = pandas.DataFrame()
    x=0
    while x<5:
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
    
def plot_HS(plot1=False, plot2=True, plot3=False, varieties=['Tremie12','Tremie13']):
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
        df_sim_var = simulation_HS(name=n)
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
            '''if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]'''
            if var == 'Tremie12' :
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
                    for f in [-0.4,-0.2,0.2,0.4,0.5,1,1.5,2,2.5]:
                        df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                df_fin = df_fin.sort(['ntop_cur','HS'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups) 

                if date == 8.7 or date == 10.98:
                    rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
                else :
                    rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','c','c','c','m','m','m','m','m','m','m','m','m','m','g','g','g','g','g','g','g','g','g','g','r','r','r','r','r','r','r','r','r','r','b','b','b','b','b','b','b','b','b','b'])
   
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
                if var=='Tremie12':
                    axes[x].text(0.4, 6.5, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                else:
                    axes[x].text(0.4, 6.5, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                fig.suptitle('Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[10], rects1[20], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})
    
    # PLOT 2
    if plot2 is True :
        for var in varieties:
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
                    rects2 = axes[x].bar(index + bar_width, df_fin['area'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','c','c','c','m','m','m','m','m','m','m','m','m','m','g','g','g','g','g','g','g','g','g','g','r','r','r','r','r','r','r','r','r','r','b','b','b','b','b','b','b','b','b','b'])
                
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
                if var=='Tremie12':
                    axes[x].text(0.4, 28.5, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                else:
                    axes[x].text(0.4, 28.5, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                fig.suptitle('Surface scan/sim par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects2[0], rects2[10], rects2[20], rects1[0]), ( 'Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )
            
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
            '''if var=='Mercia':
                val = [[0,9.74],[1,12.8]]
            elif var == 'Rht3' :
                val = [[0,9.15],[1,12.48]]'''
            if var == 'Tremie12':
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
                    for f in [-0.4,-0.2,0.2,0.4,0.5,1,1.5,2,2.5]:
                        df_fin = df_fin.append(df_all[df_all['HS']==round(date+f,2)])
                df_fin = df_fin.sort(['ntop_cur','HS'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)   
                
                #rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','m','g','r','b'])
                if date == 10.98 or date == 8.7:
                    rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','c','c','c','c','m','m','m','m','m','g','g','g','g','g','r','r','r','r','r','b','b','b','b','b'])
                else :
                    rects1 = axes[x].bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','c','c','c','c','c','c','c','c','c','m','m','m','m','m','m','m','m','m','m','g','g','g','g','g','g','g','g','g','g','r','r','r','r','r','r','r','r','r','r','b','b','b','b','b','b','b','b','b','b'])
                    
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
                if var=='Tremie12':
                    axes[x].text(0.4, 0.32, 'TREMIE 12 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                else:
                    axes[x].text(0.4, 0.32, 'TREMIE 13 - HS = '+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=14)
                fig.suptitle('Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur pour '+var, fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[10], rects1[20], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )
            
    return df_scan, df_obs, df_sim
    
'''
#old - a supprimer?
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
'''
