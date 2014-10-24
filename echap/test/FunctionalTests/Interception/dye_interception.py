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
    do_record(g, application_data, recorder, header={'TT':age, 'HS':hs})
    df =  recorder.get_records()
    print 'repartition_at_application df.columns before ', df.columns
    return df

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
    
def treatment(name='Tremie12', sim='T1', nplants=30, axis='MS', to_csv=False): 

    df = repartition_at_application(name,sim,nplants=nplants)
    
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
             
    return dfmoy, dfsd
    
#lancement simulation et concatenation des sorties pour obtenir un fichier de donnees simulees
def simulation():
    #lst = [['Mercia','T1'], ['Mercia','T2'], ['Rht3','T1'], ['Rht3','T2'], ['Tremie12','T1'], ['Tremie12','T2'], ['Tremie13','T1'], ['Tremie13','T2']]
    lst = [['Tremie12','date1'], ['Tremie12','date2'], ['Tremie12','T1'], ['Tremie12','T2'], ['Tremie13','date1'], ['Tremie13','T1'],['Tremie13','date2'], ['Tremie13','T2']]
    df_sim = pandas.DataFrame()
    for var, stade in lst :
        dfmoy, dfsd = treatment(name=var, sim=stade, nplants=30, axis='MS', to_csv=False)
        dfmoy['var'] = var
        df_sim = df_sim.append(dfmoy)
    return df_sim

#plot des donnees obs
def plot(plot1=False, plot2=True):
    '''
    Plot 1 = 
    Plot 2 = Plot de controle 'area' versus scan data (avec la colone area de dfmoy)
    '''
    
    #base
    df_scan_tremie12 = archidb.treatment_scan('Tremie12'); df_scan_tremie12['var']='Tremie12'
    df_scan_tremie13 = archidb.treatment_scan('Tremie13'); df_scan_tremie13['var']='Tremie13'
    df_scan = df_scan_tremie12.append(df_scan_tremie13)
    df_obs = idata.dye_interception()
    df_sim = simulation()
    
    # PLOT 2
    if plot2 is True :
        for var in ['Tremie12', 'Tremie13']:
            df_scan_var = df_scan[df_scan['var']==var]
            df_sim_var = df_sim[df_sim['var']==var]
            df_scan_var.HS=map(str,df_scan_var.HS)
            df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(df_scan_var.ix[:,['HS','ntop_cur','Area A_bl']])
            print df_all
            #plot
            bar_width = 0.4; opacity = 0.4
            fig, axes = plt.subplots(nrows=1, ncols=len(df_scan_var['HS'].unique()))
            val = df_scan_var['HS'].unique()
            for x, HS in enumerate(val):
                print HS
                df_fin = df_all[df_all['HS']==HS]
                n_groups = len(df_fin['ntop_cur'].unique())
                index = numpy.arange(n_groups)      
                rects1 = axes[x].bar(index, df_fin['Area A_bl'], bar_width,
                             alpha=opacity,
                             color='b')
                rects2 = axes[x].bar(index + bar_width, df_fin['area'], bar_width,
                             alpha=opacity,
                             color='r')   
                # Mise en forme
                axes[x].set_ylim(0, 35)
                axes[x].set_xlim(0, df_all['ntop_cur'].max())          
                axes[x].set_xticks(index+bar_width)
                axes[x].set_xticklabels( ('1', '2', '3', '4', '5', '6') )
                if x == 0:
                    axes[x].set_xlabel('ntop_cur'); axes[x].set_ylabel('area (cm2)')
                axes[x].text(0.4, 33, 'HS : '+str(HS), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10})
                fig.suptitle('Surface scan/sim pour '+var, fontsize=10)
            axes[x].legend((rects1[0], rects2[0]), ('Scan', 'Sim') )
    
    return df_scan, df_obs, df_sim

'''
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
