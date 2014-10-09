"""
 Ploting function companion module for architectural reconstructions
"""
import matplotlib.pyplot as plt
import pandas
import numpy as np

def dimension_plot(dimension_data, fits):
    plt.ion()
    
    fig, axes = plt.subplots(nrows=2, ncols=3)
    ax0, ax1, ax2, ax3, ax4, ax5 = axes.flat
    
    varieties = [['Mercia','r',12],['Rht3','g',11],['Tremie12','b',13],['Tremie13','m',11]]
    for name,color,nff in varieties :
        d_data = dimension_data[name]
        dim_data = d_data.where((pandas.notnull(d_data)), None)
        ax0.plot(dim_data['index_phytomer'], dim_data['L_blade'], '-^'+color, label = '_nolegend_')
        ax1.plot(dim_data['index_phytomer'], dim_data['L_sheath'], '--o'+color, label = '_nolegend_')
        ax2.plot(dim_data['index_phytomer'], dim_data['L_internode'], ':<'+color, label=name)
        #ax0.plot(dim_data['index_phytomer'], dim_data['W_blade'], style='--p'+color, label='W_blade ' +name)
        
        dim_fit = fits[name][nff]
        ax3.plot(dim_fit['index_phytomer'], dim_fit['L_blade'], '-^'+color, label = '_nolegend_')
        ax4.plot(dim_fit['index_phytomer'], dim_fit['L_sheath'], '--o'+color, label = '_nolegend_')
        ax5.plot(dim_fit['index_phytomer'], dim_fit['L_internode'], ':<'+color, label = '_nolegend_')
        #ax1.plot(dim_fit['index_phytomer'], dim_fit['W_blade'], style='--p'+color, label='W_blade ' +name)
        
        ax0.set_title('L_blade', fontsize=10)
        ax1.set_title('L_sheath', fontsize=10)
        ax2.set_title('L_internode', fontsize=10)
    
    ax2.legend(numpoints=1, bbox_to_anchor=(1.2, 1.2), prop={'size': 9})
    fig.suptitle("Dimension data / Dimension fits")

def dynamique_plot_nff(HS_GL_SSI_data):
    plt.ion()
    
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    for name, color in varieties :
        df_GL = HS_GL_SSI_data[name]['GL']
        df_HS = HS_GL_SSI_data[name]['HS']
        df_SSI = HS_GL_SSI_data[name]['SSI']
        if '11' in df_GL.columns:
            df_GL.plot('TT', '11', style='o'+color, label='GL ' +name+', nff 11')
            df_HS.plot('TT', '11', style='o'+color, label='HS ' +name+', nff 11')
        if '12' in df_GL.columns:
            df_GL.plot('TT', '12', style='^'+color, label='GL ' +name+', nff 12')
            df_HS.plot('TT', '12', style='^'+color, label='HS ' +name+', nff 12')
        if '13' in df_GL.columns:
            df_GL.plot('TT', '13', style='p'+color, label='GL ' +name+', nff 13')
            df_HS.plot('TT', '13', style='p'+color, label='HS ' +name+', nff 13')
  

    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    plt.xlabel("TT")
    
def dynamique_plot(HS_GL_SSI_data, converter = None):
    plt.ion()
    
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    
    for name,color in varieties :
        df_GL = HS_GL_SSI_data[name]['GL']
        df_HS = HS_GL_SSI_data[name]['HS']
        df_SSI = HS_GL_SSI_data[name]['SSI']
        
        if converter is None:
            df_GL.plot('TT', 'mean_pond', style=':o'+color, label='Moyenne pond GL ' +name)
            df_HS.plot('TT', 'mean_pond', style='--o'+color, label='Moyenne pond HS ' +name)
            plt.xlabel("TT")
        else:
            conv = converter[name]
            hsgl = conv(df_GL['TT'])
            hshs = conv(df_HS['TT'])
            plt.plot(hsgl, df_GL['mean_pond'], ':o'+color, label='Moyenne pond GL ' +name)
            plt.plot(hshs, df_HS['mean_pond'], '--o'+color, label='Moyenne pond HS ' +name)
            plt.xlabel("HS")
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
def dynamique_plot_sim(HS_GL_SSI_data, pars, converter = None):
    plt.ion()
    
    #varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b']] #pb Tremie13
    
    for name,color in varieties :
        df_GL = HS_GL_SSI_data[name]['GL']
        df_HS = HS_GL_SSI_data[name]['HS']
        df_SSI = HS_GL_SSI_data[name]['SSI']
        sim = pars[name]
        
        if converter is None:
            df_GL.plot('TT', 'mean_pond', style='o'+color, label='Moyenne pond GL ' +name)
            df_HS.plot('TT', 'mean_pond', style='^'+color, label='Moyenne pond HS ' +name)
            df_SSI.plot('TT', 'mean_pond', style='<'+color, label='Moyenne pond SSI ' +name)
            sim.plot('TT', 'GL', style='-o'+color, label='Sim GL ' +name)
            sim.plot('TT', 'HS', style='-^'+color, label='Sim HS ' +name)
            sim.plot('TT', 'SSI', style='-<'+color, label='Sim SSI ' +name)
            plt.xlabel("TT")
        else:
            conv = converter[name]
            hsgl = conv(df_GL['TT']); hshs = conv(df_HS['TT']); hsssi = conv(df_SSI['TT'])
            hssim = conv(sim['TT'])
            plt.plot(hsgl, df_GL['mean_pond'], 'o'+color, label='Moyenne pond GL ' +name)
            plt.plot(hshs, df_HS['mean_pond'], '^'+color, label='Moyenne pond HS ' +name)
            plt.plot(hsssi, df_SSI['mean_pond'], '<'+color, label='Moyenne pond SSI ' +name)
            plt.plot(hssim, sim['GL'], '-'+color, label='Sim GL ' +name)
            plt.plot(hssim, sim['HS'], '-'+color, label='Sim HS ' +name)
            plt.plot(hssim, sim['SSI'], '-'+color, label='Sim SSI ' +name)
            plt.xlabel("HS")
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
def dynamique_plot_GL_fits(HS_GL_SSI_data, HS_GL_fit, abs='TT', obs=True):
    plt.ion()
    
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    
    for name,color in varieties :
        conv = HS_GL_fit[name]['HS']
        gl_fit = HS_GL_fit[name]['GL']
        hs_obs = HS_GL_SSI_data[name]['HS']
        gl_obs = HS_GL_SSI_data[name]['GL']
        gl_mod = gl_fit.curve()
        #
        hs_obs['HS'] = conv(hs_obs['TT'])
        gl_obs['HS'] = conv(gl_obs['TT'])
        gl_mod['TT'] = conv.TT(gl_mod['HS'])
        #
        hs_obs['TTem'] = conv.TTem(hs_obs['TT'])
        gl_obs['TTem'] = conv.TTem(gl_obs['TT'])
        gl_mod['TTem'] = conv.TTem(gl_mod['TT'])

        #---
        if obs:
            gl_obs.plot(abs, 'mean_pond', style='o'+color, label='Mediane GL ' +name)
            hs_obs.plot(abs, 'mean_pond', style='--o'+color, label='Moyenne pond HS ' +name)
        gl_mod.plot(abs,'GL', style='-'+color, label='GL_fits ' +name)
        
    plt.xlabel(abs)
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})

def density_plot(density_data, fits, HS_converter):
    plt.ion()
    grouped = density_data.groupby('Var')
    
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    
    for name,color in varieties:
        dens = grouped.get_group(name)
        dens['HS'] = HS_converter[name](dens['TT'])
        plt.errorbar(dens['HS'], dens['density'], yerr=dens['SD'], fmt='o'+color, label = name+' density')
   
    for g in fits:       
        if g=='Mercia':
            color='r'
        elif g=='Rht3':
            color='g'
        elif g=='Tremie12':
            color='b'
        else:
            color='m'
        fits[g].plot('HS', 'density', style='-'+color, label=g+' density fits')
      
    plt.title("Plant density"); plt.xlabel("HS"); plt.xlim(xmax=25); plt.ylim(ymin=0); plt.ylim(ymax=350)
    plt.legend(bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    
'''
def plot_tillering(name='Mercia', delta_stop_del=2.5):
    fits = tillering_fits(delta_stop_del)
    fit = fits[name].axis_dynamics(include_MS = False)
    fit.plot('HS', 'total', style='--r', label='Total '+name)
    fit.plot('HS', 'primary', style='-b', label='Primary '+name)
    fit.plot('HS', 'others', style=':g', label= 'Others '+name)

    obs = archidb.tillers_per_plant()
    grouped = obs.groupby('Var')
    obs = grouped.get_group(name)
    obs.plot('HS', ['TP', 'TS', 'TPS','TT3F','FT'],style=['pb','pg','pr','py','pr'])
'''  
  
def multi_plot_tillering(obs_data, fits, HS_converter, delta_stop_del):
    plt.ion()
    
    if not isinstance(delta_stop_del,dict):
        delta_stop_del = {k:delta_stop_del for k in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']}
    
    fig, axes = plt.subplots(nrows=2, ncols=2)
    ax0, ax1, ax2, ax3 = axes.flat
    
    varieties = [['Mercia',ax0,delta_stop_del['Mercia']],['Rht3',ax1,delta_stop_del['Rht3']],['Tremie12',ax2,delta_stop_del['Tremie12']],['Tremie13',ax3,delta_stop_del['Tremie13']]]
    
    for name,ax,delta_stop_del in varieties :

        fit = fits[name].axis_dynamics(include_MS = False)

        ax.plot(fit['HS'], fit['total'], '--r', label='Total')
        ax.plot(fit['HS'], fit['primary'], '-b', label='Primary')
        ax.plot(fit['HS'], fit['others'], ':g', label= 'Others')

        em = fits[name].emited_cohort_density()
        ax.stem(em['delay'],em['total_axis'], markerfmt='xr', linefmt='--r', basefmt='k', label='Total emited cohort density')
        ax.stem(em['delay'],em['other_axis'], markerfmt='xg', linefmt='--g', basefmt='k', label='Others emited cohort density')
        ax.stem(em['delay'],em['primary_axis'], markerfmt='xb', linefmt='--b', basefmt='k', label='Primary emited cohort density')

        grouped = obs_data.groupby('Var'); obs = grouped.get_group(name)
        obs['HS'] = HS_converter[name](obs['TT'])
        ax.plot(obs['HS'], obs['TP'], 'pb', label='TP')
        ax.plot(obs['HS'], obs['TS'], 'pg', label='TS')
        ax.plot(obs['HS'], obs['TPS'], 'pr', label='TPS')
        color='#FFFF00'; ax.plot(obs['HS'], obs['TT3F'], 'p', color=color, label='TT3F')
        ax.plot(obs['HS'], obs['FT'], 'pm', label='FT')
        
        hs_debreg = fits[name].hs_debreg()
        ax.annotate('', xy=(hs_debreg, 0), xytext=(hs_debreg, 1), arrowprops=dict(facecolor='#FE9A2E', shrink=0.00))

        if name is 'Mercia':
            ax.set_xlim([0, 18]); ax.set_ylim([0, 10]) 
        elif name is 'Rht3':
            ax.set_xlim([0, 18]); ax.set_ylim([0, 10]) 
        else :
            ax.set_xlim([0, 18]); ax.set_ylim([0, 7]) 
            
        ax.set_title(name+', delta_stop_del : '+str(delta_stop_del), fontsize=10)
    
    ax1.legend(numpoints=1, bbox_to_anchor=(1.2, 1.2), prop={'size': 9})
    fig.suptitle("Tillering")
   
def graph_primary_emission(archidb):
    plt.ion()
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    
    for name,color in varieties :
        if name is 'Mercia' :
            tdb_name = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Mercia']
        elif name is 'Rht3' :
            tdb_name = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Rht3']
        elif name is 'Tremie12' :
            tdb_name = archidb.Tillering_data_Tremie12_2011_2012()
        else :
            tdb_name = archidb.Tillering_data_Tremie13_2012_2013()

        primary_emission_name = {k[1:]:v for k,v in tdb_name['emission_probabilities'].iteritems() if k != 'TC'}
        df_name = pandas.DataFrame.from_dict(primary_emission_name, orient='index')
        df_name = df_name.reset_index(); df_name.columns = ['talle', 'proba']
            
        df_name = df_name.sort(columns='talle')
        df_name.plot('talle', 'proba', style='--o'+color, label=name)
    
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    plt.xlabel("Talle"); plt.ylabel("%")
    #plt.title("Emissions probabilities")