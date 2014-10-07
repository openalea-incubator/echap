"""
 Ploting function companion module for architectural reconstructions
"""
import matplotlib.pyplot as plt
import pandas

def dynamique_plot_nff(HS_GL_SSI_data, dynT_MS):
    plt.ion()
    
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    for name, color in varieties :
        df_GL = HS_GL_SSI_data[name]['GL']
        df_HS = HS_GL_SSI_data[name]['HS']
        df_SSI = HS_GL_SSI_data[name]['SSI']
        if '11' in df_GL.columns:
            df_GL.plot('TT', '11', style=':o'+color, label='GL ' +name+', nff 11')
            df_HS.plot('TT', '11', style='--o'+color, label='HS ' +name+', nff 11')
        if '12' in df_GL.columns:
            df_GL.plot('TT', '12', style=':^'+color, label='GL ' +name+', nff 12')
            df_HS.plot('TT', '12', style='--^'+color, label='HS ' +name+', nff 12')
        if '13' in df_GL.columns:
            df_GL.plot('TT', '13', style=':p'+color, label='GL ' +name+', nff 13')
            df_HS.plot('TT', '13', style='--p'+color, label='HS ' +name+', nff 13')
  
    vars = [['Others','c'],['Tremie13','y']]
    for name, color in vars :
        a = dynT_MS[name]['a_cohort']*500; b = dynT_MS[name]['a_cohort']*1300
        dynT = pandas.DataFrame({'HS':[500,1300,1900],'dyn':[a, b, b]})
        dynT.plot('HS', 'dyn', style = '-o'+color, linewidth=2, label = 'pente de la droite ' + name)
    
        #T0 = TT_col0 + n0 / a_cohort [TT pr HS=n0]
        #T1 = TT_col0 + dec_internode_elongated_nbr / a_cohort [TT pr HS=n1]
        #c = -((nff - dec_elongated_internode_nbr)-(n2 - n1)) / (TT_col_n_phytomer - T1)
        T0 = dynT_MS[name]['TT_col_0'] + dynT_MS[name]['n0'] / dynT_MS[name]['a_cohort']
        T1 = dynT_MS[name]['TT_col_0'] + (12 - 4.) / dynT_MS[name]['a_cohort']
        C = dynT_MS[name]['TT_col_N_phytomer_potential']
    
        dynT = pandas.DataFrame({'HS':[T0, T1, C],'dyn':[dynT_MS[name]['n0'], dynT_MS[name]['n1'], dynT_MS[name]['n2']]})
        dynT.plot('HS', 'dyn', style = '-^'+color, linewidth=2, label = 'Global '+name)

    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
    plt.xlabel("HS")
    
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
        else:
            conv = converter[name]
            hsgl = conv(df_GL['TT'])
            hshs = conv(df_HS['TT'])
            plt.plot(hsgl, df_GL['mean_pond'], ':o'+color)
            plt.plot(hshs, df_HS['mean_pond'], '--o'+color)
    plt.xlabel("TT")
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})

def density_plot(density_data, fits):
    plt.ion()
    grouped = density_data.groupby('Var')
    
    dens_mercia = grouped.get_group('Mercia'); dens_rht3 = grouped.get_group('Rht3')
    dens_tremie12 = grouped.get_group('Tremie12'); dens_tremie13 = grouped.get_group('Tremie13')

    plt.errorbar(dens_mercia['HS'], dens_mercia['density'], yerr=dens_mercia['SD'], fmt='or', label = 'Mercia density')
    plt.errorbar(dens_rht3['HS'], dens_rht3['density'], yerr=dens_rht3['SD'], fmt='og', label = 'Rht3 density')
    plt.errorbar(dens_tremie12['HS'], dens_tremie12['density'], yerr=dens_tremie12['SD'], fmt='ob', label = 'Tremie12 density')
    plt.errorbar(dens_tremie13['HS'], dens_tremie13['density'], yerr=dens_tremie13['SD'], fmt='om', label = 'Tremie13 density')
   
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
      
    plt.title("Plant density"); plt.xlabel("HS"); plt.ylim(ymin=0); plt.ylim(ymax=350)
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
  
def multi_plot_tillering(obs_data, fits, delta_stop_del):
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