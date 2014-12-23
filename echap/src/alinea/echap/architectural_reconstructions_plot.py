"""
 Ploting function companion module for architectural reconstructions
"""
import matplotlib.pyplot as plt
import pandas
import numpy as np

def dimension_plot_other(dimension_data, fits):
    plt.ion()
    
    fig, axes = plt.subplots(nrows=1, ncols=3)
    ax0, ax1, ax2 = axes.flat
    
    varieties = [['Mercia','r',12],['Rht3','g',11],['Tremie12','b',13],['Tremie13','m',11]]
    for name, color, nff in varieties :
        #premier nff, ligne du haut
        d_data = dimension_data[name, nff]
        dim_data = d_data.where((pandas.notnull(d_data)), None)
        ax0.plot(dim_data['index_phytomer'], dim_data['L_sheath'], 'o'+color, label = 'Obs '+name)
        ax1.plot(dim_data['index_phytomer'], dim_data['L_internode'], 'o'+color, label = 'Obs '+name)
        ax2.plot(dim_data['index_phytomer'], dim_data['W_internode'], 'o'+color, label = 'Obs '+name)
        
        dim_fit = fits[name][nff]
        ax0.plot(dim_fit['index_phytomer'], dim_fit['L_sheath'], '-'+color, label = '_nolegend_')
        ax1.plot(dim_fit['index_phytomer'], dim_fit['L_internode'], '-'+color, label = '_nolegend_')
        ax2.plot(dim_fit['index_phytomer'], dim_fit['W_internode'], '-'+color, label = 'Fit '+name)
        
        #legende
        ax0.set_title('L_sheath nff maj', fontsize=9)
        ax1.set_title('L_internode nff maj', fontsize=9)
        ax2.set_title('W_internode nff maj', fontsize=9)
        ax2.legend(numpoints=1, bbox_to_anchor=(1.3, 1.1), prop={'size': 9})

def dimension_plot(dimension_data, fits, leaf_fits, scan, scan_old):
    plt.ion()
    
    fig, axes = plt.subplots(nrows=2, ncols=3)
    ax0, ax1, ax2, ax3, ax4, ax5 = axes.flat
    
    varieties = [['Mercia','r',[11,12]],['Rht3','g',[11,12]],['Tremie12','b',[12,13]],['Tremie13','m',[11,12]]]
    for name, color, nff_lst in varieties :
        #premier nff, ligne du haut
        d_data = dimension_data[name, nff_lst[0]]
        dim_data = d_data.where((pandas.notnull(d_data)), None)
        ax0.plot(dim_data['index_phytomer'], dim_data['L_blade'], 'o'+color, label = 'Obs '+name)
        if name=='Mercia' or name=='Rht3' or name=='Tremie13':
            ax1.plot(dim_data['index_phytomer'], dim_data['W_blade'], 'o'+color, label = '_nolegend_')
        else : #pas de donnees largeur par nff pour tremie12 donc on met des triangles
            ax1.plot(dim_data['index_phytomer'], dim_data['W_blade'], '^'+color, label = '_nolegend_')
        dim_fit = fits[name][nff_lst[0]]
        ax0.plot(dim_fit['index_phytomer'], dim_fit['L_blade'], '-'+color, label = '_nolegend_')
        ax1.plot(dim_fit['index_phytomer'], dim_fit['W_blade'], '-'+color, label = '_nolegend_')
        #FF
        form_factor_bas = leaf_fits[name].form_factor()[1]
        form_factor_haut = leaf_fits[name].form_factor()[2]
        dim_fit['L*W_blade*FF'] = dim_fit['L_blade']
        if name=='Tremie12': #FF haut pour les 4 feuilles du haut, FF bas pour les autres feuilles
            dim_fit.ix[dim_fit.index_phytomer.isin([1,2,3,4,5,6,7,8]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_bas
            dim_fit.ix[dim_fit.index_phytomer.isin([9,10,11,12]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_haut   
        else:
            dim_fit.ix[dim_fit.index_phytomer.isin([1,2,3,4,5,6,7]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_bas
            dim_fit.ix[dim_fit.index_phytomer.isin([8,9,10,11]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_haut 
        ax2.plot(dim_fit['index_phytomer'], dim_fit['L*W_blade*FF'], '-x'+color, label = '_nolegend_')
        #second nff, ligne du bas
        d_data = dimension_data[name, nff_lst[1]]
        dim_data = d_data.where((pandas.notnull(d_data)), None)
        ax3.plot(dim_data['index_phytomer'], dim_data['L_blade'], 'o'+color, label = '_nolegend_')
        if name=='Mercia' or name=='Rht3' or name=='Tremie13':
            ax4.plot(dim_data['index_phytomer'], dim_data['W_blade'], 'o'+color, label = '_nolegend_')
        else: #pas de donnees largeur par nff pour tremie12 donc on met des triangles
            ax4.plot(dim_data['index_phytomer'], dim_data['W_blade'], '^'+color, label = '_nolegend_')
        dim_fit = fits[name][nff_lst[1]]
        ax3.plot(dim_fit['index_phytomer'], dim_fit['L_blade'], '-'+color, label = '_nolegend_')
        ax4.plot(dim_fit['index_phytomer'], dim_fit['W_blade'], '-'+color, label = '_nolegend_')
        #FF
        dim_fit['L*W_blade*FF'] = dim_fit['L_blade']
        if name=='Tremie12':
            dim_fit.ix[dim_fit.index_phytomer.isin([1,2,3,4,5,6,7,8,9]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_bas
            dim_fit.ix[dim_fit.index_phytomer.isin([10,11,12,13]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_haut   
        else:
            dim_fit.ix[dim_fit.index_phytomer.isin([1,2,3,4,5,6,7,8]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_bas
            dim_fit.ix[dim_fit.index_phytomer.isin([9,10,11,12]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_haut
        ax5.plot(dim_fit['index_phytomer'], dim_fit['L*W_blade*FF'], '-'+color, label = 'Wfit*Lfit*FF '+name)
        
        #scan obs - pas de donnees scan par nff pour tremie12 donc on met des triangles
        if name=='Tremie12' or name=='Tremie13':
            scan_var = scan[scan['var']==name]
            scaned = scan_var.groupby('moyenne id_Feuille')
            res = scaned['Area A_bl'].agg([np.mean, np.std])
            res = res.reset_index()
            ax2.errorbar(res['moyenne id_Feuille'], res['mean'], yerr=res['std'], fmt='o'+color, label='scan + sd '+name)
            ax5.errorbar(res['moyenne id_Feuille'], res['mean'], yerr=res['std'], fmt='o'+color)
            
        #scan_old (scan Mercia / Rht3 de 2009) -> pas ici car confusion (on ne cherche pas a passer dans ces points)
        '''if name=='Mercia' or name=='Rht3':
            scan_var = scan_old[scan_old['variety']==name]
            scaned = scan_var.groupby('nmax')
            scaned = scaned.reset_index()
            for nmaxindex in [11,12]:
                res = scaned[scaned['nmax']==nmaxindex]
                if nmaxindex==11:
                    ax2.errorbar(res['rank'], res['A_bl'], yerr=res['sd_A_bl'], fmt='o'+color, label='scan 2009 + sd '+name)
                else :
                    ax5.errorbar(res['rank'], res['A_bl'], yerr=res['sd_A_bl'], fmt='o'+color, label='scan 2009 + sd '+name)'''
            
        
        ax0.set_ylim(0, 30); ax3.set_ylim(0, 30)
        ax1.set_ylim(0, 2.5); ax4.set_ylim(0, 2.5)
        ax2.set_ylim(0, 40); ax5.set_ylim(0, 40)
        ax0.set_xlim(0, 14); ax1.set_xlim(0, 14); ax2.set_xlim(0, 14)
        ax3.set_xlim(0, 14); ax4.set_xlim(0, 14); ax5.set_xlim(0, 14)
        
        ax3.set_ylabel('longueur'); ax4.set_ylabel('largeur'); ax5.set_ylabel('surface')
        ax3.set_xlabel('n_phytomer'); ax4.set_xlabel('n_phytomer'); ax5.set_xlabel('n_phytomer')
        #legende
        ax0.set_title('L_blade nff min', fontsize=9)#; ax3.set_title('L_blade nff max', fontsize=9)
        ax1.set_title('W_blade nff min', fontsize=9)#; ax4.set_title('W_blade nff max', fontsize=9)
        ax2.set_title('L_blade*W_blade nff min + scan', fontsize=9)#; ax5.set_title('L_blade*W_blade nff max + scan', fontsize=9)
        ax0.legend(numpoints=1, bbox_to_anchor=(1.1, 1.2), prop={'size': 9})
        ax2.legend(numpoints=1, bbox_to_anchor=(1.3, 1.1), prop={'size': 9})
    '''
    #old = tracer seulement la modalite majoritaire haut=data, bas=fit
    varieties = [['Mercia','y',12],['Rht3','g',11],['Tremie12','b',13],['Tremie13','m',12]]
    for name, color, nff in varieties :   
        d_data = dimension_data[name,nff]
        dim_data = d_data.where((pandas.notnull(d_data)), None)
        ax0.plot(dim_data['index_phytomer'], dim_data['L_blade'], '-^'+color, label = '_nolegend_')
        ax1.plot(dim_data['index_phytomer'], dim_data['W_blade'], '-o'+color, label = '_nolegend_')
        dim_data['L*W_blade'] = dim_data['L_blade'] * dim_data['W_blade']
        ax2.plot(dim_data['index_phytomer'], dim_data['L*W_blade'], '-x'+color, label=name)

        dim_fit = fits[name][nff]
        ax3.plot(dim_fit['index_phytomer'], dim_fit['L_blade'], '-^'+color, label = '_nolegend_')
        ax4.plot(dim_fit['index_phytomer'], dim_fit['W_blade'], '-o'+color, label = '_nolegend_')
        form_factor_bas = leaf_fits[name].form_factor()[1]
        form_factor_haut = leaf_fits[name].form_factor()[2]
        
        dim_fit['L*W_blade*FF'] = dim_fit['L_blade']
        if name=='Tremie12' or name=='Tremie13':
            dim_fit.ix[dim_fit.index_phytomer.isin([1,2,3,4,5,6,7,8,9]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_bas
            dim_fit.ix[dim_fit.index_phytomer.isin([10,11,12,13]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_haut   
        else:
            dim_fit.ix[dim_fit.index_phytomer.isin([1,2,3,4,5,6,7,8]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_bas
            dim_fit.ix[dim_fit.index_phytomer.isin([9,10,11,12]), 'L*W_blade*FF'] = dim_fit['L_blade'] * dim_fit['W_blade'] * form_factor_haut
            
        ax5.plot(dim_fit['index_phytomer'], dim_fit['L*W_blade*FF'], '-x'+color, label = '_nolegend_')
        
        #scan obs
        if name=='Tremie12' or name=='Tremie13':
            scan_var = scan[scan['var']==name]
            if name=='Tremie12':
                color2='c'
            else:
                color2='r'
            ax5.plot(scan_var['moyenne id_Feuille'], scan_var['Area A_bl'], 'o'+color2, label='scan '+name)

        ax0.set_title('L_blade', fontsize=10); ax3.set_title('L_blade fit', fontsize=10)
        ax1.set_title('W_blade', fontsize=10); ax4.set_title('W_blade fit', fontsize=10)
        ax2.set_title('L_blade * W_blade', fontsize=10)
        ax5.set_title('L_blade * W_blade * FF + scan', fontsize=10)
        ax0.set_ylim(0, 30); ax3.set_ylim(0, 30)
        ax1.set_ylim(0, 2.5); ax4.set_ylim(0, 2.5)
        ax2.set_ylim(0, 40); ax5.set_ylim(0, 40)
        ax0.set_xlim(0, 14); ax1.set_xlim(0, 14); ax2.set_xlim(0, 14)
        ax3.set_xlim(0, 14); ax4.set_xlim(0, 14); ax5.set_xlim(0, 14)
        
    ax2.legend(numpoints=1, bbox_to_anchor=(1.2, 1.2), prop={'size': 9})
    ax5.legend(numpoints=1, bbox_to_anchor=(1.2, 1.), prop={'size': 9})
    fig.suptitle("Dimension data / Dimension fits")'''

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
    '''   
    #varieties = [['Tremie12','b','o'],['Tremie13','m','^']]
    varieties = [['Tremie12','b','o']]
    for name, color, picto in varieties :
        df_GL = HS_GL_SSI_data[name]['GL']
        df_HS = HS_GL_SSI_data[name]['HS']
        df_SSI = HS_GL_SSI_data[name]['SSI']
        if '11' in df_GL.columns:
            #df_GL.plot('TT', '11', style='o'+color, label='GL ' +name+', nff 11')
            #df_HS.plot('TT', '11', style='o'+color, label='HS ' +name+', nff 11')
            df_SSI.plot('TT', '11', style=picto+'-b', label='SSI ' +name+' nff11')
        if '12' in df_GL.columns:
            #df_GL.plot('TT', '12', style='^'+color, label='GL ' +name+', nff 12')
            #df_HS.plot('TT', '12', style='^'+color, label='HS ' +name+', nff 12')
            df_SSI.plot('TT', '12', style=picto+'-g', label='SSI ' +name+' nff12')
        if '13' in df_GL.columns:
            #df_GL.plot('TT', '13', style='p'+color, label='GL ' +name+', nff 13')
            #df_HS.plot('TT', '13', style='p'+color, label='HS ' +name+', nff 13')
            df_SSI.plot('TT', '13', style=picto+'-r', label='SSI ' +name+' nff13')
    #donnees guillaume
    TT_Tremie12=[1167,1237,1351,1435,1512,1608,1736,1848,2029]
    #SSI_Tremie12_12=[6.725,7.0125,7.105,7.40833333333,7.76666666667,8.64666666667,9.31166666667,10.3366666667,11.2683333333]
    #SSI_Tremie12_13=[6.78235294118,6.9350678733,6.98733893557,7.00716666667,7.98235294118,8.99485294118,9.76882352941,10.8918382353,11.9805882353]
    #SSI_Tremie12_14=[6.86142857143,7.05,7.36833333333,7.70166666667,8.53857142857,9.21,10.0285714286,11.62,12.7683333333]
    SSI_Tremie12_12 = [6.68333333333,6.67583333333,6.70166666667,7.11361111111,7.38888888889,7.70333333333,8.84666666667,9.57333333333,10.7133333333]
    SSI_Tremie12_13 = [6.86266381766,6.65133903134,6.68070730831,7.52593162393,7.84303215303,8.69962962963,9.29888888889,10.1140740741,11.2477777778]
    
    #Tremie12_g = zip(TT_Tremie12, SSI_Tremie12_12, SSI_Tremie12_13, SSI_Tremie12_14)
    #Tremie12_g = pandas.DataFrame(data = Tremie12_g, columns=['TT','12','13','14'])
    Tremie12_g = zip(TT_Tremie12, SSI_Tremie12_12, SSI_Tremie12_13)
    Tremie12_g = pandas.DataFrame(data = Tremie12_g, columns=['TT','12','13'])
    Tremie12_g.plot('TT','12', style='--^g', label="SSI Tremie12 nff12 Guillaume")
    Tremie12_g.plot('TT','13', style='--^r', label="SSI Tremie12 nff13 Guillaume")
    #Tremie12_g.plot('TT','14', style='--^y', label="SSI Tremie12 nff14 Guillaume")
    
    TT_Tremie13=[1281,1325,1411,1542,1635,1759,1880,1995,2135]
    #SSI_Tremie13_11=[7.39739130435,7.51304347826,7.83608695652,8.39913043478,9.32782608696,10.2504347826,11.0,11.0,11.0]
    #SSI_Tremie13_12=[7.79871052632,7.87494987469,8.1102443609,8.98523809524,9.952125,10.9421052632,11.8189655172,11.9388888889,12.0]
    SSI_Tremie13_11 = [7.3646875,7.475625,7.59,8.0071875,8.419375,9.5040625,10.7171875,10.9958064516,11.0]
    SSI_Tremie13_12 = [7.75923076923,7.85022435897,8.07257326007,8.54,9.42313186813,10.35,11.3908333333,11.8983333333,12.0]
    
    Tremie13_g = zip(TT_Tremie13, SSI_Tremie13_11, SSI_Tremie13_12)
    Tremie13_g = pandas.DataFrame(data = Tremie13_g, columns=['TT','11','12'])
    Tremie13_g.plot('TT','11', style='--^b', label="SSI Tremie13 nff11 Guillaume")
    Tremie13_g.plot('TT','12', style='--^g', label="SSI Tremie13 nff12 Guillaume")
    '''
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
    
    varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b'],['Tremie13','m']]
    #varieties = [['Mercia','r'],['Rht3','g'],['Tremie12','b']] #pb Tremie13
    
    for name,color in varieties :
        df_GL = HS_GL_SSI_data[name]['GL']
        df_HS = HS_GL_SSI_data[name]['HS']
        df_SSI = HS_GL_SSI_data[name]['SSI']
        sim = pars[name]
        
        if converter is None:
            df_GL.plot('TT', 'mean_pond', style='o'+color, label='Moyenne pond GL ' +name)
            df_HS.plot('TT', 'mean_pond', style='^'+color, label='Moyenne pond HS ' +name)
            df_SSI.plot('TT', 'mean_pond', style='<'+color, label='Moyenne pond SSI ' +name)
            sim.plot('TT', 'GL', style='--'+color, label='Sim GL ' +name)
            sim.plot('TT', 'HS', style='-'+color, label='Sim HS ' +name)
            sim.plot('TT', 'SSI', style=':'+color, label='Sim SSI ' +name)
            plt.xlabel("TT")
            plt.ylim(ymin=0); plt.ylim(ymax=14)
        else:
            conv = converter[name]
            hsgl = conv(df_GL['TT']); hshs = conv(df_HS['TT']); hsssi = conv(df_SSI['TT'])
            hssim = conv(sim['TT'])
            plt.plot(hsgl, df_GL['mean_pond'], 'o'+color, label='Moyenne pond GL ' +name)
            plt.plot(hshs, df_HS['mean_pond'], '^'+color, label='Moyenne pond HS ' +name)
            plt.plot(hsssi, df_SSI['mean_pond'], '<'+color, label='Moyenne pond SSI ' +name)
            plt.plot(hssim, sim['GL'], '--'+color, label='Sim GL ' +name)
            plt.plot(hssim, sim['HS'], '-'+color, label='Sim HS ' +name)
            plt.plot(hssim, sim['SSI'], ':'+color, label='Sim SSI ' +name)
            plt.xlabel("HS")
            plt.ylim(ymin=0); plt.ylim(ymax=14)
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
        ax.plot(fit['HS'], fit['3F'], '--y', label= 'TT3F')

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
            ax.set_xlim([0, 18]); ax.set_ylim([0, 10]) 
            
        ax.set_title(name+', delta_stop_del : '+str(delta_stop_del), fontsize=10)
    
    ax1.legend(numpoints=1, bbox_to_anchor=(1.2, 1.2), prop={'size': 9})
    fig.suptitle("Tillering")
    
def tillering_primary(obs_data, fits, HS_converter, delta_stop_del):
    plt.ion()
    
    if not isinstance(delta_stop_del,dict):
        delta_stop_del = {k:delta_stop_del for k in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']}

    varieties = [['Mercia',delta_stop_del['Mercia']],['Rht3',delta_stop_del['Rht3']],['Tremie12',delta_stop_del['Tremie12']],['Tremie13',delta_stop_del['Tremie13']]]
    
    for name,delta_stop_del in varieties :
    
        if name=='Mercia':
            color='r'
        elif name=='Rht3':
            color='g'
        elif name=='Tremie12':
            color='b'
        elif name=='Tremie13':
            color='m'

        fit = fits[name].axis_dynamics(include_MS = False)

        if name=='Mercia':
            plt.plot(fit['HS'], fit['primary'], '-'+color, label='Primary')
        else :
            plt.plot(fit['HS'], fit['primary'], '-'+color, label='_nolegend_')           

        grouped = obs_data.groupby('Var'); obs = grouped.get_group(name)
        obs['HS'] = HS_converter[name](obs['TT'])
        
        if name=='Mercia':
            plt.plot(obs['HS'], obs['TP'], 'o'+color, label='TP', markersize=9)
            plt.plot(obs['HS'], obs['FT'], '^'+color, label='FT', markersize=9)
        else :
            plt.plot(obs['HS'], obs['TP'], 'o'+color, label='_nolegend_', markersize=9)
            plt.plot(obs['HS'], obs['FT'], '^'+color, label='_nolegend_', markersize=9)
        
    plt.xlim(xmax=20)
    plt.xlabel('HS')
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
    
def tillering_tot(obs_data, fits, HS_converter, delta_stop_del):
    plt.ion()
    
    if not isinstance(delta_stop_del,dict):
        delta_stop_del = {k:delta_stop_del for k in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']}

    varieties = [['Mercia',delta_stop_del['Mercia']],['Rht3',delta_stop_del['Rht3']],['Tremie12',delta_stop_del['Tremie12']],['Tremie13',delta_stop_del['Tremie13']]]
    
    for name,delta_stop_del in varieties :
    
        if name=='Mercia':
            color='r'
        elif name=='Rht3':
            color='g'
        elif name=='Tremie12':
            color='b'
        elif name=='Tremie13':
            color='m'

        fit = fits[name].axis_dynamics(include_MS = False)

        if name=='Mercia':
            plt.plot(fit['HS'], fit['total'], '-'+color, label='Total')
        else :
            plt.plot(fit['HS'], fit['total'], '-'+color, label='_nolegend_')        

        grouped = obs_data.groupby('Var'); obs = grouped.get_group(name)
        obs['HS'] = HS_converter[name](obs['TT'])
        
        if name=='Mercia':
            plt.plot(obs['HS'], obs['TPS'], 'p'+color, label='TPS', markersize=9)
            plt.plot(obs['HS'], obs['FT'], '^'+color, label='FT', markersize=9)
        else :
            plt.plot(obs['HS'], obs['TPS'], 'p'+color, label='_nolegend_', markersize=9)
            plt.plot(obs['HS'], obs['FT'], '^'+color, label='_nolegend_', markersize=9)
        
    plt.xlim(xmax=20)
    plt.xlabel('HS')
    plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})
   
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