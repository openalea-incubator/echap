""" Ploting function for the echap diseaese progress paper
"""

import matplotlib.pyplot as plt
import pandas
import numpy
from scipy.interpolate import interp1d


import alinea.echap.architectural_data as archidb
from alinea.echap.architectural_reconstructions import EchapReconstructions

fit = EchapReconstructions()

def varieties():
    return ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')
    
def vcolors():
    return {'Mercia':'r', 'Rht3':'g', 'Tremie12':'b', 'Tremie13':'m'}
    
def colors_nff():
    return {10:'m', 11:'g', 12:'r', 13:'c', 14:'y'}
    
def markers_nff():
    return {10: 'd', 11:'o', 12:'^', 13:'s'}

def markers_source():
    return {'tagged':'s', 'sampled':'^', 'symptom':'o', 'estimated_sampled':'*'}


def axe_density_plot(what='primary_axe_density'):

    cols = vcolors()
    fig, ax = plt.subplots()
    obs = archidb.axe_density_aggregated()
    plant_density = fit.density_fits
    axepop = fit.axepop_fits
    hs = fit.HS_fit

    for variety in varieties():
        color = cols[variety]
        df = obs[obs['label'] == variety]
        df = df[df['tag_event'] != 'sowing']
        marker_src = {'plot_data': 's', 'tillering_data': 'd'}
        for src in marker_src:
            if src in df['source'].tolist():
                df_src = df[df['source'] == src]
                ax.errorbar(df_src['TT'].values, df_src[what + '_mean'].values,
                            yerr=df_src[what + '_nanstd'].values, linestyle='',
                            color=color,
                            markerfacecolor=color, markeredgecolor=color,
                            marker=marker_src[src], markersize=7)

        start, end = hs[variety]([0, 2500])
        pd = plant_density[variety]['density_table']
        pd.ix[0, 'HS'] = start
        pd.ix[len(pd) - 1, 'HS'] = end
        pd_fun = interp1d(pd['HS'], pd['density'])
        sim = axepop[variety].axis_dynamics(pd_fun, start=start, end=end)
        sim['TT'] = hs[variety].TT(sim['HS'])
        if what.startswith('primary'):
            ax.plot(sim['TT'], sim['primary'], color=color)
            ax.set_ylim((0, 1500))
            ax.set_ylabel('Axe density (m-2)', fontsize=17)
        else:
            ax.plot(sim['TT'], sim['total'], color=color)
            ax.set_ylim((0, 2000))
            ax.set_ylabel('Axe density (m-2)', fontsize=17)

    ax.set_xlabel('Thermal time (dd)', fontsize=17)
    return ax


def _dim_plot(ax, what, data,model, along='rank', legend=True, xlab='', ylab=''):
    markers = markers_source()
    cols = vcolors()
    vars = varieties()
    
    for variety in vars:
        color = cols[variety]
        fit_dim = model[variety].dimT(nff = None)
        fit_dim.loc[:,'ranktop'] = numpy.max(fit_dim.loc[:,'rank']) - fit_dim.loc[:,'rank'] + 1
        df_dim_var = data[data['label']==variety]
        ax.plot(fit_dim[along], fit_dim.loc[:, what], color = color)
        for src in numpy.unique(df_dim_var['Source']):
            df_src = df_dim_var[df_dim_var['Source']==src]
            ax.errorbar(df_src[along], df_src[what+'_mean'], 
                        yerr=df_src[what+'_std'], linestyle='', color = color,
                        markerfacecolor=color, markeredgecolor=color,
                        marker=markers[src], markersize=7)
    
    ax.set_ylabel(ylab, fontsize = 18)
    ax.set_xlabel(xlab, fontsize = 18)
    
    proxys = []
    labels = []
    for variety, col in cols.iteritems():
        proxys += [plt.Line2D((0,1),(0,0), linestyle='', color = col,
                    markerfacecolor=col, markeredgecolor=col, marker = markers[src])]
        labels += [variety]
    if legend == True:
        ax.legend(proxys, labels, loc=1)
    
    
def dimension_plot(along='ranktop', xlab='Leaf Position',
                   what = ['L_blade', 'H_col']):
    """ Plot comparing blade lengths / collar heights for the differents varieties
        phytomers numbered from top
    """
    
    data = archidb.Dim_data()
    data = data[pandas.notnull(data['nff'])]
    data.loc[:,'ranktop'] = data.loc[:,'nff'] - data.loc[:,'rank'] + 1
    df_dim,df_dim_nff = archidb.dimensions_aggregated(data, along)
    model = fit.dimension_fits
    
    fig, axs = plt.subplots(1, 2)

    labels = {'L_blade':'Leaf length (cm)', 'H_col':'Leaf height (cm)',
              'W_blade': 'Leaf width (cm)', 'L_sheath': 'Sheath length (cm)'}
    for i, ax in enumerate(axs.flat):
        _dim_plot(ax, what[i], df_dim, model,along, legend=False, xlab=xlab, ylab=labels[what[i]])
        
def pheno_plot():
    """ Compare emergence / senescence for the different varieties
    """
    model = fit.GL_fits
    
    markers = markers_nff()
    vcols = vcolors()
    vars = varieties()
    markersize = 7
    
    fig, axs = plt.subplots(1, 2)
    
    # HS plot
    ax = axs[0]
    obs_nff, obsM = archidb.haun_stage_aggregated()
    
    for variety in vars:
    
        dfobsM = obsM[obsM['label'] == variety]
        dfobs_nff = obs_nff[obs_nff['label'] == variety]
    
        hsfit = model[variety].hsfit
        ranks = numpy.arange(0, hsfit.HSflag(),0.1)
        
        nffs = fit.axepop_fits[variety].sorted_nff()[:2]
        ranks_nff = map(lambda x: numpy.arange(0,x+0.1,0.1), nffs)
        top_ranks_nff = map(lambda x: max(x) - x + 1, ranks_nff)
        for i in range(2):
            ranks_nff[i] = ranks_nff[i][numpy.where(top_ranks_nff[i] <= min(nffs))]
            top_ranks_nff[i] = top_ranks_nff[i][numpy.where(top_ranks_nff[i] <= min(nffs))]
        TThs_nff = map(lambda x: hsfit.TTem(hsfit.TT(ranks_nff[x], nff=nffs[x])), range(2))
        
        color = vcols[variety]
        ax.plot(max(ranks) - ranks + 1, hsfit.TTem(hsfit.TT(ranks)), color=color)
        ax.fill_between(top_ranks_nff[0], TThs_nff[0], TThs_nff[1], color=color, alpha=0.2)
        
        #ax.errorbar(max(ranks) - dfobsM['HS_mean'] + 1, hsfit.TTem(dfobsM['TT']), xerr=dfobsM['HS_std'], color = color, 
        #            linestyle='', markerfacecolor='None', markeredgecolor=color,
        #            marker='o', markersize=markersize)
        for i in range(2):
            nff = nffs[i]
            df_nff = dfobs_nff[dfobs_nff['nff']==nff].reset_index()
            ax.errorbar(nff - df_nff['HS_mean'] + 1, hsfit.TTem(df_nff['TT']) , xerr=df_nff['HS_std'], color = color, 
                        linestyle='', markerfacecolor=color, markeredgecolor=color,
                        marker=markers[nff], markersize=markersize)
        
        ax.set_xlabel('Leaf position', fontsize = 18)
        ax.set_ylabel('Thermal time at full expansion (dd)', fontsize = 18)
    
    
    # Green leaf duration plot
    ax = axs[1]
    obs_nff, obsM = archidb.ssi_aggregated()
    
    for variety in vars:
    
        dfobsM = obsM[obsM['label'] == variety]
        dfobs_nff = obs_nff[obs_nff['label'] == variety]
    
        hsfit = model[variety].hsfit
        gl = model[variety]
        ranks = numpy.arange(0, hsfit.HSflag(),0.1)
    
        TTsen = gl.TTsen()(ranks)
        TTem = hsfit.TTemleaf(ranks)
        
        nffs = fit.axepop_fits[variety].sorted_nff()[:2]
        ranks_nff = map(lambda x: numpy.arange(0,x+0.1,0.1), nffs)
        top_ranks_nff = map(lambda x: max(x) - x + 1, ranks_nff)
        for i in range(2):
            ranks_nff[i] = ranks_nff[i][numpy.where(top_ranks_nff[i] <= min(nffs))]
            top_ranks_nff[i] = top_ranks_nff[i][numpy.where(top_ranks_nff[i] <= min(nffs))]
        
        TTsen_nff = map(lambda x: gl.TTsen(nff=nffs[x])(ranks_nff[x]), range(2))
        TTem_nff = map(lambda x: hsfit.TTemleaf(ranks_nff[x], nff=nffs[x]), range(2))
        
        color = vcols[variety]
        ax.plot(max(ranks) - ranks + 1, TTsen - TTem, color = color)    
        ax.fill_between(top_ranks_nff[0],TTsen_nff[0] - TTem_nff[0], TTsen_nff[1] - TTem_nff[1], color=color, alpha=0.2)
        
        #ax.errorbar(max(ranks) - dfobsM['SSI_mean'] + 1, dfobsM['TT'] - hsfit.TTemleaf(dfobsM['SSI_mean']), xerr=dfobsM['SSI_std'], color = color, 
        #            linestyle='', markerfacecolor='None', markeredgecolor=color,
        #            marker='o', markersize=markersize)
        for i in range(2):
            nff = nffs[i]
            df_nff = dfobs_nff[dfobs_nff['nff']==nff].reset_index()
            ax.errorbar(nff - df_nff['SSI_mean'] + 1, df_nff['TT'] - hsfit.TTemleaf(df_nff['SSI_mean'], nff), xerr=df_nff['SSI_std'], color = color, 
                        linestyle='', markerfacecolor=color, markeredgecolor=color,
                        marker=markers[nff], markersize=markersize)
                        
        ax.set_xlabel('Leaf position', fontsize = 18)
        ax.set_ylabel('Green life span (dd)', fontsize = 18)


