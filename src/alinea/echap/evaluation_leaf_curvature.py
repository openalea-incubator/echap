""" Evaluation of leaf curvature in ADEL"""
from alinea.echap.architectural_data import xydb_reader
from alinea.adel.leaf.curvature import *
from alinea.echap.cache_simulation import build_canopies, get_midribs
import numpy
import pandas
import matplotlib.pyplot as plt
from math import ceil
import collections
from openalea.deploy.shared_data import shared_data
import alinea.echap

from alinea.echap.architectural_reconstructions import echap_reconstructions
from functools import reduce

def get_file_name(variety = 'Tremie12', nplants = 30,):
    return 'curvature_'+variety.lower() + '_' + str(nplants) + 'pl.csv'

def get_file_path(variety = 'Tremie12', nplants = 30):
    filename = get_file_name(variety = variety, nplants = nplants)
    return str(shared_data(alinea.echap)/'curvature_simulations'/filename)                                                           
                                                        
# Run Simulations
#
def simulation_curvatures(variety='Tremie12', nplants=30, tag='reference',
                          rep=1):
    harvests = {'Tremie12': {1: '2012-03-09', 2: '2012-05-09', 3: '2012-04-11',
                             4: '2012-06-12'}, 'Tremie13': {1: 11.1},
                'Mercia': {1: 4, 2: 7, 3: 8, 4: 13},
                'Rht3': {1: 4, 2: 6, 3: 8, 5: 15.5}}
    harvests = harvests[variety]
    build_canopies(variety=variety, nplants=nplants, tag=tag, rep=rep,
                   at=list(harvests.values()))
    sims = {
    h: get_midribs(variety=variety, nplants=nplants, daydate=harvests[h],
                   tag=tag, rep=rep) for h in harvests}
    # sims = {h:adel.get_midribs(adel.setup_canopy(age=HSconv.TT(harvests[h])),resample=True) for h in harvests}
    ldf = []
    for h,df_sim in sims.items():
        df = df_sim[df_sim['axe'] == 'MS']
        df.loc[:,'variety_name'] = variety
        df.loc[:,'side'] = 1
        df.loc[:,'harvest'] = h
        df.rename(columns={'metamer':'rank', 'ntop':'ranktop','vid':'inerv'},inplace=True)
        rmax = df.set_index('plant').groupby(level=0).max()['rank']
        df.loc[:,'relative_ranktop'] = numpy.array([rmax[p] for p in df['plant']]) - df['rank'] + 1
        ldf.append(df)
    df = pandas.concat(ldf)
    df = df.reset_index(drop=True)
    return df
    
def save_simulation_curvatures(variety = 'Tremie12', nplants = 30, 
                                reset = False, reset_data = False):
    df = simulation_curvatures(variety=variety, nplants=nplants)
    df.to_csv(get_file_path(variety = variety, nplants = nplants), index = False)
  
def load_simulation_curvature(variety = 'Tremie12', nplants = 30):
    return pandas.read_csv(get_file_path(variety = variety, nplants = nplants))  
    
# Read curvature data ##############################################################################
def is_iterable(obj):
    return isinstance(obj, collections.Iterable)

def plot_leaf_curvature(name = 'Tremie12', nplants = 30, numbering = 'relative_ranktop', at=None,
                        alternate = True, hins_mean = False, 
                        fixed_color = None, alpha = 1., axs = None, 
                        fixed_xlims = [-15, 15], fixed_ylims = [-10, 90], set_lims = True, add_sim=False, annotate=True):

    df = xydb_reader(name = name)

    if at is None:
        if axs is None:
            fig, axs = plt.subplots(1, len(numpy.unique(df['harvest'])))
        if is_iterable(axs):
            iter_ax = axs.flat
        else:
            iter_ax = iter([axs])
    else:
        fig, axs = plt.subplots()
        iter_ax = iter([axs])
        what = at

    for i, h in enumerate(numpy.unique(df['harvest'])):
        if at is None:
            what = i
        if i == what:
            ax = next(iter_ax)
            ax.set_color_cycle(None)
            colors = ax._get_lines.color_cycle
            df_h = df[df['harvest'] == h]
            if hins_mean == True:
                df_grouped = df_h.groupby(numbering).mean()
            for leaf in numpy.unique(df_h[numbering]):
                if fixed_color is not None:
                    color = fixed_color
                else:
                    color = next(colors)
                for pl in numpy.unique(df_h['plant']):
                    df_pl = df_h[(df_h['plant'] == pl) & (df_h[numbering] == leaf)]
                    if len(df_pl)>0:
                        if alternate == True:
                            side = df_pl.loc[:,'side'] * numpy.power(-1, leaf%2)
                        else:
                            side = 1

                        if hins_mean == True:
                            ax.plot(df_pl.loc[:,'x']*side, df_pl.loc[:,'y']+df_grouped.loc[leaf,'hins'],
                                    color = color, alpha = alpha)
                        else:
                            ax.plot(df_pl.loc[:,'x']*side, df_pl.loc[:,'y']+df_pl.loc[:,'hins'],
                                    color = color, alpha = alpha)
            if add_sim:
                dfsim = load_simulation_curvature(name, nplants=nplants)
                ax.set_color_cycle(None)
                colors = ax._get_lines.color_cycle
                df_h = dfsim[dfsim['harvest'] == h]
                if hins_mean == True:
                    df_grouped = df_h.groupby(numbering).mean()
                for leaf in numpy.unique(df_h[numbering]):
                    if fixed_color is not None:
                        color = fixed_color
                    else:
                        color = next(colors)
                    for pl in numpy.unique(df_h['plant']):
                        df_pl = df_h[(df_h['plant'] == pl) & (df_h[numbering] == leaf)]
                        if len(df_pl)>0:
                            if alternate == True:
                                side = df_pl.loc[:,'side'] * numpy.power(-1, leaf%2)
                            else:
                                side = 1

                            if hins_mean == True:
                                ax.plot(df_pl.loc[:,'x']*side, df_pl.loc[:,'y']+df_grouped.loc[leaf,'hins'],
                                        color = color, alpha = alpha,linestyle = ':')
                            else:
                                ax.plot(df_pl.loc[:,'x']*side, df_pl.loc[:,'y']+df_pl.loc[:,'hins'],
                                        color = color, alpha = alpha,linestyle = ':')
            if annotate:
                ax.annotate('Harvest %d' % h, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    
    if is_iterable(axs):
        if set_lims == True:
            xlims = [min([ax.get_xlim()[0] for ax in axs]), max([ax.get_xlim()[1] for ax in axs])]
            ylims = [min([ax.get_ylim()[0] for ax in axs]), max([ax.get_ylim()[1] for ax in axs])]
            for ax in axs:
                ax.set_xlim(xlims)
                ax.set_ylim(ylims)
        else:
            for ax in axs:
                ax.set_xlim(fixed_xlims)
                ax.set_ylim(fixed_ylims)
        
def plot_mean_by_point_leaf_curvature(name = 'Tremie12', nplants = 30,
                                      numbering = 'relative_ranktop', 
                                      hins_mean = True, fixed_xlims = [-15, 15],
                                      fixed_ylims = [-10, 90], set_lims = False):
    df = xydb_reader(name = name)
    df.loc[:, 'i_point'] = [float(i)%20 for i in df.index]
    fig, axs = plt.subplots(1, len(numpy.unique(df['harvest'])))
    
    plot_leaf_curvature(name = name, nplants = nplants, numbering = numbering, 
                        alternate = True, hins_mean = hins_mean, alpha = 0.2, axs = axs,
                        fixed_xlims = fixed_xlims,  fixed_ylims = fixed_ylims, set_lims = False)
    if is_iterable(axs):
        iter_ax = axs.flat
    else:
        iter_ax = iter([axs])
    
    for h in numpy.unique(df['harvest']):
        ax = next(iter_ax)
        ax.set_color_cycle(None)
        colors = ax._get_lines.color_cycle
        df_h = df[df['harvest'] == h]
        df_grouped = df_h.groupby(numbering).mean()
        for leaf in numpy.unique(df_h[numbering]):
            color = next(colors)
            df_lf = df_h[(df_h[numbering] == leaf)]
            side = df_lf.loc[:,'side'] * numpy.power(-1, leaf%2)
            df_lf.loc[:, 'x'] *= side
            df_mean = df_lf.groupby('i_point').mean()
            ax.plot(df_mean.loc[:,'x'], df_mean.loc[:,'y']+df_grouped.loc[leaf,'hins'],
                        color = color, linewidth = 2)
        ax.annotate('Harvest %d' % h, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)

    if is_iterable(axs):
        if set_lims == True:
            xlims = [min([ax.get_xlim()[0] for ax in axs]), max([ax.get_xlim()[1] for ax in axs])]
            ylims = [min([ax.get_ylim()[0] for ax in axs]), max([ax.get_ylim()[1] for ax in axs])]
            for ax in axs:
                ax.set_xlim(xlims)
                ax.set_ylim(ylims)
        else:
            for ax in axs:
                ax.set_xlim(fixed_xlims)
                ax.set_ylim(fixed_ylims)
        
def mean_leaf(leaves):
    """ leaves is a xy list pf dataframes
    """
    s=[curvilinear_abscisse(x['x'],x['y']) for x in leaves]
    curvs=[curvature_xys(x[0]['x'].values,x[0]['y'].values,x[1]) for x in zip(leaves,s)]
    sum_curv = reduce(lambda x,y: [(x[0][0]+y[0][0],x[0][1] + y[0][1]),x[1] + y[1],x[2] + y[2], x[3]+ y[3]],curvs)
    w = 1./len(leaves)
    p,t,ss,dt=sum_curv
    mean_curv = [(p[0]*w, p[1]*w),t*w,ss*w,dt*w]
    x,y = curvature2xy(*mean_curv)
    return pandas.DataFrame({'x':x,'y':y})
    
def plot_mean_leaf_curvature(name = 'Tremie12', nplants = 30, numbering = 'relative_ranktop', 
                              hins_mean = True, fixed_xlims = [-15, 15],
                              fixed_ylims = [-10, 90], set_lims = False, add_simu = True):
    # TODO :
    # Reflechir a une moyenne qui tienne compte de si la feuille est entiere ou non
    # Ou filtrer les outlayers
    df = xydb_reader(name = name)
    fig, axs = plt.subplots(1, len(numpy.unique(df['harvest'])))
    
    plot_leaf_curvature(name = name, nplants = nplants, alternate = True, 
                        numbering = numbering, hins_mean = hins_mean, 
                        alpha = 0.2, axs = axs, set_lims = False, add_sim=add_simu)
    if is_iterable(axs):
        iter_ax = axs.flat
    else:
        iter_ax = iter([axs])
    
    for h in numpy.unique(df['harvest']):
        ax = next(iter_ax)
        ax.set_color_cycle(None)
        colors = ax._get_lines.color_cycle
        df_h = df[df['harvest'] == h]
        df_grouped = df_h.groupby(numbering).mean()
        for leaf in numpy.unique(df_h[numbering]):
            color = next(colors)
            df_lf = df_h[(df_h[numbering] == leaf)]
            side = df_lf.loc[:,'side']
            df_lf.loc[:, 'x'] *= side
            df_mean = mean_leaf([lf.loc[:,['x','y']] for g, lf in df_lf.groupby('inerv')])
            df_mean.loc[:, 'x'] *= numpy.power(-1, leaf%2)
            ax.plot(df_mean.loc[:,'x'], df_mean.loc[:,'y']+df_grouped.loc[leaf,'hins'],
                        color = color, linewidth = 2, linestyle='-')
        if add_simu:
            dfsim = load_simulation_curvature(name)                
            ax.set_color_cycle(None)
            colors = ax._get_lines.color_cycle                
            df_h = dfsim[dfsim['harvest'] == h]
            df_grouped = df_h.groupby(numbering).mean()
            for leaf in numpy.unique(df_h[numbering]):
                color = next(colors)
                df_lf = df_h[(df_h[numbering] == leaf)]
                side = df_lf.loc[:,'side']
                df_lf.loc[:, 'x'] *= side
                df_mean = mean_leaf([lf.loc[:,['x','y']] for g, lf in df_lf.groupby('inerv')])
                df_mean.loc[:, 'x'] *= numpy.power(-1, leaf%2)
                ax.plot(df_mean.loc[:,'x'], df_mean.loc[:,'y']+df_grouped.loc[leaf,'hins'],
                            color = color, linewidth = 2, linestyle='--')                
                        
        ax.annotate('Harvest %d' % h, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
        
    if is_iterable(axs):
        if set_lims == True:
            xlims = [min([ax.get_xlim()[0] for ax in axs]), max([ax.get_xlim()[1] for ax in axs])]
            ylims = [min([ax.get_ylim()[0] for ax in axs]), max([ax.get_ylim()[1] for ax in axs])]
            for ax in axs:
                ax.set_xlim(xlims)
                ax.set_ylim(ylims)
        else:
            for ax in axs:
                ax.set_xlim(fixed_xlims)
                ax.set_ylim(fixed_ylims)
                
    plt.text(0.5, 0.93, name, fontsize=18, 
             transform=fig.transFigure, horizontalalignment='center')