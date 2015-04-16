""" Evaluation of leaf curvature in ADEL"""
from alinea.echap.architectural_data import xydb_reader
from alinea.adel.leaf.curvature import *
import numpy
import pandas
import matplotlib.pyplot as plt
from math import ceil
import collections

# Read curvature data ##############################################################################
def is_iterable(obj):
    return isinstance(obj, collections.Iterable)

def plot_leaf_curvature(name = 'Tremie12', numbering = 'relative_ranktop', 
                        alternate = True, hins_mean = False, 
                        fixed_color = None, alpha = 1., axs = None, 
                        fixed_xlims = [-15, 15], fixed_ylims = [-10, 90], set_lims = True):
    df = xydb_reader(name = name)
    if axs is None:
        fig, axs = plt.subplots(1, len(numpy.unique(df['harvest'])))
    if is_iterable(axs):
        iter_ax = axs.flat
    else:
        iter_ax = iter([axs])
    
    for h in numpy.unique(df['harvest']):
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
        
def plot_mean_by_point_leaf_curvature(name = 'Tremie12', numbering = 'relative_ranktop', 
                                      hins_mean = True, fixed_xlims = [-15, 15],
                                      fixed_ylims = [-10, 90], set_lims = False):
    df = xydb_reader(name = name)
    df.loc[:, 'i_point'] = [float(i)%20 for i in df.index]
    fig, axs = plt.subplots(1, len(numpy.unique(df['harvest'])))
    
    plot_leaf_curvature(name = name, numbering = numbering, 
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
    s=map(lambda x: curvilinear_abscisse(x['x'],x['y']),leaves)
    curvs=map(lambda x:curvature_xys(x[0]['x'].values,x[0]['y'].values,x[1]),zip(leaves,s))
    sum_curv = reduce(lambda x,y: [(x[0][0]+y[0][0],x[0][1] + y[0][1]),x[1] + y[1],x[2] + y[2], x[3]+ y[3]],curvs)
    w = 1./len(leaves)
    p,t,ss,dt=sum_curv
    mean_curv = [(p[0]*w, p[1]*w),t*w,ss*w,dt*w]
    x,y = curvature2xy(*mean_curv)
    return pandas.DataFrame({'x':x,'y':y})
    
def plot_mean_leaf_curvature(name = 'Tremie12', numbering = 'relative_ranktop', 
                              hins_mean = True, fixed_xlims = [-15, 15],
                              fixed_ylims = [-10, 90], set_lims = False):
    df = xydb_reader(name = name)
    fig, axs = plt.subplots(1, len(numpy.unique(df['harvest'])))
    
    plot_leaf_curvature(name = name, alternate = True, 
                        numbering = numbering, hins_mean = hins_mean, 
                        alpha = 0.2, axs = axs, set_lims = False)
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
                
    plt.text(0.5, 0.93, name, fontsize=18, 
             transform=fig.transFigure, horizontalalignment='center')