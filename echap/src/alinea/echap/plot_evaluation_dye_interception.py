import numpy


import alinea.echap.interception_data as idata
import alinea.echap.evaluation_dye_interception as idye

import matplotlib.pyplot as plt
plt.ion()

run_plot=False

'''
# Pour visualiser couvert
%run dye_interception
adel, g, df = repartition_at_application('Tremie13', 'T2')
adel.plot(g,'ntop')
'''


def colors_variety():
    return {'Mercia': 'r', 'Rht3': 'g', 'Tremie12': 'c', 'Tremie13': 'm'}


def xleaf_by():
    return {'metamer': range(7, 14), 'ntop_cur': range(0, 6),
            'ntop_lig': range(-1, 5)}


def prefix_xbar():
    return {'metamer': 'L', 'ntop_cur': 'F', 'ntop_lig': 'Fl'}


def barplot_leaf(ax, obs, sim, loc, xleaf=range(1, 5), prefix='F', o_color='y',
                 s_color='b', opacity=0.4, bar_width=0.4):
    leaves = [prefix + str(x) for x in xleaf]
    xbar = dict(zip(leaves, xleaf))

    obs_bar = None
    if obs is not None:
        if loc in obs.index:
            obs = obs.loc[loc, :]
            x = [xbar.get(leaf, -10) for leaf in obs['xbar']]
            obs_bar = ax.bar(x, obs['ybar'], bar_width, alpha=opacity,
                             color=o_color, yerr=obs['yerr'], ecolor=o_color)

    sim = sim.loc[loc, :]
    x = [xbar.get(leaf, -10) + bar_width for leaf in sim['xbar']]
    sim_bar = ax.bar(x, sim['ybar'], bar_width, alpha=opacity, color=s_color,
                     yerr=sim['yerr'], ecolor=s_color)

    # Mise en forme
    ax.set_xlim(min(xleaf) - bar_width, max(xleaf) + 3 * bar_width)
    ax.set_xticks(numpy.array(xleaf) + bar_width)
    ax.set_xticklabels(leaves, rotation=90, fontsize='small')
    return obs_bar, sim_bar


def fig_observe_simule(obs, sim, treatments=['T1', 'T2'], xleaf=range(1, 5),
                       prefix='F', s_color='b', ylim=None, ylab=None,
                       title=None, add_obs=None, add_sim=None):
    nt = len(treatments)
    lg = max(1, nt / 2)
    fig, axes = plt.subplots(nrows=lg, ncols=nt / lg + (nt - nt / lg * lg),
                             sharey=True)
    axlist = fig.get_axes()
    for ifig, treatment in enumerate(treatments):
        ax = axlist[ifig]
        obs_bar, sim_bar = barplot_leaf(ax, obs, sim, treatment, xleaf=xleaf,
                                        prefix=prefix, s_color=s_color)
        if add_sim is not None:
            obs_bar, sim_bar = barplot_leaf(ax, add_obs, add_sim, treatment,
                                            xleaf=xleaf, prefix=prefix,
                                            s_color=s_color)

        if ifig == 0 and ylab is not None:
            ax.set_ylabel(ylab)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.text(min(xleaf) + 0.5, .8 * max(ax.get_ylim()), '' + str(treatment),
                bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.15)
    if title is not None:
        fig.text(0.5, 0.05, title, size='large', ha='center')
    return fig


def deposit_observe_simule(variety='Tremie12', nplants=30, nrep=1, axis='MS',
                           by='ntop_cur', xleaf=None, ylim=(0, 6.5),
                           simulation='reference', treatments=None, reset=False,
                           reset_data=False, frac=1, model='adel'):
    if treatments is None:
        treatments = idata.tag_treatments()[variety]['application']
    # obs
    df_obs = idata.dye_interception()
    if by in df_obs.columns:
        obs = idye.leaf_statistics(df_obs.loc[variety, :], what='deposit',
                                   by=by, axis=axis)
    else:
        obs = None
    # sim
    df_sim = idye.dye_interception(variety, nplants=nplants, nrep=nrep,
                                   simulation=simulation, treatments=treatments,
                                   reset=reset, reset_reconstruction=reset_data)
    df_sim['deposit_Tartrazine'] *= frac
    if model == 'adel':
        sim = idye.leaf_statistics(df_sim, what='deposit_Tartrazine', by=by,
                               axis=axis)
    elif model == 'miller':
        sim = idye.dye_interception_miller(variety=variety, nplants=nplants,
                                           tag=simulation, rep=1, at=treatments)
    else:
        raise ValueError('unknown model: ' + model)
    # plot
    color = colors_variety()[variety]
    prefix = prefix_xbar()[by]
    if xleaf is None:
        xleaf = xleaf_by()[by]
    ylab = 'Deposit (g per g.cm-2)'
    title = variety + ' ' + simulation
    fig = fig_observe_simule(obs, sim, treatments=treatments, xleaf=xleaf,
                             prefix=prefix, s_color=color, ylim=ylim, ylab=ylab,
                             title=title)
    return fig


if run_plot:
    fig = deposit_observe_simule('Tremie12', nplants=30)
    fig.savefig('deposit_Tremie12.png')


def hs_observe_simule(variety='Tremie12', nplants=30, nrep=1, axis='MS',
                      simulation='reference', treatments=None, reset=False,
                      reset_data=False):
    if treatments is None:
        treatments = idata.tag_treatments()[variety]['hs']
    # obs
    df_obs = idata.haun_stages()
    df_obs = df_obs.loc[variety, :]
    obs = df_obs.set_index('treatment')
    # sim
    df_sim = idye.dye_interception(variety, nplants=nplants, nrep=nrep,
                                   simulation=simulation, treatments=treatments,
                                   reset=reset, reset_data=reset_data)
    df_sim['HStarget'] = numpy.minimum(df_sim['HS'], df_sim['nff'])
    sim = idye.axis_statistics(df_sim, what='haun_stage', axis=axis)

    bar_width = 0.2;
    opacity = 0.4
    colors = {'Tremie12': 'c', 'Tremie13': 'm'}
    fig, ax = plt.subplots(nrows=1, ncols=1)
    xbar = 0.
    xstage = []
    sim_bars = {}
    var = variety
    for t in treatments:
        # target_bar = ax.bar(xbar, target, bar_width, alpha=opacity, color='orange')
        # xbar += bar_width
        if t in obs.index:
            obs_bar = ax.bar(xbar, obs.ix[t, 'HS_mean'], bar_width,
                             alpha=opacity, color='y',
                             yerr=obs.ix[t, 'HS_conf'], ecolor='y')
        xbar += bar_width
        xstage.append(xbar)
        sim_bars = ax.bar(xbar, sim.ix[t, 'ybar'], bar_width, alpha=opacity,
                          color=colors[var], yerr=sim.ix[t, 'yerr'],
                          ecolor=colors[var])
        xbar += bar_width * 1.5
    ax.set_ylim(0, 14)
    ax.set_xlim(-bar_width, xbar)
    ax.set_xticks(xstage)
    ax.set_xticklabels(treatments, fontsize='small')
    ax.text(0.1, 13.1, '' + str(var),
            bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10}, fontsize=12)
    fig.text(0.5, 0.05, 'HS', size='large', ha='center')
    return fig, ax


def scan_observe_simule(variety='Tremie12', nplants=30, nrep=1, axis='MS',
                        by='metamer', xleaf=None, ylim=(0, 30),
                        simulation='reference', treatments=None, reset=False,
                        reset_data=False):
    if treatments is None:
        treatments = idata.tag_treatments()[variety]['scan']
    # obs
    df_obs = idata.scan_data()
    if variety in df_obs.index:
        obs = idye.leaf_statistics(df_obs.loc[variety, :], what='A_bl', by=by,
                                   axis=axis)
        add_obs = idye.leaf_statistics(df_obs.loc[variety, :],
                                       what='A_bl_green', by=by, axis=axis)
    else:
        obs = None
        add_obs = None
    # sim
    df_sim = idye.dye_interception(variety, nplants=nplants, nrep=nrep,
                                   simulation='reference',
                                   treatments=treatments)
    sim = idye.leaf_statistics(df_sim, what='area', by=by, axis=axis)
    add_sim = idye.leaf_statistics(df_sim, what='green_area', by=by, axis=axis)
    # plot
    color = colors_variety()[variety]
    prefix = prefix_xbar()[by]
    if xleaf is None:
        xleaf = xleaf_by()[by]
    ylab = 'Leaf Area (cm2)'
    title = variety + ' ' + simulation
    fig = fig_observe_simule(obs, sim, treatments=treatments, xleaf=xleaf,
                             prefix=prefix, s_color=color, ylim=ylim, ylab=ylab,
                             title=title, add_obs=add_obs, add_sim=add_sim)
    return fig


def sil_observe_simule(variety='Tremie12', nplants=30, nrep=1, axis='MS',
                       by='ntop_cur', xleaf=None, ylim=(0, 1),
                       what='h_projection', simulation='reference',
                       treatments=None, reset=False, reset_data=False):
    if treatments is None:
        treatments = idata.tag_treatments()[variety]['silhouette']
    df_obs = idata.silhouettes()
    if variety in df_obs.index:
        obs = idye.leaf_statistics(df_obs.loc[variety, :], what=what, by=by,
                                   axis=axis)
    else:
        obs = None
    # sim
    df_sim = idye.dye_interception(variety, nplants=nplants, nrep=nrep,
                                   simulation=simulation, treatments=treatments,
                                   reset=reset, reset_data=reset_data)
    sim = idye.leaf_statistics(df_sim, what=what, by=by, axis=axis)
    # plot
    color = colors_variety()[variety]
    prefix = prefix_xbar()[by]
    if xleaf is None:
        xleaf = xleaf_by()[by]
    ylab = 'Projection factor'
    title = variety + ' ' + simulation
    fig = fig_observe_simule(obs, sim, treatments=treatments, xleaf=xleaf,
                             prefix=prefix, s_color=color, ylim=ylim, ylab=ylab,
                             title=title)
    return fig


def petri_observe_simule(variety='Tremie12', nplants=30, nrep=1, axis='MS',
                         simulation='reference', treatments=None, reset=False,
                         reset_data=False):
    if treatments is None:
        treatments = idata.tag_treatments()[variety]['application']
    # obs
    df_obs = idata.petri_dye_interception()
    df_obs = df_obs.loc[variety, :]
    obs = df_obs.set_index('treatment')
    # sim
    df_sim = idye.dye_interception(variety, nplants=nplants, nrep=nrep,
                                   simulation=simulation, treatments=treatments,
                                   reset=reset, reset_data=reset_data)
    df_sim['HStarget'] = numpy.minimum(df_sim['HS'], df_sim['nff'])
    sim = idye.axis_statistics(df_sim, what='haun_stage', axis=axis)

    bar_width = 0.2;
    opacity = 0.4
    colors = {'Tremie12': 'c', 'Tremie13': 'm'}
    fig, ax = plt.subplots(nrows=1, ncols=1)
    xbar = 0.
    xstage = []
    sim_bars = {}
    var = variety
    for t in treatments:
        # target_bar = ax.bar(xbar, target, bar_width, alpha=opacity, color='orange')
        # xbar += bar_width
        if t in obs.index:
            obs_bar = ax.bar(xbar, obs.ix[t, 'HS_mean'], bar_width,
                             alpha=opacity, color='y',
                             yerr=obs.ix[t, 'HS_conf'], ecolor='y')
        xbar += bar_width
        xstage.append(xbar)
        sim_bars = ax.bar(xbar, sim.ix[t, 'ybar'], bar_width, alpha=opacity,
                          color=colors[var], yerr=sim.ix[t, 'yerr'],
                          ecolor=colors[var])
        xbar += bar_width * 1.5
    ax.set_ylim(0, 14)
    ax.set_xlim(-bar_width, xbar)
    ax.set_xticks(xstage)
    ax.set_xticklabels(treatments, fontsize='small')
    ax.text(0.1, 13.1, '' + str(var),
            bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10}, fontsize=12)
    fig.text(0.5, 0.05, 'HS', size='large', ha='center')
    return fig, ax


